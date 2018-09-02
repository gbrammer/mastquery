"""
Utilities
"""
import os
import warnings
import numpy as np

from astropy.table import Table

try: # Python 3.x
    from urllib.parse import quote as urlencode
    from urllib.request import urlretrieve
except ImportError:  # Python 2.x
    from urllib import pathname2url as urlencode
    from urllib import urlretrieve

try: # Python 3.x
    import http.client as httplib 
except ImportError:  # Python 2.x
    import httplib

################
## mastQuery and mastJson2Table from 
## https://mast.stsci.edu/api/v0/MastApiTutorial.html
def mastQuery(request):

    server='mast.stsci.edu'

    # Grab Python Version 
    version = ".".join(map(str, sys.version_info[:3]))

    # Create Http Header Variables
    headers = {"Content-type": "application/x-www-form-urlencoded",
               "Accept": "text/plain",
               "User-agent":"python-requests/"+version}

    # Encoding the request as a json string
    requestString = json.dumps(request)
    requestString = urlencode(requestString)
    
    # opening the https connection
    conn = httplib.HTTPSConnection(server)

    # Making the query
    conn.request("POST", "/api/v0/invoke", "request="+requestString, headers)

    # Getting the response
    resp = conn.getresponse()
    head = resp.getheaders()
    content = resp.read().decode('utf-8')

    # Close the https connection
    conn.close()

    return head,content

def mastJson2Table(jsonObj):

    dataTable = Table()

    for col,atype in [(x['name'],x['type']) for x in jsonObj['fields']]:
        if atype=="string":
            atype="str"
        if atype=="boolean":
            atype="bool"
        dataTable[col] = np.array([x.get(col,None) for x in jsonObj['data']],dtype=atype)
        
    return dataTable
###############

def table_from_info(info):
    """
    Generate a query-like table based on header keywords parsed by 
    `~grizli.pipeline.auto_script.parse_visits`.
    
    """
    from astropy.table import Table
    import astropy.wcs as pywcs
    import astropy.io.fits as pyfits
    
    from . import query
    
    tab = Table()
    tab['visit_duration'] = info['EXPTIME']
    tab['ra'] = info['RA_TARG']
    tab['dec'] = info['DEC_TARG']
    tab['exptime'] = info['EXPTIME']
    tab['detector'] = info['DETECTOR']
    
    swap_detector = {}
    for k in query.INSTRUMENT_DETECTORS:
        swap_detector[query.INSTRUMENT_DETECTORS[k]] = k
    
    tab['instdet'] = [swap_detector[det] for det in info['DETECTOR']]
    tab['aperture'] = info['DETECTOR']
    
    tab['observation_id'] = [os.path.basename(file).split('_')[0].lower() for file in info['FILE']]
    tab['file_type'] = [os.path.basename(file).split('_')[1].upper() for file in info['FILE']]

    tab['artifact_id'] = info['FILE']
    tab['target'] = info['TARGNAME']
    tab['filter'] = info['FILTER']
    
    # Footprints
    footprints = []
    proposal_id = []
    for file in info['FILE']:
        im = pyfits.open(file)
        wcs = pywcs.WCS(im['SCI',1].header, fobj=im)
        fp = wcs.calc_footprint()
        footstr = '{'+', '.join([list(fpi).__repr__()[1:-1] for fpi in fp])+'}'
        footprints.append(footstr)
        proposal_id.append(im[0].header['PROPOSID'])
    
    tab['proposal_id'] = proposal_id
    tab['footprint'] = footprints
    tab['stc_s_tailored'] = tab['footprint']
    
    return tab
    
def set_warnings(numpy_level='ignore', astropy_level='ignore'):
    """
    Set global numpy and astropy warnings
    
    Parameters
    ----------
    numpy_level : {'ignore', 'warn', 'raise', 'call', 'print', 'log'}
        Numpy error level (see `~numpy.seterr`).
        
    astropy_level : {'error', 'ignore', 'always', 'default', 'module', 'once'}
        Astropy error level (see `~warnings.simplefilter`).
    
    """
    from astropy.utils.exceptions import AstropyWarning
    
    np.seterr(all=numpy_level)
    warnings.simplefilter(astropy_level, category=AstropyWarning)

def radec_to_targname(ra=0, dec=0, scl=10000):
    """Turn decimal degree coordinates into a string
    
    Example:

        >>> from grizli.utils import radec_to_targname
        >>> print(radec_to_targname(ra=10., dec=-10.))
        j004000-100000
    
    Parameters
    -----------
    ra, dec : float
        Sky coordinates in decimal degrees
    
    Returns
    --------
    targname : str
        Target name like jHHMMSS[+-]DDMMSS.
    
    """
    import astropy.coordinates 
    import astropy.units as u
    import re
    import numpy as np
    
    dec_scl = int(np.round(dec*scl))/scl
    cosd = np.cos(dec_scl/180*np.pi)
    ra_scl = int(np.round(ra*scl/cosd))/(scl/cosd)
    
    coo = astropy.coordinates.SkyCoord(ra=ra_scl*u.deg, dec=dec_scl*u.deg)
    
    cstr = re.split('[hmsd.]', coo.to_string('hmsdms', precision=2))
    targname = ('j{0}{1}'.format(''.join(cstr[0:3]), ''.join(cstr[4:7])))
    targname = targname.replace(' ', '')
    
    return targname
    
def get_irsa_dust(ra, dec, type='SandF'):
    """
    Get Galactic dust reddening from NED/IRSA at a given position
    http://irsa.ipac.caltech.edu/applications/DUST/docs/dustProgramInterface.html
    
    Parameters
    ----------
    ra, dec : float
        RA/Dec in decimal degrees.
        
    type : 'SFD' or 'SandF'
        Dust model, with        
            SandF = Schlafly & Finkbeiner 2011 (ApJ 737, 103) 
              SFD = Schlegel et al. 1998 (ApJ 500, 525)
    
    Returns
    -------
    ebv : float
        Color excess E(B-V), in magnitudes
    
    """
    import os
    import tempfile   
    import urllib.request
    from astropy.table import Table
    from lxml import objectify
    
    query = 'http://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr={0:.4f}+{1:.4f}+equ+j2000'.format(ra, dec)
    
    req = urllib.request.Request(query)
    response = urllib.request.urlopen(req)
    resp_text = response.read().decode('utf-8')
    
    root = objectify.fromstring(resp_text)
    stats = root.result.statistics

    if type == 'SFD':
        return float(str(stats.refPixelValueSFD).split()[0])
    else:
        return float(str(stats.refPixelValueSandF).split()[0])
        