"""
Utilities
"""
import os
import json

import warnings
import numpy as np

from astropy.table import Table

# character to skip clearing line on STDOUT printing
NO_NEWLINE = '\x1b[1A\x1b[1M'

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
    import sys
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

NIRCam_LW_APERTURES = ['NIRCam/'+a for a in ['NRCB5_FULL', 'NRCA5_FULL']]
NIRCam_SW_APERTURES = ['NIRCam/'+a for a in ['NRCA1_FULL', 'NRCA2_FULL', 'NRCA3_FULL', 'NRCA4_FULL', 'NRCB1_FULL', 'NRCB2_FULL', 'NRCB3_FULL', 'NRCB4_FULL']]
MIRI_APERTURES = ['MIRI/MIRIM_ILLUM']
NIRISS_APERTURES = ['NIRISS/NIS_CEN']
NIRSpec_APERTURES = ['NIRSpec/'+a for a in ['NRS_FULL_MSA1', 'NRS_FULL_MSA2','NRS_FULL_MSA3','NRS_FULL_MSA4', 'NRS_FULL_IFU']]

JWST_APERTURES = NIRCam_LW_APERTURES + NIRCam_SW_APERTURES + NIRSpec_APERTURES + NIRISS_APERTURES + MIRI_APERTURES

def get_jwst_colors():
    import matplotlib.pyplot as plt
    # colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    # for i in enumerate(['NIRCam','NIRSpec','NIRISS','MIRI']):
    #     ap_colors[ap] = colors[i]
    
    ap_colors = {'NIRCam':'#17becf','NIRSpec':'#bcbd22', 'NIRISS':'#1f77b4','MIRI':'#d62728', 'NIRCam-SW':'#e377c2', 'NIRCam-LW':'#9467bd'}
    
    return ap_colors
    
def get_jwst_apertures(ra=0, dec=0, pa_v3=0, ref_aper='NIRISS/NIS_CEN', aper_list=JWST_APERTURES, patch_kwargs={'alpha':0.5}, regions_file=None, ax=None):
    """
    
    Get JWST aperture footprints using `pysiaf`.
    
    """
    import matplotlib.pyplot as plt
    
    from collections import OrderedDict
    import pysiaf
    
    ap_colors = get_jwst_colors()
    
    instruments = np.unique([ap.split('/')[0] for ap in aper_list+[ref_aper]])
    siaf = {}
    for inst in instruments:
        siaf[inst] = pysiaf.Siaf(inst)
    
    spl = ref_aper.split('/')
    ref_ap = siaf[spl[0]][spl[1]]
    
    theta = pa_v3/180*np.pi
    _mat = np.array([[np.cos(theta), -np.sin(theta)],
                     [np.sin(theta), np.cos(theta)]])
    
    aps = OrderedDict()
    cosd = np.cos(dec/180*np.pi)
    
    for ap in aper_list:
        spl = ap.split('/')
        aper = siaf[spl[0]][spl[1]]
        xy = np.array(aper.corners('tel')).T - ref_ap.reference_point('tel')
        
        # Close
        xy = xy[list(range(4))+[0],:]
        
        xyr = np.dot(xy, _mat)
        xydeg = xyr/3600.
        xydeg[:,0] /= cosd
        xysky = xydeg + np.array([ra, dec])
        
        if spl[0] == 'NIRCam':
            if ap in NIRCam_LW_APERTURES:
                inst = 'NIRCam-LW'
            else:
                inst = 'NIRCam-SW'
        else:
            inst = spl[0]
            
        aps[ap] = {'instrument':inst, 'aperture':spl[1], 
                   'footprint':xysky.T, 'tel':xy.T, 'tel-rot':xyr.T, 
                   'reference_point':ref_ap.reference_point('tel'),
                   'patch':plt.Polygon(xysky, fc=ap_colors[inst],
                                       ec=ap_colors[inst], **patch_kwargs)}
        
        pstr =  ', '.join(['{0:.5f}'.format(c) for c in xysky.flatten()])
        aps[ap]['reg'] = 'polygon({0}) # '.format(pstr) + 'color={'+  ap_colors[inst] + '}'
    
    if regions_file is not None:
        fp = open(regions_file,'w')
        for ap in aps:
            fp.write(aps[ap]['reg']+'\n')
        fp.close()
    
    if ax is not None:
        for ap in aps:
            ax.add_patch(aps[ap]['patch'])
            
    return aps
    
        
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

def polygon_to_sregion(poly):
    """
    Convert `shapely.Polygon` vertices to an SREGION string
    """
    try:
        xy = np.array(poly.boundary.xy).T.flatten()
    except:
        xy = np.array(poly.convex_hull.boundary.xy).T.flatten()
    
    pstr = 'POLYGON({0})'.format(','.join(['{0:.6f}'.format(c) for c in xy]))
    return pstr

def sregion_to_polygon(pstr):
    """
    Convert `shapely.Polygon` vertices to an SREGION string
    """
    from shapely.geometry import Polygon
    coo = np.cast[float](pstr.lower().strip('polygon(').strip(')').split(','))
    poly = Polygon(coo.reshape(-1,2))
    return poly
    
def radec_to_targname(ra=0, dec=0, round_arcsec=(4, 60), precision=2, targstr='j{rah}{ram}{ras}{sign}{ded}{dem}', header=None):
    """Turn decimal degree coordinates into a string with rounding.
    
    Example:
        
        # Test dec: -10d10m10.10s
        >>> dec = -10. - 10./60. - 10.1/3600
        
        # Test ra: 02h02m02.20s
        >>> cosd = np.cos(dec/180*np.pi)
        >>> ra = 2*15 + 2./60*15 + 2.2/3600.*15
        
        # Round to nearest arcmin
        >>> from mastquery.utils import radec_to_targname
        >>> print(radec_to_targname(ra=ra, dec=dec, round_arcsec=(4,60),
                               targstr='j{rah}{ram}{ras}{sign}{ded}{dem}'))
        j020204m1010 # (rounded to 4 arcsec in RA)
        
        # Full precision
        >>> targstr = 'j{rah}{ram}{ras}.{rass}{sign}{ded}{dem}{des}.{dess}'
        >>> print(radec_to_targname(ra, dec,round_arcsec=(0.0001, 0.0001),
                                    precision=3, targstr=targstr))
        j020202.200m101010.100
        
    Parameters
    -----------
    ra, dec : float
        Sky coordinates in decimal degrees
    
    round_arcsec : (scalar, scalar) 
        Round the coordinates to nearest value of `round`, in arcseconds.
    
    precision : int
        Sub-arcsecond precision, in `~astropy.coordinates.SkyCoord.to_string`.
        
    targstr : string
        Build `targname` with this parent string.  Arguments 
        `rah, ram, ras, rass, sign, ded, dem, des, dess` are computed from the 
        (rounded) target coordinates (`ra`, `dec`) and passed to 
        `targstr.format`.
    
    header : `~astropy.io.fits.Header`, None
        Try to get `ra`, `dec` from header keywords, first `CRVAL` and then
        `RA_TARG`, `DEC_TARG`.
        
    Returns
    --------
    targname : str
        Target string, see the example above.
    
    """
    import astropy.coordinates 
    import astropy.units as u
    
    import re
    import numpy as np
    
    if header is not None:
        if 'CRVAL1' in header:
            ra, dec = header['CRVAL1'], header['CRVAL2']
        else:
            if 'RA_TARG' in header:
                ra, dec = header['RA_TARG'], header['DEC_TARG']
    
    cosd = np.cos(dec/180*np.pi)
    scl = np.array(round_arcsec)/3600*np.array([360/24, 1])
    
    dec_scl = int(np.round(dec/scl[1]))*scl[1]
    ra_scl = int(np.round(ra/scl[0]))*scl[0]
    
    coo = astropy.coordinates.SkyCoord(ra=ra_scl*u.deg, dec=dec_scl*u.deg)
    
    cstr = re.split('[hmsd.]', coo.to_string('hmsdms', precision=precision))
    # targname = ('j{0}{1}'.format(''.join(cstr[0:3]), ''.join(cstr[4:7])))
    # targname = targname.replace(' ', '').replace('+','p').replace('-','m')

    rah, ram, ras, rass = cstr[0:4]
    ded, dem, des, dess = cstr[4:8]
    sign = 'p' if ded[1] == '+' else 'm'
    
    targname = targstr.format(rah=rah, ram=ram, ras=ras, rass=rass,
                              ded=ded[2:], dem=dem, des=des, dess=dess,
                              sign=sign)
        
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
        