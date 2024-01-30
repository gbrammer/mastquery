"""
Utilities
"""
import os
import inspect
import json

import warnings
import logging

import numpy as np

from astropy.table import Table

LOGFILE = 'mastquery.log'

INSTRUMENT_AREAS = {'WFC3/IR':4.5,
                    'WFPC2/WFC':8.,
                    'WFPC2/PC':8.,
                    'ACS/WFC':11.3,
                    'WFC3/UVIS':7.3,
                    'NIRSPEC/IFU':0.00318, 
                    'NIRSPEC/MSA':16,
                    'NIRSPEC/SLIT':0.0002, 
                    'NIRCAM/IMAGE':12.9,
                    'NIRCAM/GRISM':11.8,
                    'NIRISS/IMAGE':5.8,
                    'MIRI/IMAGE':2.8,
                    'MIRI/SLIT':0.0007,
                    'MIRI/IFU':0.0037, 
                    'MIRI':2.8*0.001, # Todo - parse mode e.g IFU
                    'NIRISS':5.8*0.001, 
                    'NIRCAM':4.8/2.*0.001,
                    'NIRSPEC':12*0.00000001}

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
        dataTable[col] = np.array([x.get(col,None) for x in jsonObj['data']])
        
    return dataTable


def new_mast_query(request):
    """Perform a MAST query
    
    From https://mast.stsci.edu/api/v0/MastApiTutorial.html
    
    Parameters
    ----------
    request : dict
        The MAST request json object

    Returns
    -------
    head : str
        Response HTTP header
    
    content : str
        Returned data
    """
    import sys
    import json
    import requests
    from urllib.parse import quote as urlencode
    
    # Base API url
    request_url='https://mast.stsci.edu/api/v0/invoke'    
    
    # Grab Python Version 
    version = ".".join(map(str, sys.version_info[:3]))

    # Create Http Header Variables
    headers = {"Content-type": "application/x-www-form-urlencoded",
               "Accept": "text/plain",
               "User-agent":"python-requests/"+version}

    # Encoding the request as a json string
    req_string = json.dumps(request)
    req_string = urlencode(req_string)
    
    # Perform the HTTP request
    resp = requests.post(request_url, data="request="+req_string, headers=headers)
    
    # Pull out the headers and response content
    head = resp.headers
    content = resp.content.decode('utf-8')

    return head, content


def new_mastJson2Table(query_content):
    """
    Convert json query content to table
    """
    import json
    from astropy.table import Table
    
    json_data = json.loads(query_content)
    #print('xxx', json_data.keys())
    
    tabs = []
    for jdata in json_data['Tables']:
        cols = [c['name'] for c in jdata['Fields']]
        tab = Table(names=cols, rows=jdata['Rows'])
        tabs.append(tab)
    
    if len(tabs) == 1:
        return tabs[0]
    else:
        return tabs


def split_rateints(file, verbose=True, overwrite=False, sigma_clip=4):
    """
    Combine integrations in a ``rateints`` file with sigma clipping
    """
    import astropy.io.fits as pyfits
    from astropy.table import Table

    ints = pyfits.open(file)
    
    if 'INT_TIMES' not in ints:
        msg = f'mastquery.utils.split_rateints: INT_TIMES not found in {file}, skip'
        log_comment(LOGFILE, msg, verbose=verbose)
        return None
        
    N = ints[0].header['NINTS']
    
    msg = f'mastquery.utils.split_rateints: {file} NINTS={N} sigma_clip={sigma_clip}'
    log_comment(LOGFILE, msg, verbose=verbose)
    
    out = file.replace('_rateints.fits', '_rate.fits')
    
    # Array copies
    sci = ints['SCI'].data*1
    dq = ints['DQ'].data*1
    err = ints['ERR'].data*1
    var_poisson = ints['VAR_POISSON'].data*1
    var_rnoise = ints['var_rnoise'].data*1
    
    ivar = 1/err**2
    
    mask = (dq & (1+1024+4096)) == 0
    sci[~mask] = np.nan
    ivar[~mask] = 0.
    
    ### sigma-clipping
    med = np.nanmedian(sci, axis=0)
    clip = np.abs(sci - med)*np.sqrt(ivar) > sigma_clip
    ivar[clip] = 0
    sci[~mask | clip] = 0
    
    ### Weighted combination
    num = np.nansum(sci*ivar, axis=0)
    den = np.nansum(ivar, axis=0)
    avg = num/den
    err = 1/np.sqrt(den)
    
    ### Update variances
    var_poisson[~mask | clip] = np.nan
    var_rnoise[~mask | clip] = np.nan
    var_poisson = 1/np.nansum(1./var_poisson, axis=0)
    var_rnoise = 1/np.nansum(1./var_rnoise, axis=0)
    
    ### Combined DQ 
    dqe = dq[0,:,:]*(~mask | clip)[0,:,:]
    for i in range(N):
        dqe |= dq[i,:,:]*(~mask | clip)[i,:,:]
    
    ### valid data
    valid = np.isfinite(avg + err + var_poisson + var_rnoise)
    avg[~valid] = 0
    err[~valid] = 0
    var_poisson[~valid] = 0
    var_rnoise[~valid] = 0
    dqe[~valid] |= 1
    
    ### Update FITS arrays and header
    ints['SCI'].data = avg
    ints['ERR'].data = err
    ints['DQ'].data = dqe
    ints['VAR_POISSON'].data = var_poisson
    ints['VAR_RNOISE'].data = var_rnoise
    
    ints[0].header['DATAMODL'] = 'ImageModel'
    ints[0].header['FILENAME'] = out
    
    ints.writeto(out, overwrite=True)
    ints.close()
    
    return [out]


def download_recalibrated_rate(rate_file, bucket="s3://grizli-v2/reprocess_rate/", path=None, verbose=True, overwrite=True, **kwargs):
    """Get recalibrated rate file if available.
    
    These will have been processed with `grizli.aws.recalibrate` and `snowblind`.
    """
    import logging

    try:
      import boto3
    except ImportError: 
      return None
    
    from botocore.exceptions import ClientError
    
    log = logging.getLogger()

    if path is not None:
        if not os.path.exists(path):
            log.info(f'mkdir {path}')
            os.makedirs(path)
            local_file = os.path.join(path, rate_file)
    else:
        path = './'
        local_file = rate_file
    
    if (not overwrite) & os.path.exists(local_file):
        msg = f'download_recalibrated_rate: {local_file} file found'
        _ = log_comment(LOGFILE, msg, verbose=verbose)
        return local_file
    
    # Check if requested file is in the s3 bucket
    index_file = os.path.join(bucket.replace('s3://', 'https://s3.amazonaws.com/'),
                              'index.csv')
    try:
        filelist = Table.read(index_file, format='csv')
    except FileNotFoundError:
        msg = f'download_recalibrated_rate: index {index_file} not found'
        _ = log_comment(LOGFILE, msg, verbose=verbose)
        return None
    
    if rate_file not in filelist['file']:
        msg = f'download_recalibrated_rate: {rate_file} not in {index_file}'
        _ = log_comment(LOGFILE, msg, verbose=verbose)
        return None
    
    # Try to download it
    s3 = boto3.resource('s3')
    
    bkt_name = bucket.replace('s3://','').split('/')[0]
    bkt_prefix = '/'.join(bucket.replace('s3://','').split('/')[1:])
    
    bkt = s3.Bucket(bkt_name)
    
    s3_prefix = os.path.join(bkt_prefix, rate_file)
    
    msg = f"download_recalibrated_rate: download s3://{bkt_name}/{s3_prefix} > {local_file}"
    log_comment(LOGFILE, msg, verbose=verbose)
                
    try:
        bkt.download_file(s3_prefix, local_file, Callback=None, Config=None)
    except ClientError:
        msg = f"download_recalibrated_rate: download failed with ClientError"
        log_comment(LOGFILE, msg, verbose=verbose)
        return None
    
    if not os.path.exists(local_file):
        msg = f'download_recalibrated_rate: boto3 download but {local_file} not found'
        log_comment(LOGFILE, msg, verbose=verbose)
        return None
    else:
        return local_file


def download_from_mast(tab, path=None, verbose=True, overwrite=False, use_token=True, base_url=None, cloud_only=False, force_rate=False, rate_ints=False, get_recalibrated_rate=False, recalibrated_kwargs={}, **kwargs):
    """
    Download files from MAST API with `astroquery.mast.Observations.download_file`
    
    Parameters
    ----------
    tab : table, list
        Table from MAST API query with minimal column ``dataURI``
        or ``dataURL``.  If a `list`, then assumes they are strings of 
        ``dataURL``/``dataURI`` names.
    
    path : str
        Output path, defaults to current working diretory
    
    verbose : bool
        Print status messages [deprecated, now uses `logging`]
    
    overwrite : bool
        Overwrite existing files
    
    use_token : bool
        Try to use a MAST token defined in a ``MAST_TOKEN`` environment 
        variable.  See https://auth.mast.stsci.edu/info for info about the 
        token authentication.
    
    base_url : str
        A base url to use when downloading.  Default is the MAST Portal API
    
    cloud_only : bool
        See `astroquery.mast.Observations.download_file`
    
    force_rate : bool
        Replace 'cal' with 'rate' in filenames.
    
    rate_ints : bool
        Fetch ``rateints`` files and split into individual ``rate`` files
    
    Returns
    -------
    resp : dict
        Dict of responses from `astroquery.mast.Observations.download_file`
        
    """
    from astroquery.mast import Observations
    log = logging.getLogger()
    
    if (os.getenv('MAST_TOKEN') is None) | (not use_token):
        session = Observations
    else:
        try:
            session = Observations(mast_token=os.getenv('MAST_TOKEN'))
        except:
            session = Observations(token=os.getenv('MAST_TOKEN'))
    
    # Input is a list
    if isinstance(tab, list):
        tab = Table(names=['dataURI'], rows=[[u] for u in tab])
    
    if path is not None:
        if not os.path.exists(path):
            log.info(f'mkdir {path}')
            os.makedirs(path)
            
    kws = {'local_path': path, 
           'cache':(not overwrite), 
           'base_url':base_url,
           'cloud_only':cloud_only}
    
    resp = {}
    for row in tab:     
        if 'dataURI' in row.colnames:
            _uri = row['dataURI']
        else:
            _uri = row['dataURL']

        _file = os.path.basename(_uri)
        rate_file = _file.replace('_cal','_rate')
        rate_file = rate_file.replace('_uncal','_rate')
        
        if force_rate:
            _file = rate_file
            _uri = _uri.replace('_cal','_rate')
        
        if rate_ints:
            _file = rate_file.replace('_rate.fits','_rateints.fits')
            _uri = _uri.replace('_cal','_rate')
            _uri = _uri.replace('_rate.fits','_rateints.fits')
            
        if path is None:
            out_file = _file
        else:
            out_file = os.path.join(path, os.path.basename(_file))
            kws['local_path'] = out_file
        
        if get_recalibrated_rate & rate_file.endswith('_rate.fits'):
            status = download_recalibrated_rate(rate_file,
                                                path=path,
                                                verbose=verbose,
                                                overwrite=overwrite)
            if (status is not None):
                resp[out_file] = ('EXISTS', None, None)
                continue
        
        if os.path.exists(out_file) & (not overwrite):
            log.info(f'{out_file} exists, skip')
            resp[out_file] = ('EXISTS', None, None)
            
            if rate_ints == 1:
                out_rates = split_rateints(out_file, verbose=True)
                for f in out_rates:
                    resp[f] = ('EXISTS', None, None)
            
                _ = resp.pop(out_file)
            
            continue
                
        # Download the data
        if _uri.startswith('jw'):
            _uri = 'mast:JWST/product/'+_uri
            
        log.info(f"Download: {out_file}")
             
        payload = {"uri":_uri}        
        resp[out_file] = session.download_file(_uri, **kws)
        
        if (rate_ints == 1) & os.path.exists(out_file):
            out_rates = split_rateints(out_file, verbose=True)
            for f in out_rates:
                resp[f] = ('COMPLETE', None, None)
            
            _ = resp.pop(out_file)
            
    return resp


def old_download_from_mast(tab, path='./', verbose=True, overwrite=True, min_size=1, delete_small=True):
    """
    Download files from MAST API
    
    Parameters
    ----------
    tab : table, list
        Table from MAST API query with minimal column ``dataURI``
        or ``dataURL``.  If a `list`, then assumes they are strings of 
        ``dataURL``/``dataURI`` names.
    
    path : str
        Output path
    
    verbose : bool
        Print status messages [deprecated, now uses `logging`]
    
    overwrite : bool
        Overwrite existing files
    
    min_size : float
        Minimum size in megabits to check if the downloaded file is rather
        an `Access Denied` file
      
    delete_small : bool
        Remove the file if it failed the `min_size` test
    
    Returns
    -------
    outlist : list
        List of valid downloaded files
        
    """
    import requests
    import os
    log = logging.getLogger()
    
    download_url = 'https://mast.stsci.edu/api/v0.1/Download/file?'
    
    # Input is a list
    if isinstance(tab, list):
        tab = Table(names=['dataURI'], rows=[[u] for u in tab])
        
    if not os.path.exists(path):
        log.info(f'mkdir {path}')
        os.makedirs(path)
    
    outlist = []
    
    for row in tab:     
        if 'dataURI' in row.colnames:
            _uri = row['dataURI']
        else:
            _uri = row['dataURL']

        if 'filename' in row.colnames:
            _file = row['filename']
        else:
            _file = os.path.basename(_uri)
            
        out_file = os.path.join(path, os.path.basename(_file))
        
        if os.path.exists(out_file) & (not overwrite):
            log.info(f'{out_file} exists, skip')
            outlist.append(out_file)
            continue
                
        # Download the data   
        log.info(f"Download: {out_file}")
             
        payload = {"uri":_uri}        
        resp = requests.get(download_url, params=payload)
        
        # save to file        
        with open(out_file,'wb') as FLE:
            FLE.write(resp.content)
        
        # check for file 
        if not os.path.isfile(out_file):
            # No file
            msg = f'{out_file} failed to download.'
            log.warning(msg)
        else:
            fs = os.path.getsize(out_file)/1e6
            if fs < min_size:
                msg = f"{out_file} is {fs:.1} Mb so is probably "
                msg += "'Access Denied'"
                log.warning(msg)
                
                if delete_small:
                    os.remove(out_file)
            else:
                log.info(f"Complete: {out_file} ({fs:.1} Mb)")
                outlist.append(out_file)
    
    return outlist


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


#########
# Logging
#########
def log_function_arguments(LOGFILE, frame, func='func', verbose=True):
    """
    Log local variables, e.g., parameter arguements to a file
    
    Parameters
    ----------
    LOGFILE : str or None
        Output file.  If `None`, then force `verbose=True`.
    
    frame : `~inspect.currentframe()`
        Namespace object.
    
    func : str
        Function name to use
        
    verbose : bool
        Print messaage to stdout.
    
    """
    args = inspect.getargvalues(frame).locals
    args.pop('frame')
    for k in list(args.keys()): 
        if hasattr(args[k], '__builtins__'):
            args.pop(k)
    
    if func is not None:
        logstr = '\n{0}(**{1})\n'
    else:
        logstr = '\n{1}'
        
    logstr = logstr.format(func, args)
    msg = log_comment(LOGFILE, logstr, verbose=verbose, show_date=True)
    return msg


def log_comment(LOGFILE, comment, verbose=False, show_date=False, mode='a'):
    """
    Log a message to a file, optionally including a date tag
    """
    import time
        
    if show_date:
        msg = '\n# ({0})\n'.format(time.ctime())
    else:
        msg = ''
        #fp.write('\n# ({0})\n'.format(time.ctime()))
    
    msg += '{0}\n'.format(comment)
    
    if LOGFILE is not None:
        fp = open(LOGFILE, mode)
        fp.write(msg)
        fp.close()
    
    if verbose:
        print(msg[:-1])


def log_exception(LOGFILE, traceback, verbose=True, mode='a'):
    """
    Log exception information to a file, or print to screen
    
    Parameters
    ----------
    LOGFILE : str or None
        Output file.  If `None`, then force `verbose=True`.
    
    traceback : builtin traceback module
        Exception traceback, from global `import traceback`.
    
    verbose : bool
        Print exception to stdout.
    
    mode : 'a', 'w'
        File mode on `open(LOGFILE, mode)`, i.e., append or write.
    
    """
    import time
    
    trace = traceback.format_exc(limit=2)
    log = '\n########################################## \n# ! Exception ({0})\n'.format(time.ctime())
    log += '#\n# !'+'\n# !'.join(trace.split('\n'))
    log += '\n######################################### \n\n'
    if verbose | (LOGFILE is None):
        print(log)
        
    if LOGFILE is not None:
        fp = open(LOGFILE, mode)
        fp.write(log)
        fp.close()

##########
# SREGIONs
##########
        
def polygon_to_sregion(poly):
    """
    Convert `shapely.Polygon` vertices to an SREGION string
    """
    try:
        xy = np.array(poly.exterior.xy).T.flatten()
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
        
        # Round to nearest arcmin (4 seconds in RAh)
        >>> from mastquery.utils import radec_to_targname
        >>> print(radec_to_targname(ra=ra, dec=dec, round_arcsec=(4,60),
        ...                    targstr='j{rah}{ram}{ras}{sign}{ded}{dem}'))
        j020204m1010
        
        # Full precision
        >>> targstr = 'j{rah}{ram}{ras}.{rass}{sign}{ded}{dem}{des}.{dess}'
        >>> print(radec_to_targname(ra, dec,round_arcsec=(0.0001, 0.0001),
        ...                         precision=3, targstr=targstr))
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
    
def get_mw_dust(ra, dec, **kwargs):
    """
    Wrapper around functions to try to query for the MW E(B-V)
    """
    try:
        ebv = get_dustmaps_dust(ra, dec, web=True)
        return ebv
    except:
        pass
        
    try:
        ebv = get_dustmaps_dust(ra, dec, web=False)
        return ebv
    except:
        pass
    
    try:
        ebv = get_irsa_dust(ra, dec, **kwargs)
        return ebv
    except:
        pass
    
    # All failed
    return 0.00
    
def get_dustmaps_dust(ra, dec, web=True, **kwargs):
    "Use https://github.com/gregreen/dustmaps"
    
    from dustmaps.sfd import SFDQuery, SFDWebQuery
    from astropy.coordinates import SkyCoord
    
    coords = SkyCoord(ra, dec, unit='deg', frame='icrs')
    
    if web:
        sfd = SFDWebQuery()
    else:
        sfd = SFDQuery()
        
    ebv = sfd(coords)
    return ebv
    
def get_irsa_dust(ra, dec, type='SandF', **kwargs):
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
    import urllib.request as requester
        
    from astropy.table import Table
    from lxml import objectify
    
    query = 'http://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr={0:.4f}+{1:.4f}+equ+j2000'.format(ra, dec)
    
    req = requester.Request(query)
    response = requester.urlopen(req)
    resp_text = response.read().decode('utf-8')
    
    root = objectify.fromstring(resp_text)
    stats = root.result.statistics

    if type == 'SFD':
        return float(str(stats.refPixelValueSFD).split()[0])
    else:
        return float(str(stats.refPixelValueSandF).split()[0])
        
