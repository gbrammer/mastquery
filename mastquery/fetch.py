"""
Fetch data directly from ESA Hubble Science Archive
"""

import os

# DEFAULT_PRODUCTS = {'WFC3-IR':['RAW'],
#                     'WFPC2':['C0M','C1M'],
#                     'ACS-WFC':['FLC'],
#                     'WFC3-UVIS':['FLC']}
                    
DEFAULT_PRODUCTS = {'WFC3/IR':['RAW'],
                    'WFPC2/WFC':['C0M','C1M'],
                    'WFPC2/PC':['C0M','C1M'],
                    'ACS/WFC':['FLC'],
                    'WFC3/UVIS':['FLC']}


def make_curl_script(table, level=None, script_name=None, inst_products=DEFAULT_PRODUCTS, skip_existing=True, s3_sync=False, output_path='./'):
    """
    Generate a "curl" script to fetch products from the ESA HSA
    
    Parameters
    ----------
    table : `~astropy.table.Table`
        Table output from `~hsaquery.query` scripts.
        
    level : str
        Specific data type to retrieve (e.g., 'FLT', 'RAW', 'SPT', etc.).  
        If `None`, then retrieve the following:
            'RAW'       for WFC3/IR
            'FLC'       for ACS/WFC and WFC3/UVIS
            'COM'+'C1M' for WFPC2
            
    script_name : str or None
        If a string, then save the curl commands to a file.
    
    Returns
    -------
    curl_list : list
        List of curl commands.
    
    """
    import tempfile
    import glob
    import os
    
    BASE_URL = 'http://archives.esac.esa.int/ehst-sl-server/servlet/data-action?ARTIFACT_ID=' #J6FL25S4Q_RAW.FITS'
    BASE_URL = 'http://hst.esac.esa.int/ehst-sl-server/servlet/data-action?ARTIFACT_ID=' #J6FL25S4Q_RAW.FITS'
    
    if s3_sync:
        # s3://stpubdata/hst/public/icwb/icwb1iu5q/icwb1iu5q_raw.fits    
        BASE_URL = 's3://stpubdata/hst/public/'
        #BASE_URL = 's3://stpubdata/hst/public/'
        
    if level is None:
        # Get RAW for WFC3/IR, FLC for UVIS and ACS
        curl_list = []
        for i in range(len(table)):
            #inst_det = '{0}/{1}'.format(table['instrument'][i], table['detector'][i]) 
            inst_det = table['instrument_name'][i] #'{0}/{1}'.format(table['instrument'][i], table['detector'][i]) 
            
            if inst_det in inst_products:
                products = inst_products[inst_det]
            else:
                products = ['RAW']
        
            dataset = table['observation_id'][i]
            for product in products:
                if skip_existing:
                    path = '{2}/{0}_{1}.fits*'
                    
                    path_files = glob.glob(path.format(dataset.lower(),
                                           product.lower(), output_path))
                    
                    # Don't fetch RAW if FLT exists
                    if products == ['RAW']:
                        path_files += glob.glob(path.format(dataset.lower(),
                                               'flt', output_path))
                        
                    if len(path_files) > 0:
                        skip = True
                    else:
                        skip = False
                else:
                    skip = False
                    
                if not skip:
                    if s3_sync:
                        curl_list.append(make_s3_command(dataset, product, output_path=output_path, s3_sync=s3_sync))                            
                    else:
                        curl_list.append('curl {0}{1}_{2}.FITS -o {5}/{3}_{4}.fits.gz'.format(BASE_URL, dataset.upper(), product.upper(), dataset.lower(), product.lower(), output_path))
            
    else:
        if s3_sync:
            curl_list = [make_s3_command(dataset, product, output_path=output_path, s3_sync=s3_sync) for dataset in table['observation_id']]
                
        else:
            curl_list = ['curl {0}{1}_{2}.FITS -o {5}/{3}_{4}.fits.gz'.format(BASE_URL, dataset.upper(), level.upper(), dataset.lower(), level.lower(), output_path) for dataset in table['observation_id']]
    
    if script_name is not None:
        fp = open(script_name, 'w')
        fp.write('\n'.join(curl_list))
        fp.close()
    
    return curl_list


def make_s3_command(dataset, product, output_path='./', s3_sync=True):
    """
    Generate a path to the public "STPUBDATA" S3 Hubble data mirror
    
    Parameters
    ----------
    dataset : str
        Dataset name, like 'IDNCM1AGQ' (or 'idncm1agq').
        
    product : type
        File extension, like 'RAW' (or 'raw').
        
    output_path : type
        Path where to put the file
        
    s3_sync : str, bool
        If 'cp' then run `aws s3 cp`.  If anything else, run `aws s3 sync` 
        excluding all files but including the desired `product`.
    
    Returns
    -------
    cmd : str
        Sync/copy command for s3, e.g., 
        'aws s3 cp --request-payer requester s3://stpubdata/hst/public/idnc/idncm1agq/idncm1agq_raw.fits ./'.
        
    .. warning::
    
    Copying from the STPublic S3 bucket outside of AWS can incur significant 
    charges to an AWS account!
    
    """
    BASE_URL = 's3://stpubdata/hst/public/'
    
    if s3_sync == 'cp':
        cmd = 'aws s3 cp --request-payer requester {0}{1}/{2}/{2}_{3}.fits {4}/'.format(BASE_URL, dataset[:4].lower(), dataset.lower(), product.lower(), output_path)
    else:
        cmd = 'aws s3 sync --request-payer requester --exclude="*.*" --include="*{3}.fits" {0}{1}/{2}/ {4}/'.format(BASE_URL, dataset[:4].lower(), dataset.lower(), product.lower(), output_path)
    
    return cmd


def persistence_products(tab):
    import numpy as np
    try:
        tab['instrument_name'].fill_value = '0'
    except:
        pass
        
    wfc3 =  tab['instrument_name'] == 'WFC3/IR'
    progs = np.unique(tab[wfc3]['proposal_id'])
    persist_files = []
    for prog in progs:
        p = wfc3 & (tab['proposal_id'] == prog)
        visits = np.unique(tab['visit'][p])
        print(visits)
        for visit in visits:
            if isinstance(visit, str):
                vst = visit
            else:
                vst = '{0:02d}'.format(visit)
        
            file_i = 'https://archive.stsci.edu/pub/wfc3_persist/{0}/Visit{1}/{0}.Visit{1}.tar.gz'.format(prog, vst.upper())
            persist_files.append(file_i)

    return persist_files
    

### Direct download from MAST
# https://mast.stsci.edu/api/v0/pyex.html#downloadReq
def directDownload(tab, inst_products=DEFAULT_PRODUCTS, path='./', skip_existing=True):

    import sys
    import os
    import time
    import re
    import json
    
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
        
    # collect the data products
    #result = getCaomProducts()
    #data = result['data']

    # setting up the https connection
    server='mast.stsci.edu'
    conn = httplib.HTTPSConnection(server)

    # dowload the first products
    for i in range(len(tab)):
        products = inst_products[tab[i]['instrument_name']]
        for prod in products:
            ext = '_{0}.fits'.format(prod.lower())
            
            # make file path
            #outPath = "mastFiles/HST/"+tab[i]['obs_id']
            #if not os.path.exists(outPath):
            #    os.makedirs(outPath)
            outPath = os.path.join(path, tab[i]['productFilename'].split('_')[0]+ext)
        
            if os.path.exists(outPath) & skip_existing:
                print('File {0} exists.'.format(outPath))
                continue
                
            # Download the data
            uri = tab[i]['dataURI'].split('_')[0]+ext
            conn.request("GET", "/api/v0/download/file?uri="+uri)
            resp = conn.getresponse()
            fileContent = resp.read()
    
            # save to file
            with open(outPath,'wb') as FLE:
                FLE.write(fileContent)
        
            # check for file 
            if not os.path.isfile(outPath):
                print("ERROR: " + outPath + " failed to download.")
            else:
                print("COMPLETE: ", outPath)

    conn.close()


def fetch_wfpc2_calib(file='g6q1912hu_r4f.fits', path=os.getenv('uref')):
    from stsci.tools import convertwaiveredfits
    try: # Python 3.x
        import http.client as httplib 
    except ImportError:  # Python 2.x
        import httplib
    
    server='mast.stsci.edu'
    conn = httplib.HTTPSConnection(server)
    outPath = os.path.join(path, file)
    uri = 'mast:HST/product/'+file
    
    conn.request("GET", "/api/v0/download/file?uri="+uri)
    resp = conn.getresponse()
    fileContent = resp.read()
    
    # save to file
    with open(outPath,'wb') as FLE:
        FLE.write(fileContent)

    # check for file 
    if not os.path.isfile(outPath):
        print("ERROR: " + outPath + " failed to download.")
        status = False
    else:
        print("COMPLETE: ", outPath)
        status = True
        
    conn.close()
    
    if status:
        try:
            hdu = convertwaiveredfits.convertwaiveredfits(outPath)
            hdu.writeto(outPath.replace('.fits','_c0h.fits'))
        except:
            return True
            
        while 'HISTORY' in hdu[0].header:
            hdu[0].header.remove('HISTORY')

    