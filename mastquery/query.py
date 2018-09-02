"""
Demo of querying the ESA HST archive
"""

import numpy as np

MASTER_COLORS = {'G102':'#1f77b4',
'F125W':'#ff7f0e',
'F160W':'#2ca02c',
'G141':'#d62728',
'F140W':'#9467bd',
'F105W':'#8c564b',
'F775W':'#8c564b'}

# DEFAULT_RENAME = {'OPTICAL_ELEMENT_NAME':'FILTER',
#                   'EXPOSURE_DURATION':'EXPTIME',
#                   #'INSTRUMENT_NAME':'INSTRUMENT', 
#                   #'DETECTOR_NAME':'DETECTOR', 
#                   'STC_S':'FOOTPRINT', 
#                   #'SPATIAL_RESOLUTION':'PIXSCALE', 
#                   'TARGET_NAME':'TARGET', 
#                   'SET_ID':'VISIT'}

DEFAULT_RENAME = {'t_exptime':'exptime',
                  'target_name':'target',
                  's_region':'footprint', 
                  's_ra':'ra',
                  's_dec':'dec',
                  'filters':'filter'}
                  
DEFAULT_COLUMN_FORMAT = {'t_min':'.4f',
           't_max':'.4f',
           'exptime':'.0f',
           'ra':'.6f',
           'dec':'.6f'}

# Don't get calibrations.  Can't use "INTENT LIKE 'SCIENCE'" because some 
# science observations are flagged as 'Calibration' in the ESA HSA.
DEFAULT_QUERY = {'project':['HST'], 'intentType':['science'], 'mtFlag':[False]}

# DEFAULT_EXTRA += ["TARGET.TARGET_NAME NOT LIKE '{0}'".format(calib) 
#                  for calib in ['DARK','EARTH-CALIB', 'TUNGSTEN', 'BIAS',
#                                'DARK-EARTH-CALIB', 'DARK-NM', 'DEUTERIUM',
#                                'INTFLAT', 'KSPOTS', 'VISFLAT']]

INSTRUMENT_DETECTORS = {'WFC3/UVIS':'UVIS', 'WFC3/IR':'IR', 'ACS/WFC':'WFC', 'ACS/HRC':'HRC', 'WFPC2':'1', 'STIS/NUV':'NUV-MAMA', 'STIS/ACQ':'CCD'}

def run_query(box=None, proposid=[13871], instruments=['WFC3/IR'], filters=[], extensions=['RAW','C1M'], base_query=DEFAULT_QUERY, maxitems=100000, timeout=300, rename_columns=DEFAULT_RENAME, lower=True, sort_column=['OBSERVATION_ID'], remove_tempfile=True, get_query_string=False, quiet=True):
    """
    
    Optional position box query:
        box = [ra, dec, radius] with ra and dec in decimal degrees and radius
        in arcminutes.
    
    Some science observations are flagged as INTENT = Calibration, so may have
    to run with extra=[] for those cases and strip out true calibs another
    way.
    
    """
    import os
    import tempfile   
    import urllib.request
    import time
    
    from astropy.table import Table
    from . import utils
    
    if quiet:
        utils.set_warnings(numpy_level='ignore', astropy_level='ignore')
        
    request = {'params':{"columns":"*"},
               'format':'json',
               'pagesize':maxitems,
               'page':1,
               'removenullcolumns':True,
               'timeout':timeout,
               'removecache':True}

    # Box search around position
    if (box is not None):
        ra, dec, radius = box
        box_str = "{0}, {1}, {2}".format(ra, dec, radius/60.)
        request['params']['position'] = box_str
        request['service'] = 'Mast.Caom.Filtered.Position'
    else:
        request['service'] = 'Mast.Caom.Filtered'

    query = {}
    for k in base_query:
        query[k] = base_query[k]

    if len(proposal_id) > 0:
        query['proposal_id'] = ['{0}'.format(p) for p in proposal_id]
        
    if len(filters) > 0:
        query['filters'] = [filters]

    if len(instruments) > 0:
        query['instrument_name'] = [instruments]
    
    # Translate to filter format
    query_list = []
    for k in query:
        # Range
        if k.startswith('!'):
            query_list.append({"paramName":k[1:], 'min':query[k][0], 'max':query[k][1]})
        elif k.startswith('*'):
            query_list.append({"paramName":k[1:], values:[], 'freeText':query[k]})
        else:
            query_list.append({"paramName":k[1:], values:query[k]})
            
    request['params']['filters'] = query_list
            
    if get_query_string:
        return request
    
    # Run the query
    headers, outString = utils.mastQuery(request)
    
    try:
        outData = json.loads(outString)
        tab = utils.mastJson2Table(outData)
    except:
        print('Failed to initialize the table (result={0})'.format(outString))
        return False
        
    tab.meta['query'] = json.dumps(request), 'Query input'
    tab.meta['qtime'] = time.ctime(), 'Query timestamp'
    
    if len(tab) == 0:
        return tab
    
    # Sort
    tab.sort('proposal_id')
    
    # Add coordinate name
    if 'ra' in tab.colnames:
        jtargname = [utils.radec_to_targname(ra=tab['ra'][i], dec=tab['dec'][i], scl=6) for i in range(len(tab))]
        tab['jtargname'] = jtargname
    
    fix_byte_columns(tab)
        
    for c in rename_columns:
        if c in tab.colnames:
            tab.rename_column(c, rename_columns[c])
        
    # cols = tab.colnames
    # if ('instrument' in cols) & ('detector' in cols):
    #     tab['instdet'] = ['{0}/{1}'.format(tab['instrument'][i], tab['detector'][i]) for i in range(len(tab))]
            
    #tab['OBSERVATION_ID','orientat'][so].show_in_browser(jsviewer=True)
    
    set_default_formats(tab)
    
    return tab

def fix_byte_columns(tab):
    for col in tab.colnames:
        try:
            if tab[col].dtype == np.dtype('O'):
                strcol = [item.decode('utf-8') for item in tab[col]]
                tab.remove_column(col)
                tab[col] = strcol
        except:
            pass
            
def add_postcard(table, resolution=256):
    
   url = ['http://archives.esac.esa.int/ehst-sl-server/servlet/data-action?OBSERVATION_ID={0}&RETRIEVAL_TYPE=POSTCARD&RESOLUTION={1}'.format(o, resolution) for o in table['observation_id']]
   
   img = ['<a href="{0}"><img src="{0}"></a>'.format(u) for u in url]
   table['postcard'] = img
   
   return True
   
   if False:
       tab = grizli.utils.GTable(table)
       tab['observation_id','filter','orientat','postcard'][tab['visit'] == 1].write_sortable_html('tab.html', replace_braces=True, localhost=True, max_lines=10000, table_id=None, table_class='display compact', css=None)
        
def old_parse_polygons(polystr):
    if 'UNION' in polystr.upper():
        spl = polystr[:-1].split('Polygon')[1:]
    else:
        spl = polystr.replace('FK5','ICRS').split('ICRS')[1:]
        
    poly = [np.cast[float](p.split()).reshape((-1,2)) for p in spl]
    return poly

def parse_polygons(polystr):
    if hasattr(polystr, 'decode'):
        spl = polystr.decode('utf-8').strip()[1:-1].split(',')
    else:
        spl = polystr.strip()[1:-1].split(',')

    poly = [np.cast[float](spl).reshape((-1,2))] #[np.cast[float](p.split()).reshape((-1,2)) for p in spl]
    return poly
    
    

def set_default_formats(table, formats=DEFAULT_COLUMN_FORMAT):
    """
    Set default print formats
    """
    DEFAULT_FORMATS = {'start_time_mjd':'.4f',
               'end_time_mjd':'.4f',
               'exptime':'.0f',
               'ra':'.6f',
               'dec':'.6f',
               'ecl_lat':'.6f',
               'ecl_lon':'.6f',
               'gal_lat':'.6f',
               'gal_lon':'.6f',
               'fov_size':'.3f'}#,
               #'pixscale':'.3f'}
    
    for f in formats:
        if f in table.colnames:
            table[f].format = formats[f]
            
def set_orientat_column(table):
    """
    Make a column in the `table` computing each orientation with
    `get_orientat`.
    """
    table['orientat'] = [query.get_orientat(p) for p in table['footprint']]
    table['orientat'].format = '.1f'
    
def get_orientat(polystr='Polygon ICRS 127.465487 18.855605 127.425760 18.853486 127.423118 18.887458 127.463833 18.889591'):
    """
    
    Compute the "ORIENTAT" position angle (PA of the detector +y axis) from an 
    ESA archive polygon, assuming that the first two entries of the polygon 
    are the LL and UL corners of the detector.
    
    """
    from astropy.coordinates import Angle
    import astropy.units as u
    
    try:
        p = parse_polygons(polystr)[0]
    except:
        p = old_parse_polygons(polystr)[0]
        
    dra = (p[1,0]-p[0,0])*np.cos(p[0,1]/180*np.pi)
    dde = p[1,1] - p[0,1]
    
    orientat = 90+np.arctan2(dra, dde)/np.pi*180
    orientat -= 0.24 # small offset to better match header keywords
    
    orientat = Angle.wrap_at(orientat*u.deg, 180*u.deg).value
    
    return orientat
    
def show_footprints(tab, ax=None):
    """
    Show pointing footprints in a plot
    """
    import matplotlib.pyplot as plt
    
    # Show polygons
    mpl_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    filters = np.unique(tab['filter'])

    colors = {}
    
    if ax is None:
        ax = plt
        
    for i, f in enumerate(filters):
        if f in MASTER_COLORS:
            colors[f] = MASTER_COLORS[f]
        else:
            colors[f] = mpl_colors[i % len(mpl_colors)]
        
    for i in range(len(tab)):
        poly = parse_polygons(tab['footprint'][i])#[0]

        for p in poly:
            pclose = np.vstack([p, p[0,:]]) # repeat the first vertex to close
            
            ax.plot(pclose[:,0], pclose[:,1], alpha=0.1, color=colors[tab['filter'][i]])
            
            # Plot a point at the first vertex, pixel x=y=0.
            ax.scatter(pclose[0,0], pclose[0,1], marker='.', color=colors[tab['filter'][i]], alpha=0.1)
    
    return colors
    

