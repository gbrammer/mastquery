"""
JWST queries

https://mast.stsci.edu/api/v0/_jwst_inst_keywd.html

https://mast.stsci.edu/api/v0/_services.html#MastScienceInstrumentKeywordsNircam

"""
import os
import logging
import numpy as np
import yaml

import json

from tqdm import tqdm
import astropy.table
import astropy.time

from shapely.geometry import Point

from sregion import SRegion, patch_from_polygon

try:
    import pysiaf
except ImportError:
    print('`import pysiaf` failed, necessary for JWST processing')
    
ALL_COLUMNS = ['ArchiveFileID', 'filename', 'fileSetName', 'productLevel',
               'act_id', 'apername', 'asnpool', 'asntable', 'bartdelt',
               'bendtime', 'bkgdtarg', 'bkglevel', 'bkgsub', 'bmidtime',
               'bstrtime', 'category', 'cont_id', 'datamode', 'dataprob',
               'date', 'date_mjd', 'date_end', 'date_end_mjd', 'date_obs',
               'date_obs_mjd', 'detector', 'drpfrms1', 'drpfrms3', 'duration',
               'effexptm', 'effinttm', 'eng_qual', 'exp_type', 'expcount',
               'expend', 'expmid', 'exposure', 'expripar', 'expstart',
               'fastaxis', 'filter', 'frmdivsr', 'gdstarid',
               'groupgap', 'gs_dec', 'gs_mag', 'gs_order', 'gs_ra',
               'gsendtim', 'gsendtim_mjd', 'gsstrttm', 'gsstrttm_mjd',
               'gs_udec', 'gs_umag', 'gs_ura', 'helidelt', 'hendtime',
               'hga_move', 'hmidtime', 'hstrtime', 'instrume', 'intarget',
               'is_psf', 'lamp', 'mu_dec', 'mu_epoch', 'mu_epoch_mjd',
               'mu_ra', 'nexposur', 'nextend', 'nframes', 'ngroups', 'nints',
               'nresets', 'nrststrt', 'nsamples', 'numdthpt', 'nwfsest',
               'obs_id', 'observtn', 'obslabel', 'origin', 'pcs_mode',
               'pi_name', 'pps_aper', 'prd_ver', 'program', 'prop_dec',
               'prop_ra', 'pwfseet', 'readpatt', 'sca_num', 'scicat',
               'sdp_ver', 'selfref', 'seq_id', 'slowaxis', 'subarray',
               'subcat', 'subsize1', 'subsize2', 'substrt1', 'substrt2',
               'targ_dec', 'targ_ra', 'targname', 'targoopp', 'targprop',
               'targtype', 'targudec', 'targura', 'telescop', 'template',
               'tframe', 'tgroup', 'timesys', 'title', 'tsample', 'tsovisit',
               'visit', 'visit_id', 'visitend', 'visitend_mjd', 'visitgrp',
               'visitsta', 'visitype', 'vststart', 'vststart_mjd', 'zerofram',
               'errtype', 'wtype', 'datamodl', 'exp_only', 'exsegnum',
               'exsegtot', 'intstart', 'intend', 'date_beg', 'date_beg_mjd',
               'obsfoldr', 'sctarate', 'opmode', 'osf_file', 'masterbg',
               'scatfile', 'srctyapt', 'tcatfile', 'texptime', 'patt_num',
               'pattsize', 'patttype', 'pridtpts', 'subpxpts', 'crowdfld',
               'engqlptg', 'oss_ver', 'noutputs', 'gs_v3_pa', 'dirimage',
               'pixfrac', 'pxsclrt', 'segmfile', 'va_dec', 'va_ra',
               's_region', 'cal_ver', 'cal_vcs', 'crds_ctx', 'crds_ver',
               'focuspos', 'fwcpos', 'pupil', 'pwcpos', 'nrimdtpt',
               'fileSize', 'checksum', 'ingestStartDate',
               'ingestStartDate_mjd', 'ingestCompletionDate',
               'ingestCompletionDate_mjd', 'FileTypeID', 'publicReleaseDate',
               'publicReleaseDate_mjd', 'isRestricted', 'isItar', 'isStale',
               'FileSetId', 'dataURI']

DEFAULT_COLUMNS = {
  'NIS':['ArchiveFileID', 'obs_id', 'filename', 'fileSetName',
         'visit', 'visit_id', 'observtn', 'productLevel', 'apername',
         'gs_v3_pa', 'pps_aper', 'template',
         'bkglevel', 'date_mjd', 'expstart', 'expend', 'exposure',
         'telescop', 'instrume', 'detector', 'subarray',
         'filter', 'pupil', 'fwcpos', 'pwcpos', 'lamp',
         'readpatt', 'nexposur', 'nframes', 'ngroups', 'tframe',
         'nints', 'nresets', 'effexptm as exptime',
         'pi_name', 'program', 'title', 'category',
         'targname', 'prop_ra', 'prop_dec',
         'patttype', 'patt_num', 'pattsize',
         'pridtpts', 'subpxpts',
         'xoffset', 'yoffset',
         's_region', 'crowdfld', 'pcs_mode', 'expripar',
         'publicReleaseDate_mjd as ReleaseDate',
         'FileSetId', 'dataURI','fileSize'],
  'NRC':['ArchiveFileID', 'obs_id', 'filename', 'fileSetName',
         'visit', 'visit_id', 'observtn', 'productLevel', 'apername',
         'gs_v3_pa', 'pps_aper', 'template',
         'bkglevel', 'date_mjd', 'expstart', 'expend', 'exposure',
         'telescop', 'instrume', 'module', 'channel', 'detector', 'subarray',
         'filter','pupil', 'lamp',
         'readpatt', 'nexposur', 'nframes', 'ngroups', 'tframe',
         'nints', 'nresets', 'effexptm as exptime',
         'pi_name', 'program', 'title', 'category',
         'targname', 'prop_ra', 'prop_dec', 
         'patttype', 'patt_num', 'pattsize',
         'pridtype', 'pridtpts', 'subpxpat', 'subpxpts', 'smgrdpat', 
         'xoffset', 'yoffset',
         's_region', 'crowdfld', 'pcs_mode', 'expripar', 
         'publicReleaseDate_mjd as ReleaseDate',
         'FileSetId', 'dataURI','fileSize'],
  'MIR':['ArchiveFileID', 'obs_id', 'filename', 'fileSetName',
         'visit', 'visit_id', 'observtn', 'productLevel', 'apername',
         'gs_v3_pa', 'pps_aper', 'template',
         'bkglevel', 'date_mjd', 'expstart', 'expend', 'exposure',
         'telescop', 'instrume', 'channel', 'band', 'detector', 'subarray',
         'filter', 'lamp',
         'readpatt', 'nexposur', 'nframes', 'ngroups', 'tframe',
         'nints', 'nresets', 'effexptm as exptime',
         'mirngrps', 'mirnfrms', 
         'pi_name', 'program', 'title', 'category',
         'targname', 'prop_ra', 'prop_dec', 
         'patttype', 'patt_num', 'pattsize',
         'pridtpts', 'subpxpts',
         'xoffset', 'yoffset',
         'dithdirc', 'dithopfr', 'dithpnts',
         's_region', 'crowdfld', 'pcs_mode', 'expripar', 
         'publicReleaseDate_mjd as ReleaseDate',
         'FileSetId', 'dataURI','fileSize'],
  'NRS':['ArchiveFileID', 'obs_id', 'filename', 'fileSetName',
         'visit', 'visit_id', 'observtn', 'productLevel', 'apername',
         'gs_v3_pa', 'pps_aper', 'template',
         'bkglevel', 'date_mjd', 'expstart', 'expend', 'exposure',
         'telescop', 'instrume', 'detector', 'subarray',
         'nrs_norm', 'nrs_ref', 
         'filter', 'fxd_slit', 'grating', 'lamp',
         'msaconid', 'msametfl', 'msametid', 'msastate', 'preimage',
         'readpatt', 'nexposur', 'nframes', 'ngroups', 'tframe',
         'nints', 'nresets', 'effexptm as exptime',
         'pi_name', 'program', 'title', 'category',
         'targname', 'prop_ra', 'prop_dec', 
         'patttype', 'patt_num', 'pattsize',
         'pridtpts', 'subpxpat', 'subpxpts', 'nod_type', 
         'xoffset', 'yoffset',
         's_region', 'crowdfld', 'pcs_mode', 'expripar', 
         'publicReleaseDate_mjd as ReleaseDate',
         'FileSetId', 'dataURI','fileSize'],
}

REQUEST = {'service': 'Mast.Jwst.Filtered.Niriss',
                'params': {'columns': '*', 'filters': []},
                'format': 'json',
                'pagesize': 100000}
                
SERVICES = {'NIS':'Mast.Jwst.Filtered.Niriss', 
            'NRC':'Mast.Jwst.Filtered.Nircam',
            'MIR':'Mast.Jwst.Filtered.Miri',
            'NRS':'Mast.Jwst.Filtered.Nirspec', 
            'FGS':'Mast.Jwst.Filtered.Fgs',
            'GS':'Mast.Jwst.Filtered.GuideStar'}

FILTER_EXAMPLES = [{'paramName':
                    'productLevel', 'values': ['2a','2b']},
                   # {'paramName': 'filetype', 
                   #  'values': ['countrate','calibrated']},
                   {'paramName': 'program', 'values': ['1076']},
                   {'paramName': 'expstart', 
                    'values': [{"min":59672.0, "max":59690.}]},
                   {'paramName': 'filename',
                    'values': [], 'freeText':"%rate.fits"} 
                   ]

CALIB_FILTERS = [{'paramName':
                    'productLevel', 'values': ['2a','2b']},
                 # {'paramName': 'filetype', 
                 #    'values': ['countrate','calibrated']},
                 {'paramName': 'filename', 'values': [],
                    'freeText':"%_[cr]a%[le].fits"}
                   ]

FULL_SUBARRAY = [{'paramName': 'subarray', 'values': ['FULL']}]

FINE_GUIDE = [{'paramName': 'pcs_mode', 'values': ['FINEGUIDE']}]


def make_query_filter(column='filters', values=['F200W'], text=None, range=None):
    """
    Generate a query filter for the MAST API
    
    
    """
    filt = {'paramName': column}
    if text is not None:
        filt['freeText'] = text
        filt['values'] = []
    elif range is not None:
        filt['values'] = [{'min':range[0], 'max':range[1]}]
    else:
        filt['values'] = values
    
    return [filt]
    
    
def make_program_filter(progs):
    """
    Helper for generating a filter on a list of program ids
    
    Parameters
    ----------
    progs : list
        `program` ID numbers
    
    Returns
    -------
    pfilter : dict
        MAST CAOM filter dict for program IDs.
        
    """
    return make_query_filter(column='program', 
                        values=[f'{int(p)}' for p in progs])


def query_all_jwst(instruments=['NRC','NIS','NRS','MIR'], columns=None, rd=None, sort='expstart', fix=True, **kwargs):
    """
    Combined query for all instruments
    
    Parameters
    ----------
    instruments : list
        Instruments to query (``NRC``, ``NIS``, ``NRS``, ``MIR``)
    
    columns : list, None
        Columns to query.  If None, then get all (`'*`') columns
    
    rd : (float, float), (float, float, float), None
        If a 2-tuple, then return exposures where `s_region` contains `rd`.  
        If a 3-tuple, then find fields where `(targ_ra, targ_dec)` within 
        `rd[2]` arcmin of `rd[0,1]`
    
    sort : float
        Column for sorting the output table
    
    fix : bool
        Fixes to make make instrument queries more compatible with the 
        high-level MAST query results (e.g., HST).
        
        >>> mastquery.jwst.set_missing_footprints(jw)
        >>> mastquery.jwst.fix_jwst_sregions(jw)
        >>> mastquery.jwst.set_footprint_centroids(jw)
        >>> mastquery.jwst.match_dataportal_columns(jw)
    
    kwargs : dict
        Keyword args passed to instrument-level `mastquery.jwst.query_jwst` 
        queries, e.g., `filters`, `recent_days`.
        
    Returns
    -------
    jw : table
        Query result
        
    """                        
        
    log = logging.getLogger()
    log.debug('kwargs')
    log.debug(yaml.dump(kwargs))
    
    res = []
    for inst in instruments:
        ri = query_jwst(instrument=inst, columns=columns, **kwargs)
        if len(ri) > 0:
            res.append(ri)
    
    if len(res) == 0:
        return None
        
    for t in res:
        empty = []
        for c in t.colnames:
            if t[c].dtype in [object]:
                fill = np.in1d(t[c], [None])
                if fill.sum() == len(t):
                    empty.append(c)
                    continue
                
                t[c][fill] = t[c][~fill][0]*0
                t[c] = np.ma.masked_array(t[c], mask=fill,
                                  dtype=t[c][~fill][0].__class__)

            if hasattr(t[c], 'mask'):
                if t[c].mask.sum() == len(t):
                    empty.append(c)
                
        t.remove_columns(empty)
    
    jw = astropy.table.vstack(res)
    
    if rd is not None:
        from grizli import utils
        
        if len(rd) == 2:
            # Contains point
            pt = Point(*rd)
            
            tab = utils.GTable()
            tab['ra'], tab['dec'] = [rd[0]], [rd[1]]
            idx, dr = tab.match_to_catalog_sky(jw, 
                                         other_radec=('prop_ra', 'prop_dec'))
            
            sel = dr.value < 60*20
            ix = np.where(sel)[0]
            for i in ix:
                fp = jw['s_region'][i]
                try:
                    sr = utils.SRegion(fp)
                    sel[i] = sr.shapely[0].contains(pt)
                except ValueError:
                    #print('err', fp)
                    sel[i] = False
            
            sel = np.array(sel)
            
        else:
            # Separation in arcmin
            tab = utils.GTable()
            tab['ra'], tab['dec'] = [rd[0]], [rd[1]]
            idx, dr = tab.match_to_catalog_sky(jw, 
                                         other_radec=('prop_ra', 'prop_dec'))
            
            sel = dr.value < rd[2]*60
        
        jw = jw[sel]
    
    if sort is not None:
        if sort in jw.colnames:
            log.debug(f"Sort on column {sort}")
            
            so = np.argsort(jw[sort])
            jw = jw[so]
        else:
            log.warning(f"sort='{sort}' column not found")
    
    if fix:
        log.info('Apply fixes to JWST query')
        set_missing_footprints(jw)
        jw = fix_jwst_sregions(jw)
        set_footprint_centroids(jw)
        match_dataportal_columns(jw)
    
    if len(jw) > 0:
        log.info(f'Full query on {instruments} N={len(jw)}')
    else:
        log.warning('Nothing found!')
        
    return jw


def fix_jwst_sregions(res):
    """
    Fix s_region polygons based on pointing keywords, mostly for MIRIM_FULL
    """
    import pysiaf

    siaf = {}
    for ins in ['MIRI','NIRISS','NIRCAM','NIRSPEC','FGS']:
        siaf[ins] = pysiaf.Siaf(ins)

    s_region = []
    
    if hasattr(res['s_region'], 'mask'):
        print(f"Remove {res['s_region'].mask.sum()} rows with missing s_region")
        res = res[~res['s_region'].mask]
        
    res['orig_s_region'] = res['s_region']
    
    _iter = range(len(res))
    if len(res) > 2000:
        _iter = tqdm(_iter)
        
    for i in _iter:
        if not res['xoffset'][i]:
            s_region.append(res['s_region'][i])
            continue
        
        aper_name = res['apername'][i]
                
        if res['pps_aper'][i] not in siaf[res['instrume'][i]].apertures:
            s_region.append(res['s_region'][i])
            continue
            
        ap_ref = siaf[res['instrume'][i]].apertures[res['pps_aper'][i]]

        if aper_name == 'MIRIM_FULL':
            apers = ['MIRIM_ILLUM','MIRIM_CORONLYOT']
        else:
            apers = [aper_name]

        sreg = []

        for aper_i in apers:
            ap = siaf[res['instrume'][i]].apertures[aper_i]

            idl = np.array(ap.corners('idl', rederive=False)).T
            idl += np.array([-res['xoffset'][i], -res['yoffset'][i]])
            asec = np.array(ap.idl_to_tel(*idl.T)).T
            asec -= np.array([ap_ref.V2Ref, ap_ref.V3Ref])

            theta = np.deg2rad(res['gs_v3_pa'][i] + 0)
            _mat = np.array([[np.cos(theta), -np.sin(theta)],
                             [np.sin(theta), np.cos(theta)]])

            arot = asec.dot(_mat)[::-1,:]

            cosd = np.array([np.cos(res['targ_dec'][i]/180*np.pi), 1])

            asky = arot / 3600. / cosd 
            asky += np.array([res['targ_ra'][i], res['targ_dec'][i]])

            sr = SRegion(asky)
            sreg.append(sr.s_region)

        s_region.append(' '.join(sreg)) #sr.s_region

    res['s_region'] = s_region
    return res


def set_missing_footprints(res):
    """
    Set missing s_region entries, mostly NIRISS
    """
    log = logging.getLogger()
    log.debug(f'N = {len(res)}')
    
    # Set Grism footprints
    sprev = {}

    for i, sn in enumerate(res['s_region']):
        try:
            sr = SRegion(sn)
            sprev[res['instrume'][i]] = sn
        except:
            test = 'GR' in res['filter'][i]
            if 'pupil' in res.colnames:
                test |= ('GR' in res['pupil'][i])

            if test:
                #use previous footprint for GRISM
                try:
                    res[i]['s_region'] = sprev[res['instrume'][i]]
                except:
                    # Can fail e.g. for PASSAGE with grism-only visits
                    # make a circle with radius 1 arcmin
                    _ra, _dec = res['targ_ra'][i], res['targ_dec'][i]
                    dummy = f"CIRCLE {_ra:.6f} {_dec:.6f} 0.017"
                    res[i]['s_region'] = dummy
                    
                pass

def set_footprint_centroids(res):
    """
    Set ra/dec to s_region centroid
    """
    log = logging.getLogger()
    log.debug(f'N = {len(res)}')

    res['ra'] = 0.
    res['dec'] = 0.

    for i, s in enumerate(res['s_region']):
        sr = SRegion(s)
        res['ra'][i], res['dec'][i] = sr.centroid[0]


def match_dataportal_columns(res):
    """
    Match some columns in the jwst query to the output from the dataportal 
    query
    """
    log = logging.getLogger()
    log.debug(f'N = {len(res)}')
    
    # Fill columns with single value
    fill_cols = {'obs_collection':'JWST',
                 'wavelength_region':'IR', 
                 'project':'JWST',
                 'provenance_name':'CAOM'}

    for c in fill_cols:
        if c in res.colnames:
            #print(c)
            res.remove_column(c)

        res[c] = fill_cols[c]
    
    res['orig-filter'] = res['filter']
    
    # Translate columns
    col_translate = {'targname':'target', 
                     'filter-pupil':'filter',
                     'pi_name': 'proposal_pi',
                     'program': 'proposal_id',
                     'instrume': 'instrument_name', 
                     's_region': 'footprint',
                     'expstart':'t_min',
                     'expend':'t_max',
                     'effexptm':'exptime',
                     'dataURI':'dataURL',
                     'title':'obs_title',
                     'publicReleaseDate_mjd':'t_obs_release',
                     'category':'proposal_type'}

    for c in col_translate:
        if c in res.colnames:
            res[col_translate[c]] = res[c]

    res.remove_column('filter')
    res['filter'] = res['inst-mode']

    for i, (inst, f) in enumerate(zip(res['instrume'], res['filter-pupil'])):
        f_i = f.replace('-','.').replace('.CLEAR','').replace('CLEAR.','')
        # Combine grisms and filters for WFSS so that they are
        # associated together
        f_i = f_i.replace('GR150R.','').replace('GR150C.','')
        f_i = f_i.replace('.GRISMR','').replace('.GRISMC','')

        inst_i = {'NIRCAM':'NC', 
                  'NIRISS':'NI',
                  'MIRI':'MI',
                  'NIRSPEC':'NS',
                  'FGS':'FGS'}[inst]
                  
        if 'GR' in f_i:
            res['filter'][i] = f'{f_i}'
        else:
            res['filter'][i] = f'{inst_i}.{f_i}'


def query_jwst(instrument='NIS', columns='*', filters=CALIB_FILTERS+FULL_SUBARRAY+FINE_GUIDE, extra={'format':'json', 'pagesize': 100000}, recent_days=None, rates_and_cals=False, extensions=['rate', 'cal'], verbose=True):
    """
    Query 
    """    
    log = logging.getLogger()
    
    from .utils import mastQuery, mastJson2Table
    from .utils import new_mastJson2Table, new_mast_query
        
    if recent_days is not None:
        now = astropy.time.Time.now().mjd
        filters = [f for f in filters]
        filters.append({'paramName': 'expstart', 
                        'values': [{"min":float(now-recent_days),
                                    "max":float(now+1)}]})
    
    if instrument.upper() not in SERVICES:
        msg = "Valid options for `instrument` are 'NIS', 'NRC', 'NRS', 'MIR' "
        raise ValueError(msg)
    
    if columns is None:
        columns = ','.join(DEFAULT_COLUMNS[instrument.upper()])
        
    request = {'service': SERVICES[instrument.upper()]}
    request['params'] = {'columns':columns, 
                         'filters':filters}
    for k in extra:
        request[k] = extra[k]
    
    try:

        log.info(f'Query JWST {instrument}')
        log.debug(f'Request: ')
        log.debug(yaml.dump(request))
        
        head, content = new_mast_query(request)
    except:
        head, content = mastQuery(request)
    
    try:
        tab = new_mastJson2Table(content)
    except KeyError:
        json_data = json.loads(content)
        if json_data['status'] == 'ERROR':
            raise ValueError(json_data['msg'])
        
        tab = mastJson2Table(json_data)
        
    if 'filename' in tab.colnames:
        so = np.argsort(tab['filename'])
        tab = tab[so]
    elif 'fileSetName' in tab.colnames:
        so = np.argsort(tab['fileSetName'])
        tab = tab[so]
    
    if (extensions is not None) & ('filename' in tab.colnames):
        ext = [f.split('_')[-1].split('.')[0] for f in tab['filename']]
        in_ext = np.in1d(ext, extensions)
        if in_ext.sum() == 0:
            log.warning(f'No files with extensions {extensions} found.')
            log.warning(f'Available extensions are {np.unique(ext).tolist()}')
            
        tab = tab[in_ext]
        
    if (recent_days is not None) & ('expstart' in tab.colnames):
        tab['dt'] = now - tab['expstart']
        tab['dt'].format = '.2f'
    
    if ('program' in tab.colnames) & ('pi_name' in tab.colnames):
        tab['prog_pi'] = [f'{ca:>4}-{pr} {pi}'
                          for ca, pr, pi in zip(tab['category'], 
                                                tab['program'],
                                                tab['pi_name'])]
    
    if ('filter' in tab.colnames) & ('pupil' in tab.colnames):
        tab['filter-pupil'] = [f'{f}-{p}'.replace('---','')
                               for f, p in zip(tab['filter'], tab['pupil'])]
        
        tab['inst-mode'] = [f'{ii}-{f}-{p}'.replace('---','').replace('-CLEAR','')
                               for ii, f, p in
                            zip(tab['instrume'], tab['filter'], tab['pupil'])]
    elif 'filter' in tab.colnames:
        tab['filter-pupil'] = [f'{f}'.replace('---','')
                               for f in tab['filter']]
        
        tab['inst-mode'] = [f'{ii}-{f}'.replace('---','').replace('-CLEAR','')
                               for ii, f in
                            zip(tab['instrume'], tab['filter'])]
    
    log.info(f'Query JWST {instrument} N={len(tab)}\n')
    
    return tab


def query_guidestar_log(mjd=None, program=None, exp_type=['FGS_FINEGUIDE']):
    """
    Query GuideStar log files
    """
    filters = []
    filters += make_query_filter('exp_type', values=exp_type)
    
    if program is not None:
        filters += make_program_filter([program])
        
    if mjd is not None:
        filters += make_query_filter('datastrt', range=mjd)
    
    gs = query_jwst(instrument='GS',
                         columns='*',
                         filters=filters,
                         extensions=['cal'])
    
    if len(gs) > 0:
        if ('expstart' not in gs.colnames) & ('datastrt' in gs.colnames):
            gs['expstart'] = gs['datastrt']

        if ('expend' not in gs.colnames) & ('dataend' in gs.colnames):
            gs['expend'] = gs['dataend']
            
    return gs

