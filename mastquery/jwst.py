"""
JWST queries

https://mast.stsci.edu/api/v0/_jwst_inst_keywd.html

https://mast.stsci.edu/api/v0/_services.html#MastScienceInstrumentKeywordsNircam

"""
import os
import numpy as np

ALL_COLUMNS = ['ArchiveFileID', 'filename', 'fileSetName', 'productLevel',
               'act_id', 'apername', 'asnpool', 'asntable', 'bartdelt',
               'bendtime', 'bkgdtarg', 'bkglevel', 'bkgsub', 'bmidtime',
               'bstrtime', 'category', 'cont_id', 'datamode', 'dataprob',
               'date', 'date_mjd', 'date_end', 'date_end_mjd', 'date_obs',
               'date_obs_mjd', 'detector', 'drpfrms1', 'drpfrms3', 'duration',
               'effexptm', 'effinttm', 'eng_qual', 'exp_type', 'expcount',
               'expend', 'expmid', 'exposure', 'expripar', 'expstart',
               'fastaxis', 'filetype', 'filter', 'frmdivsr', 'gdstarid',
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
         'bkglevel', 'date_mjd', 'expstart', 'expend', 'exposure','filetype',
         'telescop', 'instrume', 'detector', 'subarray',
         'filter', 'pupil', 'fwcpos', 'pwcpos', 'lamp',
         'readpatt', 'nexposur', 'nframes', 'ngroups', 'tframe',
         'nints', 'nresets', 'effexptm as exptime',
         'pi_name', 'program', 'title',
         'targname', 'prop_ra', 'prop_dec',
         'patttype', 'patt_num', 'pattsize',
         'pridtype', 'pridtpts', 'subpxpat', 'subpxpts',
         'xoffset', 'yoffset',
         's_region', 'crowdfld', 'pcs_mode', 'expripar',
         'publicReleaseDate_mjd as ReleaseDate',
         'FileSetId', 'dataURI'],
  'NRC':['ArchiveFileID', 'obs_id', 'filename', 'fileSetName',
         'visit', 'visit_id', 'observtn', 'productLevel', 'apername',
         'gs_v3_pa', 'pps_aper', 'template',
         'bkglevel', 'date_mjd', 'expstart', 'expend', 'exposure','filetype',
         'telescop', 'instrume', 'module', 'channel', 'detector', 'subarray',
         'filter','pupil', 'lamp',
         'readpatt', 'nexposur', 'nframes', 'ngroups', 'tframe',
         'nints', 'nresets', 'effexptm as exptime',
         'pi_name', 'program', 'title',
         'targname', 'prop_ra', 'prop_dec', 
         'patttype', 'patt_num', 'pattsize',
         'pridtype', 'pridtpts', 'subpxpat', 'subpxpts', 'smgrdpat', 
         'xoffset', 'yoffset',
         's_region', 'crowdfld', 'pcs_mode', 'expripar', 
         'publicReleaseDate_mjd as ReleaseDate',
         'FileSetId', 'dataURI'],
  'MIR':['ArchiveFileID', 'obs_id', 'filename', 'fileSetName',
         'visit', 'visit_id', 'observtn', 'productLevel', 'apername',
         'gs_v3_pa', 'pps_aper', 'template',
         'bkglevel', 'date_mjd', 'expstart', 'expend', 'exposure','filetype',
         'telescop', 'instrume', 'channel', 'band', 'detector', 'subarray',
         'filter', 'lamp',
         'readpatt', 'nexposur', 'nframes', 'ngroups', 'tframe',
         'nints', 'nresets', 'effexptm as exptime',
         'mirngrps', 'mirnfrms', 
         'pi_name', 'program', 'title',
         'targname', 'prop_ra', 'prop_dec', 
         'patttype', 'patt_num', 'pattsize',
         'pridtype', 'pridtpts', 'subpxpat', 'subpxpts',
         'xoffset', 'yoffset',
         'dithdirc', 'dithopfr', 'dithpnts',
         's_region', 'crowdfld', 'pcs_mode', 'expripar', 
         'publicReleaseDate_mjd as ReleaseDate',
         'FileSetId', 'dataURI'],
  'NRS':['ArchiveFileID', 'obs_id', 'filename', 'fileSetName',
         'visit', 'visit_id', 'observtn', 'productLevel', 'apername',
         'gs_v3_pa', 'pps_aper', 'template',
         'bkglevel', 'date_mjd', 'expstart', 'expend', 'exposure','filetype',
         'telescop', 'instrume', 'detector', 'subarray',
         'nrs_norm', 'nrs_ref', 
         'filter', 'fxd_slit', 'grating', 'lamp',
         'msaconid', 'msametfl', 'msametid', 'msastate', 'preimage',
         'readpatt', 'nexposur', 'nframes', 'ngroups', 'tframe',
         'nints', 'nresets', 'effexptm as exptime',
         'pi_name', 'program', 'title',
         'targname', 'prop_ra', 'prop_dec', 
         'patttype', 'patt_num', 'pattsize',
         'pridtype', 'pridtpts', 'subpxpat', 'subpxpts', 'nod_type', 
         'xoffset', 'yoffset',
         's_region', 'crowdfld', 'pcs_mode', 'expripar', 
         'publicReleaseDate_mjd as ReleaseDate',
         'FileSetId', 'dataURI'],
}

REQUEST = {'service': 'Mast.Jwst.Filtered.Niriss',
                'params': {'columns': '*', 'filters': []},
                'format': 'json',
                'pagesize': 50000}
                
SERVICES = {'NIS':'Mast.Jwst.Filtered.Niriss', 
            'NRC':'Mast.Jwst.Filtered.Nircam',
            'MIR':'Mast.Jwst.Filtered.Miri',
            'NRS':'Mast.Jwst.Filtered.Nirspec'}

FILTER_EXAMPLES = [{'paramName':
                    'productLevel', 'values': ['2a','2b']},
                   {'paramName': 'filetype', 
                    'values': ['countrate','calibrated']},
                   {'paramName': 'program', 'values': ['1076']},
                   {'paramName': 'expstart', 
                    'values': [{"min":59672.0, "max":59690.}]},
                   {'paramName': 'filename',
                    'values': [], 'freeText':"%rate.fits"} 
                   ]

CALIB_FILTERS = [{'paramName':
                    'productLevel', 'values': ['2a','2b']},
                 {'paramName': 'filetype', 
                    'values': ['countrate','calibrated']},
                   ]

FULL_SUBARRAY = [{'paramName': 'subarray', 'values': ['FULL']}]

FINE_GUIDE = [{'paramName': 'pcs_mode', 'values': ['FINEGUIDE']}]

def programs(progs):
    """
    Helper for generating a filter on a list of program ids
    """
    return [{'paramName': 'program', 'values': [f'{int(p)}' for p in progs]}]


def query_jwst(instrument='NIS', columns='*', filters=CALIB_FILTERS+FULL_SUBARRAY+FINE_GUIDE, extra={'format':'json', 'pagesize': 100000}, recent_days=None, rates_and_cals=False, extensions=['rate', 'cal']):
    """
    Query 
    """
    import json
    import astropy.time
    
    from mastquery.utils import mastQuery, mastJson2Table
        
    if recent_days is not None:
        now = astropy.time.Time.now().mjd
        filters = [f for f in filters]
        filters.append({'paramName': 'expstart', 
                        'values': [{"min":now-recent_days, "max":now+1}]})
    
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
    
    head, content = mastQuery(request)
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
            print(f'Warning: no files with extensions {extensions} found.')
            print(f'Available extensions are {np.unique(ext).tolist()}')
            
        tab = tab[in_ext]
        
    if (recent_days is not None) & ('expstart' in tab.colnames):
        tab['dt'] = now - tab['expstart']
        tab['dt'].format = '.2f'
        
    return tab


