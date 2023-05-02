"""
Scripts to find overlapping HST data
"""

import time
import os
import yaml
import copy
import traceback
import inspect
import glob
from collections import OrderedDict

import numpy as np
from tqdm import tqdm

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from astropy.table import Table
from astropy import units as u
from astropy.coordinates import angles

from shapely.geometry import Polygon, Point

from sregion import SRegion, patch_from_polygon

from . import query, utils
from .plot_utils import draw_axis_labels, insert_legacysurveys_thumbnail

TQDM_MIN = 2000

def parse_overlap_polygons(polygons, fractional_overlap=0, verbose=2):
    """
    """
    
    match_poly = [polygons[0]]
    match_ids = [[0]]
    
    # Loop through polygons and combine those that overlap
    _iter = range(1,len(polygons))
    if len(polygons) > TQDM_MIN:
        _iter = tqdm(_iter)
        
    for i in _iter:
        # if verbose > 1:
        #     print(utils.NO_NEWLINE+'Parse {0:4d} (N={1})'.format(i, 
        #                                             len(match_poly)))
        
        has_match = False
        for j in range(len(match_poly)):
            isect = match_poly[j].intersection(polygons[i])
            #print(tab['target'][i], i, isect.area > 0)
            if fractional_overlap > 0:
                test_area = fractional_overlap*polygons[i].area
            else:
                test_area = 0.5/3600.
                
            if isect.area > test_area:
                #print(j, isect.area*3600)
                match_poly[j] = match_poly[j].union(polygons[i])
                match_ids[j].append(i)
                has_match = True
                continue
                
        if not has_match:
            match_poly.append(polygons[i])
            match_ids.append([i])
    
    ##################
    # Iterate joining polygons
    for iter in range(3):
        mpolygons = copy.deepcopy(match_poly)
        mids = copy.deepcopy(match_ids)
    
        match_poly = [mpolygons[0]]
        match_ids = [mids[0]]
        
        _iter = range(1,len(mpolygons))
        if len(mpolygons) > TQDM_MIN:
            _iter = tqdm(_iter)
        
        for i in _iter:
            # if verbose > 1:
            #     print(utils.NO_NEWLINE+'Parse, iter {0}, {1:4d}'.format(iter+1, i))
                
            has_match = False
            for j in range(len(match_poly)):
                isect = match_poly[j].intersection(mpolygons[i])
                #print(tab['target'][i], i, isect.area > 0)
                if isect.area > fractional_overlap*mpolygons[i].area:
                    #print(isect.area*3600)
                    match_poly[j] = match_poly[j].union(mpolygons[i])
                    match_ids[j].extend(mids[i])
                    has_match = True
                    continue

            if not has_match:
                match_poly.append(mpolygons[i])
                match_ids.append(mids[i])
        
        if verbose > 0:
            print('Iter #{0}, N_Patch = {1}'.format(iter+1, len(match_poly)))
        
        if len(mpolygons) == len(match_poly):
            break
    
    #np.save('overlaps.npy', [match_poly, match_ids])
    return match_poly, match_ids


def _save_poly_file(poly_file, match_poly, match_ids):
    """
    Save polygon summary
    """
    _data = {'match_poly': [SRegion(p).s_region for p in match_poly],
             'match_ids': match_ids}
             
    with open(poly_file,'w') as _fp:
        yaml.dump(_data, _fp)


def _load_poly_file(poly_file):
    """
    Load the saved polygon summary file
    """
    with open(poly_file) as _fp:
        _ = yaml.load(_fp, Loader=yaml.SafeLoader)
    
    match_poly = [SRegion(p).union(as_polygon=True) for p in _['match_poly']]
    match_ids = _['match_ids']
    return match_poly, match_ids


def find_overlaps(tab, buffer_arcmin=1., filters=[], instruments=['WFC3/IR', 'WFC3/UVIS', 'ACS/WFC'], proposal_id=[], SKIP=False, base_query=query.DEFAULT_QUERY, extra={}, close=True, use_parent=False, extensions=['FLT','C1M'], include_subarrays=False, min_area=0.2, show_parent=True, show_parent_box=True, targstr='j{rah}{ram}{ras}{sign}{ded}{dem}', prefix='', suffix='', jstr='{prefix}{jname}{suffix}', fractional_overlap=0, patch_alpha=0.1, parent_alpha=0.1, tile_alpha=0.1, verbose=2, keep_single_name=True, poly_file='overlaps.yaml', load_poly_file=False, min_radius=2, bad_area=0.3):
    """
    Compute discrete groups from the parent table and find overlapping
    datasets.
    
    Parameters
    ----------
    tab : `~astropy.table.Table`
        Parent table from which to compute the groups.
        
    buffer_arcmin : float
        Buffer, in arcminutes, to add around the parent group polygons
        
    filters : list
        List of filters to query.  If empty then return all.
        
    instruments : list
        List of instruments to query.  If empty then return all.
    
    proposid : list
        List of proposal IDs to query.  If empty then return all.
    
    SKIP : bool
        Don't recompute if a table with the same rootname already exists.
        
    extra : list
        Extra query parameters.  The available columns are described at 
        https://mast.stsci.edu/api/v0/_c_a_o_mfields.html.
        
    close : bool
        If true, close the figure objects.
    
    use_parent : bool
        Use parent table rather than performing a new query
    
    targstr, prefix, suffix, jstr : str
        The rootname is determined from a combination of these strings, 
        with `targstr` passed to `~mastquery.utils.radec_to_targname` to 
        return a string `jname` with the final root produced with `jstr` with 
        format arguments `prefix`, `suffix` and `jname.
    
    keep_single_name : bool
        If only a single polygon is found, try to use the RA/DEC keys of the 
        parent table metadata as the output key.
        
    Returns
    -------
    tables : list
        
        List of grouped tables (`~astropy.table.Table`).

    """
    # Get shapely polygons for each exposures
    polygons = []
    
    poly_buffer = buffer_arcmin/60 # ~1 arcmin, but doesn't account for cos(dec)
    #poly_buffer = 0.5/60 # ~1 arcmin, but doesn't account for cos(dec)
    
    if isinstance(tab, Polygon):
        poly_query = True
        polygons = [tab]
        use_parent = False
        show_parent = True
    else:
        #print('Parse polygons')
        poly_query = False
        
        # Fix "CLEAR" filters
        for i, filt_i in enumerate(tab['filter']):
            if 'clear' in filt_i.lower():
                spl = filt_i.lower().split(';')
                if len(spl) > 1:
                    for s in spl:
                        if 'clear' not in s:
                            #print(filt_i, s)
                            filt_i = s.upper()
                            break

                tab['filter'][i] = filt_i.upper()
        
        for i in range(len(tab)):
                        
            pshape, is_bad, poly = query.instrument_polygon(tab[i])
            
            # poly = query.parse_polygons(tab['footprint'][i])#[0]
            # pt = Point(tab['ra'][i], tab['dec'][i])
            # pt_buff = pt.buffer(20./60)
            # 
            # pshape = None
            # for pi in poly:
            #     try:
            #         psh = Polygon(pi)
            #     except:
            #         continue
            #     
            #     if psh.intersection(pt_buff).area == 0:
            #         continue
            #            
            #     if pshape is None:
            #         pshape = psh
            #     else:
            #         pshape = pshape.union(psh)
            #     
            #     # Some guide star problems
            #     if pshape.area*3600 > 20:
            #         continue
            # 
            # # Expected area, neglects subarrays
            # if tab[i]['instrument_name'] in utils.INSTRUMENT_AREAS:
            #     area = utils.INSTRUMENT_AREAS[tab[i]['instrument_name']]
            # else:
            #     area = 8.
            # 
            # msg = f"Footprint problem: i={i} {tab[i]['obs_id']}, area={area:4.1f}, {tab[i]['instrument_name']:>10}, npoly={len(poly)}, expstart={tab[i]['expstart']}"
            # 
            # # Got to the end and pshape is None, probably because doesn't 
            # # overlap with the target position 
            # # (a posteriori alignment problems?)
            # if pshape is None:
            #     print(msg)
            #     pshape = pt.buffer(np.sqrt(area)*np.sqrt(2)/60.)
            # else:
            #     if pshape.area < 0.1*area/3600:
            #         print(msg)
            #         pshape = pt.buffer(np.sqrt(area)*np.sqrt(2)/60.)
                
            polygons.append(pshape.buffer(poly_buffer))
            
    # match_poly = [polygons[0]]
    # match_ids = [[0]]
    # 
    # # Loop through polygons and combine those that overlap
    # for i in range(1,len(polygons)):
    #     if verbose > 1:
    #         print(utils.NO_NEWLINE+'Parse {0:4d}'.format(i))
    #     
    #     has_match = False
    #     for j in range(len(match_poly)):
    #         isect = match_poly[j].intersection(polygons[i])
    #         #print(tab['target'][i], i, isect.area > 0)
    #         if fractional_overlap > 0:
    #             test_area = fractional_overlap*polygons[i].area
    #         else:
    #             test_area = 0.5/3600.
    #             
    #         if isect.area > test_area:
    #             #print(isect.area*3600)
    #             match_poly[j] = match_poly[j].union(polygons[i])
    #             match_ids[j].append(i)
    #             has_match = True
    #             continue
    #             
    #     if not has_match:
    #         match_poly.append(polygons[i])
    #         match_ids.append([i])
    # 
    # ##################
    # # Iterate joining polygons
    # for iter in range(3):
    #     mpolygons = copy.deepcopy(match_poly)
    #     mids = copy.deepcopy(match_ids)
    # 
    #     match_poly = [mpolygons[0]]
    #     match_ids = [mids[0]]
    # 
    #     for i in range(1,len(mpolygons)):
    #         if verbose > 1:
    #             print(utils.NO_NEWLINE+'Parse, iter {0}, {1:4d}'.format(iter+1, i))
    #             
    #         has_match = False
    #         for j in range(len(match_poly)):
    #             isect = match_poly[j].intersection(mpolygons[i])
    #             #print(tab['target'][i], i, isect.area > 0)
    #             if isect.area > fractional_overlap*mpolygons[i].area:
    #                 #print(isect.area*3600)
    #                 match_poly[j] = match_poly[j].union(mpolygons[i])
    #                 match_ids[j].extend(mids[i])
    #                 has_match = True
    #                 continue
    # 
    #         if not has_match:
    #             match_poly.append(mpolygons[i])
    #             match_ids.append(mids[i])
    #     
    #     if verbose > 0:
    #         print('Iter #{0}, N_Patch = {1}'.format(iter+1, len(match_poly)))
    #     
    #     if len(mpolygons) == len(match_poly):
    #         break
            
    if os.path.exists(poly_file) & (load_poly_file):
        match_poly, match_ids = _load_poly_file(poly_file)
    else:
        match_poly, match_ids = parse_overlap_polygons(polygons, 
                                      fractional_overlap=fractional_overlap, 
                                      verbose=verbose)
        
        _save_poly_file(poly_file, match_poly, match_ids)                           
        #np.save(poly_file, [match_poly, match_ids])

    # Save figures and tables for the unique positions
    BLUE = '#6699cc'
    
    tables = []
    
    NPOLYS = len(match_poly)
    
    for ipo in range(NPOLYS):
        #i+=1
        p = match_poly[ipo].buffer(0.0001)
                
        #######
        # Query around central coordinate
        idx = np.array(match_ids[ipo])
        if poly_query:
            ra = float(p.centroid.x)
            dec = float(p.centroid.y)
        else:
            if (NPOLYS == 1) & keep_single_name:
                try:
                    ra, dec = tab.meta['RA'], tab.meta['DEC']
                except:
                    ra = np.mean(tab['ra'][idx])
                    dec = np.mean(tab['dec'][idx])                
            else:
                ra, dec = np.mean(tab['ra'][idx]), np.mean(tab['dec'][idx])
        
        # Get poly size
        xy = np.array(p.convex_hull.boundary.xy)
        xradius = np.abs(xy[0]-ra).max()*np.cos(dec/180*np.pi)*60
        yradius = np.abs(xy[1]-dec).max()*60
        
        box_radius = np.max([xradius*1.5, yradius*1.5, min_radius])
        box = [ra, dec, box_radius]
                
        # Build target name from RA/Dec
        jname = utils.radec_to_targname(box[0], box[1], round_arcsec=(4, 60), targstr=targstr)
        jname = jstr.format(prefix=prefix, jname=jname, suffix=suffix) #prefix+jname+suffix
        
        if (os.path.exists('{0}_footprint.pdf'.format(jname))) & SKIP:
            print('\n********** SKIP *********\n', ipo+1, jname, box[0], box[1])
            continue
        else:
            print('\n\n', ipo+1, jname, box[0], box[1])
            
        if use_parent:
            xtab = tab
        else:
            xtab = query.run_query(box=box, proposal_id=proposal_id, instruments=instruments, filters=filters, base_query=base_query, **extra)
            if isinstance(xtab, dict):
                utils.log_comment(jname+'.failed', xtab, verbose=False, 
                                  show_date=True, mode='a')
                utils.log_exception(jname+'.failed', traceback, 
                                    verbose=True, mode='a')
                continue
                
            if not include_subarrays:
                query.set_area_column(xtab)
                try:
                    subarray = (xtab['instrument_name'] == 'WFC3/IR') & (xtab['area'] < 3.8)
                except:
                    # Bug?
                    xtab['instrument_name'].fill_value = 0
                    subarray = (xtab['instrument_name'] == 'WFC3/IR') & (xtab['area'] < 3.8)
                    
                xtab = xtab[~subarray]
            
        try:
            #ebv = utils.get_irsa_dust(ra, dec, type='SandF')
            ebv = utils.get_mw_dust(ra, dec)
        except:
            print('Couldn\'t connect to IRSA dust server, setting MW_EBV=0')
            ebv = 0
            
        xtab.meta['NAME'] = jname
        xtab.meta['RA'] = ra
        xtab.meta['DEC'] = dec
        xtab.meta['MW_EBV'] = ebv
        xtab.meta['FOLAP'] = (fractional_overlap, 'Fractional overlap parameter')
        xtab.meta['MIN_AREA'] = (min_area, 'Minimum overlap fraction')
        xtab.meta['BUFFER'] = (buffer_arcmin, 'Buffer for overlaps in arcmin')
        xtab.meta['SREGION'] = (utils.polygon_to_sregion(p), 'Parent footprint')
        
        # Only include ancillary data that directly overlaps with the primary
        # polygon
        pointing_overlaps = np.zeros(len(xtab), dtype=bool)
        for j in range(len(xtab)):            
            pshape, is_bad, poly = query.instrument_polygon(xtab[j])
            # if is_bad:
            #     print('xxxx', pshape.centroid.xy, pshape.area*3600., p.centroid.xy, p.area*3600)
                
            # try:
            #     poly = query.parse_polygons(xtab['footprint'][j])#[0]
            # except:
            #     pointing_overlaps[j] = False
            #     continue
            
            # pshape = None
            # for pi in poly:
            #     try:
            #         psh = Polygon(pi)
            #     except:
            #         continue
            #     
            #     if psh.intersection(pt_buff).area == 0:
            #         continue
            #         
            #     if pshape is None:
            #         pshape = psh
            #     else:
            #         pshape = pshape.union(psh)
            # 
            # # Expected area, neglects subarrays
            # if xtab[j]['instrument_name'] in utils.INSTRUMENT_AREAS:
            #     area = utils.INSTRUMENT_AREAS[xtab[j]['instrument_name']]
            # else:
            #     area = 8.
            # 
            # msg = f"Footprint problem: i={j} {xtab[j]['obs_id']}, area={area:4.1f}, {xtab[j]['instrument_name']:>10}, npoly={len(poly)}, expstart={xtab[j]['expstart']}"
            # 
            # # Got to the end and pshape is None, probably because doesn't 
            # # overlap with the target position 
            # # (a posteriori alignment problems?)
            # if pshape is None:
            #     print(msg)
            #     pshape = pt.buffer(np.sqrt(area)*np.sqrt(2)/60.)
            # else:
            #     if pshape.area < 0.1*area/3600:
            #         print(msg)
            #         pshape = pt.buffer(np.sqrt(area)*np.sqrt(2)/60.)
            
            isect = p.intersection(pshape.buffer(0.0001))
            pointing_overlaps[j] = isect.area > min_area*pshape.area
                                    
            try:
                isect = p.intersection(pshape.buffer(0.0001))
                pointing_overlaps[j] = isect.area > min_area*pshape.area
            except:
                pointing_overlaps[j] = False
                
        if pointing_overlaps.sum() == 0:
            continue
        
        xtab = xtab[pointing_overlaps]
        
        # Unique targets
        filter_target = np.array(['{0} {1}'.format(xtab['instrument_name'][i], xtab['filter'][i]) for i in range(pointing_overlaps.sum())])
                        
        ########### 
        # Make the figure
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        # Show the parent table
        if show_parent:
            if poly_query:
                ax.add_patch(patch_from_polygon(p, alpha=parent_alpha))
            else:
                colors = query.show_footprints(tab[idx], ax=ax, alpha=parent_alpha, bad_area=bad_area)
                
        ax.scatter(box[0], box[1], marker='+', color='k', zorder=1e4, alpha=0.8)
        
        if not poly_query:
            if ('boxra' in tab.meta) & show_parent_box:
                ax.scatter(tab.meta['boxra'][0], tab.meta['boxdec'][0],
                           marker='*', color='r', zorder=1e4, alpha=0.8)
            
        try:
            colors = query.show_footprints(xtab, ax=ax, alpha=tile_alpha)
        except:
            continue
            
        if patch_alpha > 0:
            patch1 = patch_from_polygon(p, fc=BLUE, ec=BLUE,
                                  alpha=patch_alpha, zorder=2)
        
            ax.plot(xy[0], xy[1], alpha=patch_alpha, color=BLUE)
            ax.add_patch(patch1)
        
        ax.grid()
        
        # Resize for square dimensions
        xr, yr = ax.get_xlim(), ax.get_ylim()
        dx = (xr[1]-xr[0])*np.cos(yr[0]/180*np.pi)*60
        dy = (yr[1]-yr[0])*60
        ax.set_title(jname)
        ax.set_xlim(ax.get_xlim()[::-1])
        
        cosd = np.cos(dec/180*np.pi)
        ax.set_aspect(1/cosd)
        
        fig.set_size_inches(5,5*np.clip(dy/dx, 0.2, 5))
        
        # Add summary
        fp = open('{0}_info.dat'.format(jname),'w')
        dyi = 0.02*np.clip(dx/dy, 0.2, 5)
        
        title = '{0:>13.5f} {1:>13.5f}  E(B-V)={2:.3f}'.format(ra, dec, ebv)
        fp.write('# {0}\n'.format(title))
        
        ax.text(0.05, 0.97, title, ha='left', va='top', transform=ax.transAxes, fontsize=6)
        
        for i, t in enumerate(np.unique(xtab['proposal_id'])):
            fp.write('proposal_id {0} {1}\n'.format(jname, t))
            ts = np.unique(xtab['target'][xtab['proposal_id'] == t])
            
            if len(ts) > 4:
                tstr = '{0} {1}'.format(t, ' '.join(['{0}'.format(ti) for ti in ts[:4]])) + ' ...'
            else:
                tstr = '{0} {1}'.format(t, ' '.join(['{0}'.format(ti) for ti in ts]))
                
            ax.text(0.05, 0.97-dyi*(i+1), tstr, ha='left', va='top', transform=ax.transAxes, fontsize=6)
            
        for i, t in enumerate(np.unique(xtab['target'])):
            fp.write('target {0} {1}\n'.format(jname, t))
            
        print(np.unique(xtab['target']), '\n')
        
        for i, filt in enumerate(np.unique(filter_target)):
            mf = filter_target == filt
            print('filter {0}  {1:>20s}  {2:>3d}  {3:>8.1f}'.format(jname, filt, mf.sum(), xtab['exptime'][mf].sum()))
            fp.write('filter {0}  {1:>20s}  {2:>3d}  {3:>8.1f}\n'.format(jname, filt, mf.sum(), xtab['exptime'][mf].sum()))
            
            c = colors[filt.split()[1]]
            ax.text(0.95, 0.97-dyi*i, '{1:>20s}  {2:>3d}  {3:>8.1f}\n'.format(jname, filt, mf.sum(), xtab['exptime'][mf].sum()), ha='right', va='top', transform=ax.transAxes, fontsize=6, color=c)
            
        fp.close()
        
        # Timestamp        
        ax.text(0.97, 0.03, time.ctime(), fontsize=5, transform=ax.transAxes, ha='right', va='bottom')
        
        fig.tight_layout(pad=0.5)
        
        # Save figure and table
        fig.savefig('{0}_footprint.pdf'.format(jname))
        
        if close:
            plt.close()

        xtab.write('{0}_footprint.fits'.format(jname), format='fits', 
                   overwrite=True)
        
        _data = {'p':SRegion(p).s_region, 
                 'box':[float(b) for b in box]}
        
        with open('{0}_footprint.yaml'.format(jname),'w') as _fp:
            yaml.dump(_data, _fp)
            
        #np.save('{0}_footprint.npy'.format(jname), [p, box])
        
        tables.append(xtab)
    
    return tables


def summary_table(tabs=None, output='overlap_summary'):

    from astropy.table import Table
    import astropy.table
    from mastquery.overlaps import parse_overlap_table
    from mastquery.query import set_transformed_coordinates
    
    try:
        from grizli import utils
        HAS_GRIZLI = True
    except:
        HAS_GRIZLI = False
        
    if tabs is None:
        tabs = [Table.read(file) for file in glob.glob('*footprint.fits')]
    
    for tab in tabs:
        set_transformed_coordinates(tab)
        
    # for tab in tabs:
    #     tab.remove_column('moving_target')
    # #
    # for tab in tabs:
    #     prop = np.cast[int](tab['proposal_id'])
    #     tab.remove_column('proposal_id')
    #     tab['proposal_id'] = prop
        
    names, props = parse_overlap_table(tabs[0])
    props = []
    pdict = OrderedDict()
    for name in names:
        pdict[name] = []
        
    for i in range(len(tabs)):
        print('Parse table ', i)
        n_i, p_i = parse_overlap_table(tabs[i])
        for name in names:
            if name in n_i:
                pdict[name].append(p_i[n_i.index(name)])
            else:
                pdict[name].append('---')
                          
    mtab = Table(pdict) #rows=props, names=names)
    mtab['RA'].format = '.5f'
    mtab['DEC'].format = '.5f'
    #mtab.rename_column('MW_EBV', 'E(B-V)')
    mtab['MW_EBV'].format = '.2f'
    
    mtab['EclLat'].format = '.1f'
    mtab['EclLon'].format = '.1f'
    
    mtab['GalLat'].format = '.1f'
    mtab['GalLon'].format = '.1f'
    
    mtab['AreaG102'].format = '.1f'
    mtab['AreaG141'].format = '.1f'
    
    mtab['TexpG102'].format = '.1f'
    mtab['TexpG141'].format = '.1f'
    
    mtab['TperG102'].format = '.1f'
    mtab['TperG141'].format = '.1f'
    
    # Links
    mast_link = ['<a href=https://archive.stsci.edu/hst/search.php?RA={0}&DEC={1}&radius=3.&max_records=5000&sci_aec=S&action=Search>{2}</a>'.format(t['RA'], t['DEC'], 'MAST') for t in mtab]
    mtab['MAST'] = mast_link
    
    fp_link = ['<a href={0}_footprint.pdf>Footprint</a>'.format(t['NAME']) for t in mtab]
    mtab['Footprint'] = fp_link

    fp_link = ['<a href={0}_footprint.pdf><img src={0}_footprint.png height=300px></a>'.format(t['NAME']) for t in mtab]
    mtab['Footprint'] = fp_link
    
    mtab.write('{0}.fits'.format(output), overwrite=True)
    
    if HAS_GRIZLI:
        gtab = utils.GTable(mtab)
        gtab.write_sortable_html('{0}.html'.format(output),
                     replace_braces=True, localhost=False, 
                     max_lines=len(mtab)+10, table_id=None, 
                     table_class='display compact', css=None, 
                     filter_columns=['NG102', 'NG141', 'NG800L', 'NG280', 'MW_EBV'])
    
        return gtab
    else:
        return mtab


def parse_overlap_table(tab):
    """
    Compute properties of the overlap table
    
    Parameters
    ----------
    tab : `~astropy.table.Table`
        Overlap table output from `find_overlaps`.
        
    Returns
    -------
    names : list
        If `get_colnames==True`, then also return a list of the parameter
        names.

    properties : list
        List of extracted properties.
        
    """
    from . import query# as mutils
    
    query.set_transformed_coordinates(tab)
    
    # Meta attributes
    names, properties = [], []
    for k in tab.meta:
        names.append(k)
        properties.append(tab.meta[k])
    
    # Ecliptic and Galactic coords
    names.extend(['EclLat','EclLon','GalLat','GalLon'])
    properties.append(np.mean(tab['ecl_lat']))
    properties.append(np.mean(tab['ecl_lon']))
    properties.append(np.mean(tab['gal_b']))
    properties.append(np.mean(tab['gal_l']))
    
    # Unique elements of the table
    names.append('NFILT')
    properties.append(len(np.unique(tab['filter'])))
    for c in ['filter', 'target', 'target_classification', 'proposal_id']:
        if c not in tab.colnames:
            continue
        
        names.append(c)
        properties.append(' '.join(['{0}'.format(p) for p in np.unique(tab[c])]))
        if c in ['target_classification']:
            properties[-1] = properties[-1].title().replace(';','\n')
            
    for c in ['proposal_pi']:
        names.append(c.split()[0])
        properties.append(' '.join(['{0}'.format(p.split()[0].title()) for p in np.unique(tab[c])]))
    
    # By grism
    for g in ['G102', 'G141', 'G800L', 'G280']:
        m = tab['filter'] == g
        names.extend(['{0}{1}'.format(p, g) for p in ['N', 'Area', 'Texp', 'Tper', 'PA']])
        if m.sum() == 0:
            properties.extend([0,0,0,0,0])
            continue
        
        # N
        properties.append(m.sum())
        
        # Area
        PAs = []
        for i, polystr in enumerate(tab['footprint'][m]):
            try:
                poly = query.parse_polygons(polystr)
            except:
                poly = query.old_parse_polygons(polystr)
                
            PAs.append(int(np.round(query.get_orientat(polystr))))

            pshape = None
            for pi in poly:
                try:
                    psh = Polygon(pi)
                except:
                    continue
                    
                if pshape is None:
                    pshape = psh
                else:
                    pshape = pshape.union(psh)
            
            if i == 0:
                gpoly = pshape
            else:
                try:
                    gpoly = gpoly.union(pshape)
                except:
                    print(f'Poly failed ({i})')
                    
        cosd = np.cos(tab.meta['DEC']/180*np.pi)
        area = gpoly.area*3600.*cosd
        properties.append(area)
        
        # Texp
        texp = np.sum(tab['exptime'][m])
        properties.append(texp)
        
        # Tper
        properties.append(texp/3000./(area/4.4))
        
        # Number of PAs
        properties.append(len(np.unique(PAs)))

    return names, properties


ASSOC_ARGS = {'max_pa':2, 'max_sep':0.5, 'max_time':1.e4/86400., 'match_filter':True, 'match_instrument':True, 'match_program':True, 'hack_grism_pa':True, 'parse_for_grisms':True}


def split_associations(tab, force_split=False, root=None, assoc_args=ASSOC_ARGS, make_figure=True, xsize=6, nlabel=3, assoc_min=0, fill_grism=True, force_fill=False, **kwargs):
    """
    Split table by groups from `compute_associations`
    
    assoc_args passed directly to `compute_associations`.
    
    """ 
    if root is None:
        root = tab.meta['NAME']
    
    if ('assoc_idx' not in tab.colnames) | (force_split):
        compute_associations(tab, **assoc_args)
        tab['assoc_idx'] += assoc_min
        
    polys = OrderedDict()
    
    tab_poly = None
    
    for ix in np.unique(tab['assoc_idx']):
        sel = tab['assoc_idx'] == ix
        #orient = np.round(np.median(tab['orientat'][sel])/round_pa)*round_pa
        
        prog = tab['obs_id'][sel][0][1:4]
        visit = tab['obs_id'][sel][0][4:6]

        targ = '-'.join(np.unique(tab['target'][sel])).lower()
        filt = '-'.join(np.unique(tab['filter'][sel])).lower()
        inst = '-'.join(np.unique(tab['instrument_name'][sel])).lower()
        inst = inst.replace('/','')
        
        label = '{0}_{5:04d}_{1}_{2}_{3}_{4}'
        label = label.format(root, prog, targ, inst, filt, ix)
        
        poly_i = None
        
        for j in range(sel.sum()):
            p_j, is_bad, poly = query.instrument_polygon(tab[sel][j])
                
            if tab_poly is None:
                tab_poly = p_j.buffer(0.001)
            
            if not tab_poly.buffer(2).intersects(p_j.buffer(2)):
                print('Skip')
                continue
                
            if poly_i is None:
                poly_i = p_j
            else:
                if not poly_i.buffer(2).intersects(p_j.buffer(2)):
                    print('x Skip')
                    continue

                poly_i = poly_i.union(p_j)
        
        # for fp in tab['footprint'][sel]:
        #     for p in query.parse_polygons(fp):
        #         p_j = Polygon(p).buffer(0.001)
        #         if tab_poly is None:
        #             tab_poly = p_j.buffer(0.001)
        #         
        #         if not tab_poly.buffer(2).intersects(p_j.buffer(2)):
        #             print('Skip')
        #             continue
        #             
        #         if poly_i is None:
        #             poly_i = Polygon(p)
        #         else:
        #             if not poly_i.buffer(2).intersects(p_j.buffer(2)):
        #                 print('x Skip')
        #                 continue
        # 
        #             poly_i = poly_i.union(p_j)
                        
        polys[label] = poly_i
    
    if make_figure:
        fig = make_association_figure(tab, polys, root=root, xsize=xsize, nlabel=nlabel, fill_grism=fill_grism, force_fill=force_fill, **kwargs)
        return polys, fig
    else:
        return polys

LS_ARGS = dict(pixscale=1,
               layers=['ls-dr9', 'sdss', 'unwise-neo7'],
               zorder=-1000,
               alpha=0.8,
               aspect='auto',
               verbose=True,
               grayscale=False,
               grayscale_params=[99, 1.5, -0.02])

def make_association_figure(tab, polys, highlight=None, root=None, xsize=6, nlabel=3, fill_grism=True, force_fill=False, with_ls_thumbnail=False, ls_args=LS_ARGS, db_query_str=None, pad_arcmin=1, **kwargs):
    """Make a figure to show associations
    
    with_ls_thumbnail : bool
        Get a cutout from LegacySurveys
    
    ls_args : dict
        Arguments to 
        
    """    
    #cc = plt.rcParams['axes.prop_cycle']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    
    if root is None:
        root = tab.meta['NAME']
    
    handles, labels = [], []
    
    fig = plt.figure()
    ax = fig.add_subplot(111) #, projection=pywcs.WCS())

    sc = ax.scatter(tab.meta['RA'], tab.meta['DEC'], marker='+', c='k')
    title = '{0}  ({1:.4f}, {2:.4f})  E(B-V)={3:.03f}'.format(root, tab.meta['RA'], tab.meta['DEC'], tab.meta['MW_EBV'])
    ax.set_title(title)
    
    #labels.append(label)
    
    filter_list = []
    full_filter_list = [f.split('_')[-1].split('-')[-1] for f in polys]
    so = np.argsort(full_filter_list)
    keys = list(polys.keys())
    
    grism_patches = {}
    grism_colors = {'g141':'r', 'g102':'orange', 'g800l':'b'}
    
    # Fix "CLEAR" filters
    for i, filt_i in enumerate(tab['filter']):
        if 'clear' in filt_i.lower():
            spl = filt_i.lower().split(';')
            if len(spl) > 1:
                for s in spl:
                    if 'clear' not in s:
                        #print(filt_i, s)
                        filt_i = s.upper()
                        break
            
            tab['filter'][i] = filt_i.upper()
        
    xra, yra = None, None
    for i in so:
        # in polys:
        f = keys[i]
        p_i = polys[f]
        xfilt_i = f.split('_')[-1]
        filt_i = xfilt_i.split('-')[-1]
        if 'clear' in filt_i.lower():
            spl = filt_i.split(';')
            if len(spl) > 1:
                for s in spl:
                    if 'clear' not in s:
                        filt_i = s
                        break
                                
        if hasattr(p_i, 'geoms'):
            all_p = [p for p in p_i.geoms]
        else:
            all_p = [p_i]
        
        for p in all_p:
            try:
                xb, yb = np.array(p.buffer(0.001).boundary.xy)
            except:
                try:
                    xb, yb = np.array(p.convex_hull.boundary.xy)
                except:
                    continue
                
            if xra is None:
                xra = [xb.min(), xb.max()]
                yra = [yb.min(), yb.max()]
            else:
                xra = [np.minimum(xb.min(), xra[0]), 
                       np.maximum(xb.max(), xra[1])]
                   
                yra = [np.minimum(yb.min(), yra[0]), 
                       np.maximum(yb.max(), yra[1])]

            ax.set_xlim(xra[::-1])
            ax.set_ylim(yra)
        
            if 'g141' in xfilt_i:
                has_grism = True
            elif 'g102' in xfilt_i:
                has_grism = True
            elif 'g800l' in xfilt_i:
                has_grism = True    
            elif 'gr150' in xfilt_i:
                has_grism = True
                c_i = colors[hash(xfilt_i) % len(colors)]
                grism_colors[xfilt_i] = c_i
                #print('!!! ', xfilt_i, c_i, grism_colors)
            else:
                c_i = colors[hash(filt_i) % len(colors)]
                has_grism = force_fill
                grism_colors[filt_i] = c_i
                has_grism = True
            
            if filt_i in grism_colors:
                c_i = grism_colors[filt_i]
                
            fc = 'None'
            zo = -1000
            if highlight is not None:
                if f == highlight:
                    fc = c_i
                    zo = 1000
            
            alpha = 0.4**(not with_ls_thumbnail)
            
            if (fill_grism & has_grism) | (force_fill):                
                #fc = c_i
                #zo = -100000
                #alpha = 0.25
                if filt_i in grism_patches:
                    grism_patches[filt_i] = grism_patches[filt_i].union(p.buffer(0.00001))
                else:
                    grism_patches[filt_i] = p.buffer(0.00001)
                
                if hasattr(p, 'geoms'):
                    piter = p.geoms
                else:
                    piter = [p]
                
                for pi in piter:
                    patch = patch_from_polygon(pi,
                                         alpha=0.2**(not with_ls_thumbnail)/2,
                                         fc=fc, ec=c_i,
                                         label=None, zorder=zo)
            
                    pat = ax.add_patch(patch)
                
            else:        
                if filt_i not in filter_list:
                    _fmatch = tab['filter'] == filt_i.upper()
                    _expt = tab['exptime'][_fmatch].sum()/1000
                    if filt_i.lower().startswith('nc.'):
                        # NIRCAM chips
                        if filt_i.lower() < 'nc.f260':
                            _expt /= 8
                        else:
                            _expt /= 2
                            
                    label = '{0:5} {1:>8.1f}'.format(filt_i.strip(';'), _expt)
                else:
                    label = '_'.join(f.split('_')[1:])
                
                if hasattr(p, 'geoms'):
                    piter = p.geoms
                else:
                    piter = [p]
                
                for pi in piter:
                    patch = patch_from_polygon(pi, alpha=alpha, fc=fc, ec=c_i,
                                     label=label, zorder=zo)
            
                    pat = ax.add_patch(patch)
                             
                if filt_i not in filter_list:
                    handles.append(pat)
                    labels.append(label)
                    filter_list.append(filt_i)
    
    # print('xxx', grism_patches)
    
    if fill_grism | force_fill:
        for g in grism_patches:
            fc = ec = grism_colors[g]
            filt_i = g
            
            is_grism = filt_i in ['g102', 'g141', 'g800l']
            is_grism |= 'gr150' in filt_i
            
            if is_grism | force_fill:
                alpha = 0.2**(not with_ls_thumbnail)
            else:
                fc = 'None'
                alpha = 0.2**(not with_ls_thumbnail)
            
            _expt = tab['exptime'][tab['filter'] == filt_i.upper()].sum()/1000
            if filt_i.lower().startswith('nc.'):
                # NIRCAM chips
                if filt_i.lower() < 'nc.f260':
                    _expt /= 8
                else:
                    _expt /= 2
            
            label = '{0:5} {1:>8.1f}'.format(filt_i.strip(';'), _expt)
            
            if hasattr(grism_patches[g], 'geoms'):
                piter = grism_patches[g].geoms
            else:
                piter = [grism_patches[g]]
            
            for pi in piter:
                patch = patch_from_polygon(pi, alpha=alpha, fc=fc, ec=ec,
                                 label=label, zorder=100)
        
                pat = ax.add_patch(patch)
                         
            if filt_i not in filter_list:
                handles.append(pat)
                labels.append(label)
                filter_list.append(filt_i)
    
    try:
        from grizli.aws import db
    except:
        print('db_query_str specified but `from grizli.aws import db` failed!')
        db_query_str = None
    
    if db_query_str is not None:
        # "filter in ('F160W','F115W-CLEAR','F444W-CLEAR') OR instrume in ('NIRISS')"
        
        pt = f"point({tab.meta['RA']:.5f}, {tab.meta['DEC']:.5f})"
        sr = SRegion(tab.meta['SREGION'][0])

        R = np.sqrt(2*sr.sky_area()[0]).value
        
        _sql = f"""
SELECT file, filter, instrume, footprint, crval1, crval2
    FROM exposure_files
    WHERE polygon(footprint) && polygon(circle({pt}, {R/60.:.3f}))
    AND {db_query_str}"""
        
        dbfp = db.SQL(_sql)
        
        if len(dbfp) > 0:
            print(f'{len(dbfp)} rows found for query: {_sql}')
            
            for fp, insi, fi in zip(dbfp['footprint'], dbfp['instrume'], dbfp['filter']):
                sr = SRegion(fp)
                if insi in ('NIRCAM', 'NIRISS','MIRI'):
                    zo=-10
                    c = 'olive'
                    al = 0.2
                else:
                    if fi in ('F160W'):
                        c = 'lightsteelblue'
                    else:
                        c = 'skyblue'
                        
                    zo=-100
                    al = 0.05

                for p in sr.patch(ec='None', fc=c,alpha=al, zorder=zo):
                    ax.add_patch(p)
    
    ax.legend(handles, labels, fontsize=7, 
              ncol=int(np.minimum(len(labels), 4)), loc='upper right')
    
    ax.grid()

    cosd = np.cos(tab.meta['DEC']/180*np.pi)

    xra += pad_arcmin/60/cosd*np.array([-1,1])
    yra += pad_arcmin/60*np.array([-1,1])
    ax.set_xlim(xra[::-1])
    ax.set_ylim(yra)
    
    dy = yra[1] - yra[0]
    dx = xra[1] - xra[0]
    
    ax.set_aspect(1./cosd)
    
    tab.meta['XMIN'] = xra[0]
    tab.meta['XMAX'] = xra[1]
    tab.meta['YMIN'] = yra[0]
    tab.meta['YMAX'] = yra[1]
    # Boxquery for HSC
    tab.meta['boxquery'] = 'boxQuery(coord,{0:.5f},{1:.5f},{2:.5f},{3:.5f})'.format(xra[0], xra[1], yra[0], yra[1])

    fig.set_size_inches(xsize, xsize*np.clip(dy/(dx*cosd), 0.2, 5))
    
    draw_axis_labels(ax=ax, nlabel=nlabel)
    
    ax.text(0.03, 0.03, time.ctime(), fontsize=5,
            transform=ax.transAxes, ha='left', va='bottom')
    
    if with_ls_thumbnail:
        url, img = insert_legacysurveys_thumbnail(ax, **ls_args)
        
        if img is not None:
            ax.text(0.03, 0.03, time.ctime(), fontsize=5, color='w',
                    transform=ax.transAxes, ha='left', va='bottom')

    fig.tight_layout(pad=0.2)
    return fig

def compute_associations(tab, max_sep=0.5, max_pa=0.05, max_time=1e4/86400., match_filter=True, match_instrument=True, match_program=True, hack_grism_pa=True, parse_for_grisms=True, match_detector=True):
    """
    Associate visits by filter + position + PA + date
    """
    from . import query
    
    cosd = np.cos(tab['dec']/180*np.pi)
    dx = (tab['ra'] - np.median(tab['ra']))*cosd*60    
    dy = (tab['dec'] - np.median(tab['dec']))*60
    
    ori = np.array([query.get_orientat(tab['footprint'][i])
                    for i in range(len(tab))])
     
    if hack_grism_pa:
        ## At least one visit (Refsdal) had a 90 degree offset between 
        ## the grism and associated direct exposures computed this way
        is_grism = (tab['filter'] == 'G141') | (tab['filter'] == 'G102') | (tab['filter'] == 'G800L')
        if is_grism.sum() > 0:
            grism_progs = np.unique(tab['proposal_id'][is_grism])
            grism_prog = tab['proposal_id'] == grism_progs[0]
            for id in grism_progs:
                grism_prog |= tab['proposal_id'] == id
            
            ori[grism_prog] = ori[grism_prog] % 90

    tab['orientat'] = ori
    indices = np.arange(len(tab))
    assoc_idx = np.zeros(len(tab), dtype=int)-1
    
    if 'dataURL' in tab.colnames:
        visit_numbers = []
        for d in tab['dataURL']:
            if '/product/jw' in d:
                # JWST
                _file = d.split('/product/')[-1]
                _visit = _file.split('_')[-6:]
            else:
                _visit = d.split('/product/')[-1][4:6]
            visit_numbers.append(_visit)
            
    else:
        visit_numbers = np.array([o[4:6] for o in tab['obs_id']])
    
    assoc = []
    for i in range(len(tab)):
        if assoc_idx[i] >= 0:
            continue
                
        assoc_idx[i] = len(assoc)
        
        assoc_i = {'pos':Point(dx[i], dy[i]).buffer(max_sep),
                   'ori':ori[i],
                   'filter':tab['filter'][i],
                   'indices':[i],
                   'proposal_id':tab['proposal_id'][i],
                   'idx':len(assoc),
                   't_min':tab['t_min'][i],
                   't_max':tab['t_max'][i],
                   'instrument_name':tab['instrument_name'][i],
                   'visit_number':visit_numbers[i]}
        
        if 'detector' in tab.colnames:
            assoc_i['detector'] = tab['detector'][i]
            
        for j in range(i+1, len(tab)):
            
            if assoc_idx[j] >= 0:
                continue
                
            f_j = tab['filter'][j]
            pr_j = tab['proposal_id'][j]
            
            if 'detector' in assoc_i:
                det_j = tab['detector'][j]
                
            visit_j = visit_numbers[j]
            instr_j = tab['instrument_name'][j]
            dpa = assoc_i['ori'] - ori[j]
            dt =  tab['t_min'][j] - assoc_i['t_max']
            p_j = Point(dx[j], dy[j]).buffer(max_sep)
            
            # Has match
            test = (np.abs(dpa) < max_pa) & (p_j.intersects(assoc_i['pos']))
            test &= np.abs(dt) < max_time
            if not match_detector:
                test |= (visit_j == assoc_i['visit_number'])
            
            if match_instrument:
                test &= (instr_j == assoc_i['instrument_name'])
                
            if match_filter & (not parse_for_grisms):
                test &= (f_j == assoc_i['filter'])
            
            if match_program:
                test &= (pr_j == assoc_i['proposal_id'])
            
            if match_detector & ('detector' in assoc_i):
                test &= (det_j == assoc_i['detector'])
                
            if test:
                #print('Has match!', j)
                #break
                
                assoc_idx[j] = assoc_idx[i]
                assoc_i['pos'] = assoc_i['pos'].union(p_j)
                assoc_i['indices'].append(j)
                assoc_i['t_max'] = tab['t_max'][j]
        
        assoc.append(assoc_i)
    
    if match_filter & (parse_for_grisms):
        # Now split assoc that *don't* have grisms
        assoc_orig = assoc_idx*1
        new_i = 0
        assoc_idx = assoc_orig*0
        for ix in np.unique(assoc_orig):
            sel = assoc_orig == ix
            filts = list(np.unique(tab['filter'][sel]))
            is_grism = np.sum([f in filts
                               for f in ['G800L','G102','G141']]) > 0
            
            if is_grism:
                #print('XXX'); break
                assoc_idx[sel] = new_i
                new_i += 1
            else:
                for filt in filts:
                    sel_filt = sel & (tab['filter'] == filt)
                    assoc_idx[sel_filt] = new_i
                    new_i += 1
                
    tab['assoc_idx'] = assoc_idx


def muse_query(tab, make_figure=True, xsize=5, nlabel=3, min_size=4, cmap='jet_r', rerun_query=True, query_kwargs={'public':False, 'science':False, 'get_html_version':True}):
    """
    Query ESO archive around the HST data
    """
    import urllib

    from astropy.table import Table
    from astropy.time import Time
    from astropy.coordinates import SkyCoord
    
    from shapely.geometry import Polygon, Point
    #from descartes import PolygonPatch
    
    from astroquery.eso import Eso
    eso = Eso()
    
    NOW = Time.now().iso
    
    meta = tab.meta
    xr = (meta['XMIN'], meta['XMAX'])
    yr = (meta['YMIN'], meta['YMAX'])
    ra, dec = meta['BOXRA'], meta['BOXDEC']
    
    cosd = np.cos(dec/180*np.pi)
    dx = (xr[1]-xr[0])*cosd*60
    dy = (yr[1]-yr[0])*60
    
    box_width = np.maximum(dx, dy)
    query_size = np.maximum(min_size, box_width/2)
    
    coo = SkyCoord(ra, dec, unit='deg')
    
    muse_file = '{0}_muse.ecsv'.format(tab.meta['NAME'])
    if (not os.path.exists(muse_file)) | rerun_query:
        res = eso.query_instrument('muse', coord1=ra, coord2=dec, box="00 {0:02} 00".format(int(np.ceil(box_width))), dp_cat='SCIENCE')
    
        res.meta['TQUERY'] = (NOW, 'Timestamp of query execution')
        res.meta['RA'] = (ra, 'Query center, RA')
        res.meta['DEC'] = (dec, 'Query center, Dec')
        res.meta['R'] = (query_size, 'Query radius, arcmin')
        res.meta['N'] = len(res)
        res.meta['NAME'] = tab.meta['NAME']
    
        if len(res) > 0:
            res['field_root'] = res.meta['NAME'].lower()
            
        res.write(muse_file, overwrite=True, format='ascii.ecsv')
    else:
        res = Table.read(muse_file, format='ascii.ecsv')
        
    if make_figure & (len(res) > 0):
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        #ax.set_xlim(xr)
        #ax.set_ylim(yr)

        # Grow by factor of 2
        #expand = 2
        
        expand = np.maximum(min_size*2/np.minimum(dx, dy), 1)
        
        ax.set_xlim(ra+dx/60/2*expand/cosd, ra-dx/60/2*expand/cosd)
        ax.set_ylim(dec-dy/60/2*expand, dec+dy/60/2*expand)
        
        ax.scatter(ra, dec, marker='+', color='k')
        
        # HST patch
        p0 = np.array([[xr[0], yr[0]],
                       [xr[0], yr[1]],
                       [xr[1], yr[1]],
                       [xr[1], yr[0]]])

        p_hst = None
        for fph in tab['footprint']:
            for p in query.parse_polygons(fph):
                p_j = Polygon(p).buffer(0.001)
                if p_hst is None:
                    p_hst = p_j
                else:
                    p_hst = p_hst.union(p_j)

        ax.add_patch(patch_from_polygon(p_hst, ec='k', fc='None', 
                              alpha=0.8, label='HST'))
        
        band_labels = []
        
        is_public = res['Release Date'] < NOW
        
        for i in range(len(res)):
            fpstr = 'CIRCLE {0} {1} 0.008333'.format(res['RA'][i], res['DEC'][i])
            fps = query.parse_polygons(fpstr)
            
            if 'NOAO' in res['INS MODE'][i]:
                color = '#1f77b4'
            else:
                color = '#d62728'
                            
            if is_public[i]:
                linestyle='-'
            else:
                linestyle='--'
                                    
            for j, fp in enumerate(fps):
                fp_j = Polygon(fp).buffer(0.1/3600.)
                if res['INS MODE'][i] in band_labels:
                    label = None
                else:
                    band_labels.append(res['INS MODE'][i])
                    label = '{0}'.format(res['INS MODE'][i])
                    
                ax.add_patch(patch_from_polygon(fp_j, ec=color, fc='None', 
                                      alpha=0.8, label=label, 
                                      linestyle=linestyle))
                
                if is_public[i]:
                    ax.add_patch(patch_from_polygon(fp_j, ec=color, fc=color, 
                                          alpha=0.1+0.1*is_public[i]))
                                   
        ax.grid()
        ax.set_title('{0} MUSE'.format(res.meta['NAME']))
        
        #ax.set_xlim(ax.get_xlim()[::-1])
        
        ax.set_aspect(1/cosd)
        ax.legend(ncol=1, fontsize=6, loc='upper right')
        fig.set_size_inches(xsize, xsize*np.clip(dy/dx, 0.2, 5))
        
        if nlabel > 0:
            draw_axis_labels(ax=ax, nlabel=nlabel)
        
        ax.text(0.03, 0.03, NOW, fontsize=5, transform=ax.transAxes, ha='left', va='bottom')

        fig.tight_layout(pad=0.2)
        fig.savefig('{0}_muse.png'.format(meta['NAME']))
        
    else:
        fig = None
        
    return res, fig


def alma_query(tab, make_figure=True, xsize=5, nlabel=3, min_size=4, cmap='jet_r', rerun_query=True, query_kwargs={'public':False, 'science':False, 'get_html_version':True}):
    """
    Query ALMA archive around the HST data
    """
    import time
    import urllib
    import numpy as np
    import matplotlib.pyplot as plt
    
    from astropy.table import Table
    from astropy.time import Time
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    
    from shapely.geometry import Polygon
    #from descartes import PolygonPatch
    
    from astroquery.alma import Alma
    from astroquery.alma import utils as alma_utils
    
    NOW = Time.now().iso
    
    meta = tab.meta
    xr = (meta['XMIN'], meta['XMAX'])
    yr = (meta['YMIN'], meta['YMAX'])
    ra, dec = meta['BOXRA'], meta['BOXDEC']
    
    cosd = np.cos(dec/180*np.pi)
    dx = (xr[1]-xr[0])*cosd*60
    dy = (yr[1]-yr[0])*60
    
    box_width = np.maximum(dx, dy)
    query_size = np.maximum(min_size, box_width/2)
    
    coo = SkyCoord(ra, dec, unit='deg')
    
    alma_file = '{0}_alma.ecsv'.format(tab.meta['NAME'])
    if (not os.path.exists(alma_file)) | rerun_query:
        res = Alma.query_region(coo, query_size*u.arcmin, **query_kwargs)
    
        res.meta['TQUERY'] = (NOW, 'Timestamp of query execution')
        res.meta['RA'] = (ra, 'Query center, RA')
        res.meta['DEC'] = (dec, 'Query center, Dec')
        res.meta['R'] = (query_size, 'Query radius, arcmin')
        res.meta['N'] = len(res)
        res.meta['NAME'] = tab.meta['NAME']
    
        if len(res) > 0:
            res['field_root'] = res.meta['NAME'].lower()
            
        res.write(alma_file, overwrite=True, format='ascii.ecsv')
    else:
        res = Table.read(alma_file, format='ascii.ecsv')
        
    if make_figure & (len(res) > 0):
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        #ax.set_xlim(xr)
        #ax.set_ylim(yr)

        # Grow by factor of 2
        #expand = 2
        
        expand = np.maximum(min_size*2/np.minimum(dx, dy), 1)
        
        ax.set_xlim(ra+dx/60/2*expand/cosd, ra-dx/60/2*expand/cosd)
        ax.set_ylim(dec-dy/60/2*expand, dec+dy/60/2*expand)
        
        ax.scatter(ra, dec, marker='+', color='k')
        
        # HST patch
        p0 = np.array([[xr[0], yr[0]],
                       [xr[0], yr[1]],
                       [xr[1], yr[1]],
                       [xr[1], yr[0]]])

        p_hst = None
        for fph in tab['footprint']:
            for p in query.parse_polygons(fph):
                p_j = Polygon(p).buffer(0.001)
                if p_hst is None:
                    p_hst = p_j
                else:
                    p_hst = p_hst.union(p_j)

        ax.add_patch(patch_from_polygon(p_hst, ec='k', fc='None', 
                              alpha=0.8, label='HST'))
        
        band_labels = []
        so = np.argsort(res['Band'])
        
        is_public = res['Release date'] < NOW
        
        for i in so:
            fpstr = res['Footprint'][i]
            fps = query.parse_polygons(fpstr)
            is_mosaic = res['Mosaic'][i] in ['mosaic']
            try:
                color = plt.cm.get_cmap(cmap)(int(res['Band'][i])/10)
            except:
                color = 'r'
            
            if is_public[i]:
                linestyle='-'
            else:
                linestyle='--'
                                    
            for j, fp in enumerate(fps):
                fp_j = Polygon(fp).buffer(0.1/3600.)
                if res['Band'][i] in band_labels:
                    label = None
                else:
                    band_labels.append(res['Band'][i])
                    label = 'Band {0}'.format(res['Band'][i])
                    
                ax.add_patch(patch_from_polygon(fp_j, ec=color, fc='None', 
                                      alpha=0.8, label=label, 
                                      linestyle=linestyle))
                
                if is_mosaic & is_public[i]:
                    ax.add_patch(patch_from_polygon(fp_j, ec=color, fc=color, 
                                          alpha=0.1+0.1*is_public[i]))
                                   
        ax.grid()
        ax.set_title('{0} ALMA'.format(res.meta['NAME']))
        
        #ax.set_xlim(ax.get_xlim()[::-1])
        
        ax.set_aspect(1/cosd)
        ax.legend(ncol=1, fontsize=6, loc='upper right')
        fig.set_size_inches(xsize, xsize*np.clip(dy/dx, 0.2, 5))
        
        if nlabel > 0:
            draw_axis_labels(ax=ax, nlabel=nlabel)
        
        ax.text(0.03, 0.03, NOW, fontsize=5, transform=ax.transAxes, ha='left', va='bottom')

        fig.tight_layout(pad=0.2)
        fig.savefig('{0}_alma.png'.format(meta['NAME']))
    else:
        fig = None
        
    return res, fig
    
#### URLs for IPAC / Spitzer queries

# SHA web query
SHA_URL = "https://sha.ipac.caltech.edu/applications/Spitzer/SHA/#id=SearchByPosition&RequestClass=ServerRequest&DoSearch=true&SearchByPosition.field.radius={size}&SearchByPosition.field.matchByAOR=false&UserTargetWorldPt={ra};{dec};EQ_J2000&SimpleTargetPanel.field.resolvedBy=nedthensimbad&MoreOptions.field.prodtype=aor,pbcd,bcd,supermosaic,inventory&InstrumentPanel.field.irac=_all_&InstrumentPanel.field.mips=_all_&InstrumentPanel.field.irs=_none_&InstrumentPanel.field.panel=instrument&InventorySearch.radius={size}&shortDesc=Position&isBookmarkAble=true&isDrillDownRoot=true&isSearchResult=true"

# API table output
IPAC_URL = "https://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/DataService?RA={ra}&DEC={dec}&SIZE={size}&VERB=3&DATASET=ivo%3A%2F%2Firsa.ipac%2Fspitzer.level{level}"

IPAC_AOR_URL = "https://irsa.ipac.caltech.edu/applications/Spitzer/SHA/servlet/DataService?REQKEY={aor}&VERB=3&DATASET=ivo%3A%2F%2Firsa.ipac%2Fspitzer.level{level}"
        
def spitzer_query(tab, level=1, make_figure=True, cmap='Spectral', xsize=6, nlabel=3, min_size=4, maxwavelength=20):
    """
    Query the Spitzer archive around an HST pointing table
    """
    import time
    import urllib
    import numpy as np
    import matplotlib.pyplot as plt
    
    from astropy.table import Table
    from astropy.time import Time
    
    from shapely.geometry import Polygon
    #from descartes import PolygonPatch
           
    if False:
        tab = utils.read_catalog('j021744m0346_footprint.fits')
        tab = utils.read_catalog('j224324m0936_footprint.fits')
        root='j012016p2134'
        tab = utils.read_catalog('{0}_footprint.fits'.format(root))
    
    if isinstance(tab, str):
        tab = utils.read_catalog('{0}_footprint.fits'.format(tab))
            
    meta = tab.meta
    xr = (meta['XMIN'], meta['XMAX'])
    yr = (meta['YMIN'], meta['YMAX'])
    ra, dec = meta['BOXRA'], meta['BOXDEC']
    
    cosd = np.cos(dec/180*np.pi)
    dx = (xr[1]-xr[0])*cosd*60
    dy = (yr[1]-yr[0])*60
    
    box_width = np.maximum(dx, dy)
    query_size = np.maximum(min_size, box_width/2)/60.
    
    if level == -1:
        # First query Level 2 and then get by AOR
        query_url = IPAC_URL.format(ra=ra, dec=dec, size=query_size, level=2)
        
        f = urllib.request.urlopen(query_url)
        data = f.read().decode('utf-8')
        f.close()
        
        lev2 = Table.read(data, format='ascii.ipac')
        skip = np.array([l.startswith('IRS') for l in lev2['modedisplayname']])
        
        # only up to 24um
        skip |= lev2['minwavelength'] > maxwavelength
        
        aors = np.unique(lev2['reqkey'][~skip])
        print('      level=-1, query by N={0} AORs'.format(len(aors)))
        
        aor_string = ','.join(['{0}'.format(aor) for aor in aors])
        
        query_url = IPAC_AOR_URL.format(aor=aor_string, level=1)
        
        f = urllib.request.urlopen(query_url)
        data = f.read().decode('utf-8')
        f.close()
        
        level = 1
        
    else:
        query_url = IPAC_URL.format(ra=ra, dec=dec, size=query_size, level=level)
    
        f = urllib.request.urlopen(query_url)
        data = f.read().decode('utf-8')
        f.close()
    
    mode_data = {}
    modes = ['IRAC 3.6um', 'IRAC 4.5um', 'IRAC 5.8um', 'IRAC 8.0um', 
             'MIPS 24um', 'MIPS 70um', 'MIPS 160um']
    
    try:
        ipac = Table.read(data, format='ascii.ipac')
        
        if level == 2:
            #ipac['exptime'] = ipac['endtime'] - ipac['begintime']
            dt = Time(ipac['endtime']) - Time(ipac['begintime']) 
            ipac['exposuretime'] = dt.sec
             
        for mode in modes:
            m = ipac['wavelength'] == mode
            mode_i = {}
            mode_i['exptime'] = ipac['exposuretime'][m].sum()
            mode_i['N'] = m.sum()
            mstr = mode.replace(' ','').replace('.','')[:-2]
            mode_i['short'] = mstr
            mode_data[mode] = mode_i

            ipac.meta['T{0}'.format(mstr)] = int(mode_i['exptime']), 'Exposure time in {0}, seconds'.format(mode)
            ipac.meta['N{0}'.format(mstr)] = mode_i['N'], 'Number of BCD files for {0}'.format(mode)
    
    except:
        ipac = Table()
        ipac['x'] = [0]
        
        for mode in modes:
            mode_i = {}
            mode_i['exptime'] = 0.
            mode_i['N'] = 0
            mstr = mode.replace(' ','').replace('.','')[:-2]
            mode_i['short'] = mstr
            mode_data[mode] = mode_i

            ipac.meta['T_{0}'.format(mstr)] = int(mode_i['exptime']), 'Exposure time in {0}, seconds'.format(mode)
            ipac.meta['N_{0}'.format(mstr)] = mode_i['N'], 'Number of BCD files for {0}'.format(mode)
                    
    ipac.meta['URL'] = query_url
    ipac.meta['TQUERY'] = (time.ctime(), 'Timestamp of query execution')
    ipac.meta['LEVEL'] = (level, 'IPAC calibration level')
    ipac.meta['RA'] = (ra, 'Query center, RA')
    ipac.meta['DEC'] = (dec, 'Query center, Dec')
    ipac.meta['R'] = (query_size, 'Query radius, deg')
    
    ipac.write('{0}_ipac.fits'.format(meta['NAME']), overwrite=True)
    
    if 'x' in ipac.colnames:
        return ipac, None
        
    if make_figure:
        colors = list(plt.cm.get_cmap(cmap)(np.linspace(0, 0.35, 4)))
        colors += list(plt.cm.get_cmap(cmap)(np.linspace(0.6, 1, 3))[::-1])

        #modes = np.unique(ipac['wavelength'])
        fp = np.array([ipac['ra1'], ipac['dec1'], ipac['ra2'], ipac['dec2'], 
                       ipac['ra3'], ipac['dec3'], ipac['ra4'], ipac['dec4']])

        fig = plt.figure()
        ax = fig.add_subplot(111)
        #ax.set_xlim(xr)
        #ax.set_ylim(yr)

        # Grow by factor of 2
        #expand = 2
        
        expand = np.maximum(min_size*2/np.minimum(dx, dy), 1)
        
        ax.set_xlim(ra+dx/60/2*expand/cosd, ra-dx/60/2*expand/cosd)
        ax.set_ylim(dec-dy/60/2*expand, dec+dy/60/2*expand)
        
        ax.scatter(ra, dec, marker='+', color='k')
        
        # HST patch
        p0 = np.array([[xr[0], yr[0]],
                       [xr[0], yr[1]],
                       [xr[1], yr[1]],
                       [xr[1], yr[0]]])

        p_hst = None
        for fph in tab['footprint']:
            for p in query.parse_polygons(fph):
                p_j = Polygon(p).buffer(0.001)
                if p_hst is None:
                    p_hst = p_j
                else:
                    p_hst = p_hst.union(p_j)
        
        ipac['with_hst'] = False
        
        for i, mode in enumerate(modes):
            m = ipac['wavelength'] == mode
            mode_i = mode_data[mode]
            if m.sum() == 0:
                mode_i['fp'] = None
                continue
                
            fps = fp[:,m]
            fp_i = Polygon(fps[:,0].reshape((-1,2))).buffer(0.1/3600.)
            fp_i_hst = Polygon(fps[:,0].reshape((-1,2))).buffer(0.1/3600.)
            with_hst = np.zeros(m.sum(), dtype=bool)
            for j, f in enumerate(fps.T):
                fp_j = Polygon(f.reshape((-1,2))).buffer(0.1/3600.)
                with_hst[j] = p_hst.intersection(fp_j).area > 0
                if with_hst[j]:
                    fp_i_hst = fp_i_hst.union(fp_j)
                    ax.add_patch(patch_from_polygon(fp_j, ec=colors[i],
                                                    fc='None', 
                                                    alpha=0.05))
                    
                fp_i = fp_i.union(fp_j)
             
            mode_i['fp'] = fp_i
            ipac['with_hst'][m] |= with_hst
            
            label = '{0} - {1:3.1f}/{2:3.1f} hr'.format(mode_i['short'], mode_i['exptime']/3600., ipac['exposuretime'][m][with_hst].sum()/3600.)
            
            if mode_i['short'] in ['IRAC36', 'IRAC45']:
                ax.add_patch(patch_from_polygon(fp_i_hst, ec=colors[i],
                                          fc=colors[i], 
                                          alpha=0.2, label=label))
                
                ax.add_patch(patch_from_polygon(fp_i,
                                      ec=colors[i], fc=colors[i], 
                                      alpha=0.05, label=label))
                
            else:
                ax.add_patch(patch_from_polygon(fp_i,
                                      ec=colors[i], fc='None', 
                                      alpha=0.8, label=label, 
                                      linestyle='--'))
                
        ax.add_patch(patch_from_polygon(p_hst, ec='k', fc='None', 
                              alpha=0.8, label='HST'))
               
        ax.grid()
        ax.set_title('{0} - Level{1}'.format(meta['NAME'], level))
        
        #ax.set_xlim(ax.get_xlim()[::-1])
        
        ax.set_aspect(1/cosd)
        ax.legend(ncol=2, fontsize=6, loc='upper right')
        fig.set_size_inches(xsize, xsize*np.clip(dy/dx, 0.2, 5))
        
        if nlabel > 0:
            draw_axis_labels(ax=ax, nlabel=nlabel)
        
        ax.text(0.03, 0.03, time.ctime(), fontsize=5, transform=ax.transAxes, ha='left', va='bottom')

        fig.tight_layout(pad=0.2)
        fig.savefig('{0}_ipac.png'.format(meta['NAME']))
        ipac.write('{0}_ipac.fits'.format(meta['NAME']), overwrite=True)
    else:
        fig = None
        
    return ipac, fig


def make_all():
    """
    Make all IRAC queries
    """
    import glob
    import matplotlib.pyplot as plt
    from mastquery import overlaps
    
    files = glob.glob('*footprint.fits')
    files.sort()
    plt.ioff()
    failed = []
    for file in files:
        if os.path.exists(file.replace('_footprint', '_ipac')):
            print('skip {0}'.format(file))
            continue

        else:
            print(file)
        
        tab = utils.read_catalog(file)
        if 'NASSOC' in tab.meta:
            test = tab.meta['NASSOC'] < 100
        else:
            test = True
        
        if 'XMIN' not in tab.meta:
            print('Split {0}'.format(file))
            _res = overlaps.split_associations(tab)
                
        if test:
            try:
                ipac, fig = spitzer_query(tab, min_size=10, level=2)
                if len(ipac) < 50:
                    ipac, fig = spitzer_query(tab, min_size=10, level=1)
                    
                plt.close('all')
            except:
                failed.append(file)
    
    plt.ion()
    
    if False:
        import os
        import numpy as np
        import matplotlib.pyplot as plt                                                                                                         
        from mastquery.overlaps import spitzer_query
        from grizli import utils
        
        fields = np.loadtxt('fields.list', dtype=str)
        
        #for field in fields:
        
        files = glob.glob('*footprint.fits')
        files.sort()
        for file in files:
            if os.path.exists(file.replace('footprint','ipac')):
                continue
            
            field = file.split("_footprint")[0]
            tab = utils.read_catalog('{0}_footprint.fits'.format(field))
            print(field, tab.meta['BOXRAD'])
            if tab.meta['BOXRAD'] > 11:
                continue
                
            ipac, fig = spitzer_query(tab, min_size=10, level=-1)
            plt.close('all')
    