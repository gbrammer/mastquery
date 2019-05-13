"""
Scripts to find overlapping HST data
"""

from . import query, utils
  
def test():
    
    import copy
    
    import numpy as np
    import matplotlib.pyplot as plt

    from shapely.geometry import Polygon
    from descartes import PolygonPatch
    
    from hsaquery import query
    from hsaquery.query import parse_polygons

    
    # Example: high-z cluster pointings
    tab = query.run_query(box=None, proposid=[14594], instruments=['WFC3-IR', 'ACS-WFC'], extensions=['FLT'], filters=['G102','G141'], extra=[])
    
    # Ebeling
    tab = query.run_query(box=None, proposid=[15132,14098,13671,12884,12166,11103,10875,], instruments=['WFC3-IR'], extensions=['FLT'], filters=[], extra=[])
    # Relics
    tab = query.run_query(box=None, proposid=[14096], instruments=['WFC3-IR'], extensions=['FLT'], filters=[], extra=[])
    tab = query.run_query(box=None, proposid=[11591], instruments=['WFC3-iR'], extensions=['FLT'], filters=[], extra=[])
    tab = query.run_query(box=None, proposid=[13666,14148,14496], instruments=['WFC3-IR'], extensions=['FLT'], filters=[], extra=[])
    tab = tab[tab['target'] != 'ANY']
    
    # MACS 0454
    box = [73.5462181, -3.0147200, 3]
    tab = query.run_query(box=box, proposid=[], instruments=['WFC3-IR', 'ACS-WFC'], extensions=['FLT'], filters=['F110W'], extra=[])
    
def find_overlaps(tab, buffer_arcmin=1., filters=[], instruments=['WFC3/IR', 'WFC3/UVIS', 'ACS/WFC'], proposal_id=[], SKIP=False, base_query=query.DEFAULT_QUERY, extra={}, close=True, use_parent=False, extensions=['FLT','C1M'], include_subarrays=False, min_area=0.2, show_parent=True, show_parent_box=True, targstr='j{rah}{ram}{ras}{sign}{ded}{dem}', prefix='', suffix='', jstr='{prefix}{jname}{suffix}', fractional_overlap=0, patch_alpha=0.1, parent_alpha=0.1, tile_alpha=0.1, verbose=2, keep_single_name=True):
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
    import copy
    import os
    
    import numpy as np
    import matplotlib.pyplot as plt

    from astropy.table import Table

    from shapely.geometry import Polygon
    from descartes import PolygonPatch
        
    import time
    # Get shapely polygons for each exposures
    polygons = []
    
    poly_buffer = buffer_arcmin/60 # ~1 arcmin, but doesn't account for cos(dec)
    #poly_buffer = 0.5/60 # ~1 arcmin, but doesn't account for cos(dec)
    
    for i in range(len(tab)):
        poly = query.parse_polygons(tab['footprint'][i])#[0]
        pshape = [Polygon(p) for p in poly]
        for i in range(1,len(poly)):
            pshape[0] = pshape[0].union(pshape[i])
        
        polygons.append(pshape[0].buffer(poly_buffer))
    
    match_poly = [polygons[0]]
    match_ids = [[0]]
    
    # Loop through polygons and combine those that overlap
    for i in range(1,len(tab)):
        if verbose > 1:
            print(utils.NO_NEWLINE+'Parse {0:4d}'.format(i))
        
        has_match = False
        for j in range(len(match_poly)):
            isect = match_poly[j].intersection(polygons[i])
            #print(tab['target'][i], i, isect.area > 0)
            if fractional_overlap > 0:
                test_area = fractional_overlap*polygons[i].area
            else:
                test_area = 0.5/3600.
                
            if isect.area > test_area:
                #print(isect.area*3600)
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
    
        for i in range(1,len(mpolygons)):
            if verbose > 1:
                print(utils.NO_NEWLINE+'Parse, iter {0}, {1:4d}'.format(iter+1, i))
                
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
    
    np.save('overlaps.npy', [match_poly, match_ids])
            
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
        if (NPOLYS == 1) & keep_single_name:
            try:
                ra, dec = tab.meta['RA'], tab.meta['DEC']
            except:
                ra, dec = np.mean(tab['ra'][idx]), np.mean(tab['dec'][idx])                
        else:
            ra, dec = np.mean(tab['ra'][idx]), np.mean(tab['dec'][idx])
        
        # Get poly size
        xy = np.array(p.convex_hull.boundary.xy)
        xradius = np.abs(xy[0]-ra).max()*np.cos(dec/180*np.pi)*60
        yradius = np.abs(xy[1]-dec).max()*60

        box = [ra, dec, np.maximum(xradius*1.5, yradius*1.5)]
        
        # Build target name from RA/Dec
        jname = utils.radec_to_targname(box[0], box[1], round_arcsec=(4, 60), targstr=targstr)
        jname = jstr.format(prefix=prefix, jname=jname, suffix=suffix) #prefix+jname+suffix
        

        if (os.path.exists('{0}_footprint.pdf'.format(jname))) & SKIP:
            print('\n********** SKIP *********\n', i, jname, box[0], box[1])
            continue
        else:
            print('\n\n', i, jname, box[0], box[1])
            
        if use_parent:
            xtab = tab
        else:
            try:
                xtab = query.run_query(box=box, proposal_id=proposal_id, instruments=instruments, filters=filters, base_query=base_query, **extra)
            except:
                print('Failed!')
                pass
                
            if not include_subarrays:
                query.set_area_column(xtab)
                try:
                    subarray = (xtab['instrument_name'] == 'WFC3/IR') & (xtab['area'] < 3.8)
                except:
                    # Bug?
                    xtab['instrument_name'].fill_value = 0
                    subarray = (xtab['instrument_name'] == 'WFC3/IR') & (xtab['area'] < 3.8)
                    
                xtab = xtab[~subarray]
            
        ebv = utils.get_irsa_dust(ra, dec, type='SandF')
        xtab.meta['NAME'] = jname
        xtab.meta['RA'] = ra
        xtab.meta['DEC'] = dec
        xtab.meta['MW_EBV'] = ebv
        xtab.meta['FOLAP'] = (fractional_overlap, 'Fractional overlap parameter')
        xtab.meta['MIN_AREA'] = (min_area, 'Minimum overlap fraction')
        xtab.meta['BUFFER'] = (buffer_arcmin, 'Buffer for overlaps in arcmin')
        
        # Only include ancillary data that directly overlaps with the primary
        # polygon
        pointing_overlaps = np.zeros(len(xtab), dtype=bool)
        for j in range(len(xtab)):
            poly = query.parse_polygons(xtab['footprint'][j])#[0]
            pshape = [Polygon(pi) for pi in poly]
            for k in range(1,len(poly)):
                pshape[0] = pshape[0].union(pshape[k])
            
            try:
                isect = p.intersection(pshape[0].buffer(0.0001))
                pointing_overlaps[j] = isect.area > min_area*pshape[0].area
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
            colors = query.show_footprints(tab[idx], ax=ax, alpha=parent_alpha)
        
        ax.scatter(box[0], box[1], marker='+', color='k', zorder=1e4, alpha=0.8)
        
        if ('boxra' in tab.meta) & show_parent_box:
            ax.scatter(tab.meta['boxra'][0], tab.meta['boxdec'][0], marker='*', color='r', zorder=1e4, alpha=0.8)
            
        colors = query.show_footprints(xtab, ax=ax, alpha=tile_alpha)
        
        patch1 = PolygonPatch(p, fc=BLUE, ec=BLUE, alpha=patch_alpha, zorder=2)
        
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
        
        fig.set_size_inches(5,5*dy/dx)
        
        # Add summary
        fp = open('{0}_info.dat'.format(jname),'w')
        dyi = 0.02*dx/dy
        
        ax.text(0.05, 0.97, '{0:>13.5f} {1:>13.5f}  E(B-V)={2:.3f}'.format(ra, dec, ebv), ha='left', va='top', transform=ax.transAxes, fontsize=6)
        
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

        xtab.write('{0}_footprint.fits'.format(jname), format='fits', overwrite=True)
        np.save('{0}_footprint.npy'.format(jname), [p, box])
        
        tables.append(xtab)
    
    return tables
    
def summary_table(tabs=None, output='overlap_summary'):
    import glob
    from collections import OrderedDict
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
    import numpy as np
    from shapely.geometry import Polygon
    from mastquery import query# as mutils
    
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
            pshape = [Polygon(p) for p in poly]
            for j in range(1,len(poly)):
                pshape[0] = pshape[0].union(pshape[j])
            
            if i == 0:
                gpoly = pshape[0]
            else:
                gpoly = gpoly.union(pshape[0])
        
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
    
def compute_associations(tab, max_sep=0.5, max_pa=0.05, max_time=1e4, match_filter=True, match_instrument=True, match_program=True, hack_grism_pa=True, parse_for_grisms=True):
    """
    Associate visits by filter + position + PA + date
    """
    from . import query
    
    import numpy as np
    from shapely.geometry import Point
    
    cosd = np.cos(tab['dec']/180*np.pi)
    dx = (tab['ra'] - np.median(tab['ra']))*cosd*60    
    dy = (tab['dec'] - np.median(tab['dec']))*60
    
    ori = np.array([query.get_orientat(tab['footprint'][i]) for i in range(len(tab))])
     
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
    
    assoc = []
    for i in range(len(tab)):
        if assoc_idx[i] >= 0:
            continue
                
        assoc_idx[i] = len(assoc)
        
        assoc_i = {'pos':Point(dx[i], dy[i]).buffer(max_sep), 'ori':ori[i], 'filter':tab['filter'][i], 'indices':[i], 'proposal_id':tab['proposal_id'][i], 'idx':len(assoc), 't_min':tab['t_min'][i], 'instrument_name':tab['instrument_name'][i]}
        
        for j in range(i+1, len(tab)):
            
            if assoc_idx[j] >= 0:
                continue
                
            f_j = tab['filter'][j]
            pr_j = tab['proposal_id'][j]
            instr_j = tab['instrument_name'][j]
            dpa = assoc_i['ori'] - ori[j]
            dt = assoc_i['t_min'] - tab['t_min'][j]
            p_j = Point(dx[j], dy[j]).buffer(max_sep)
            
            # Has match
            test = (np.abs(dpa) < max_pa) & (p_j.intersects(assoc_i['pos']))
            
            test &= np.abs(dt) < max_time

            if match_instrument:
                test &= (instr_j == assoc_i['instrument_name'])
                
            if match_filter & (not parse_for_grisms):
                test &= (f_j == assoc_i['filter'])
            
            if match_program:
                test &= (pr_j == assoc_i['proposal_id'])
                
            if test:
                #print('Has match!', j)
                #break
                
                assoc_idx[j] = assoc_idx[i]
                assoc_i['pos'] = assoc_i['pos'].union(p_j)
                assoc_i['indices'].append(j)
        
        assoc.append(assoc_i)
    
    if match_filter & (parse_for_grisms):
        # Now split assoc that *don't* have grisms
        assoc_orig = assoc_idx*1
        new_i = 0
        assoc_idx = assoc_orig*0
        for ix in np.unique(assoc_orig):
            sel = assoc_orig == ix
            filts = list(np.unique(tab['filter'][sel]))
            is_grism = np.sum([f in filts for f in ['G800L','G102','G141']]) > 0
            
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
    
    if False:
        assoc = 48
        sel = tab['assoc_idx'] == assoc
        tabs = overlaps.find_overlaps(tab[sel], use_parent=True, buffer_arcmin=0.1, filters=['F814W'], proposal_id=[], instruments=['ACS/WFC'], close=False, suffix='-f606w-{0:02d}'.format(assoc))
        
    
    