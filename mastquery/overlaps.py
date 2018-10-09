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
    
def find_overlaps(tab, buffer_arcmin=1., filters=[], instruments=['WFC3/IR', 'WFC3/UVIS', 'ACS/WFC'], proposal_id=[], SKIP=False, base_query=query.DEFAULT_QUERY, close=True, use_parent=False, extensions=['FLT','C1M'], include_subarrays=False, min_area=0.2, show_parent=True):
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
        Extra query parameters.
        
    close : bool
        If true, close the figure objects.
    
    use_parent : bool
        Use parent table rather than performing a new query
        
    Returns
    -------
    tables : list
        
        List of grouped tables (`~astropy.table.Table`).

    """
    import copy
    import os
    
    import numpy as np
    import matplotlib.pyplot as plt

    from shapely.geometry import Polygon
    from descartes import PolygonPatch
        
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
        #print('Parse', i)
        has_match = False
        for j in range(len(match_poly)):
            isect = match_poly[j].intersection(polygons[i])
            #print(tab['target'][i], i, isect.area > 0)
            if isect.area*3600. > 0.5:
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
            #print(i)
            has_match = False
            for j in range(len(match_poly)):
                isect = match_poly[j].intersection(mpolygons[i])
                #print(tab['target'][i], i, isect.area > 0)
                if isect.area > 0:
                    #print(isect.area*3600)
                    match_poly[j] = match_poly[j].union(mpolygons[i])
                    match_ids[j].extend(mids[i])
                    has_match = True
                    continue

            if not has_match:
                match_poly.append(mpolygons[i])
                match_ids.append(mids[i])
        
        print('Iter #{0}, N_Patch = {1}'.format(iter+1, len(match_poly)))
        if len(mpolygons) == len(match_poly):
            break
    
    np.save('overlaps.py', [match_poly, match_ids])
            
    # Save figures and tables for the unique positions
    BLUE = '#6699cc'
    
    tables = []
    
    for i in range(len(match_poly)):
        #i+=1
        p = match_poly[i]
                
        #######
        # Query around central coordinate
        idx = np.array(match_ids[i])
        ra, dec = np.mean(tab['ra'][idx]), np.mean(tab['dec'][idx])
        
        # Get poly size
        xy = p.convex_hull.boundary.xy
        xradius = np.abs(xy[0]-ra).max()*np.cos(dec/180*np.pi)*60
        yradius = np.abs(xy[1]-dec).max()*60

        box = [ra, dec, np.maximum(xradius*1.5, yradius*1.5)]
        
        # Build target name from RA/Dec
        jname = utils.radec_to_targname(box[0], box[1], scl=1000)
        print('\n\n', i, jname, box[0], box[1])

        if (os.path.exists('{0}_footprint.pdf'.format(jname))) & SKIP:
            continue
                            
        if use_parent:
            xtab = tab
        else:
            try:
                xtab = query.run_query(box=box, proposal_id=proposal_id, instruments=instruments, filters=filters, base_query=base_query)
            except:
                pass
                
            if not include_subarrays:
                query.set_area_column(xtab)
                subarray = (xtab['instrument_name'] == 'WFC3/IR') & (xtab['area'] < 3.8)
                xtab = xtab[~subarray]
            
        ebv = utils.get_irsa_dust(ra, dec, type='SandF')
        xtab.meta['NAME'] = jname
        xtab.meta['RA'] = ra
        xtab.meta['DEC'] = dec
        xtab.meta['MW_EBV'] = ebv
              
        # Only include ancillary data that directly overlaps with the primary
        # polygon
        pointing_overlaps = np.zeros(len(xtab), dtype=bool)
        for j in range(len(xtab)):
            poly = query.parse_polygons(xtab['footprint'][j])#[0]
            pshape = [Polygon(pi) for pi in poly]
            for k in range(1,len(poly)):
                pshape[0] = pshape[0].union(pshape[k])
            
            try:
                isect = p.intersection(pshape[0])
                pointing_overlaps[j] = isect.area > min_area*pshape[0].area
            except:
                pointing_overlaps[j] = False
                
        xtab = xtab[pointing_overlaps]
        
        # Unique targets
        filter_target = np.array(['{0} {1}'.format(xtab['instrument_name'][i], xtab['filter'][i]) for i in range(pointing_overlaps.sum())])
                        
        ########### 
        # Make the figure
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        # Show the parent table
        if show_parent:
            colors = query.show_footprints(tab[idx], ax=ax)
        
        ax.scatter(box[0], box[1], marker='+', color='k')
        
        colors = query.show_footprints(xtab, ax=ax)
        
        patch1 = PolygonPatch(p, fc=BLUE, ec=BLUE, alpha=0.1, zorder=2)
        
        ax.plot(xy[0], xy[1])
        ax.add_patch(patch1)
        
        ax.grid()
        
        # Resize for square dimensions
        xr, yr = ax.get_xlim(), ax.get_ylim()
        dx = (xr[1]-xr[0])*np.cos(yr[0]/180*np.pi)*60
        dy = (yr[1]-yr[0])*60
        ax.set_title(jname)
        ax.set_xlim(ax.get_xlim()[::-1])
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
    
    try:
        from grizli import utils
        HAS_GRIZLI = True
    except:
        HAS_GRIZLI = False
        
    if tabs is None:
        tabs = [Table.read(file) for file in glob.glob('*footprint.fits')]
    
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
    
    