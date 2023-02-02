import os
import glob
import matplotlib.pyplot as plt

from .. import query, overlaps

QUERY = None

def test_query():
    """
    Test basic query
    """
    global QUERY
    
    _res = query.run_query(proposal_id=[11359], filters=['G141'])
    assert(len(_res) == 1)
    
    olap = overlaps.find_overlaps(_res, proposal_id=[11359], 
                                  instruments=['WFC3/IR'])
    assert(len(olap) == 1)
    
    _ = overlaps.split_associations(olap[0])
    
    plt.close('all')
    
    # cleanup
    files = glob.glob(f"{olap[0].meta['NAME']}*")
    files += glob.glob('overlaps.yaml')
    
    for file in files:
        os.remove(file)
    
    QUERY = _res

def test_products():
    """
    Get products info
    """
    global QUERY
    
    if QUERY is not None:
        _prod = query.get_products_table(QUERY, extensions=['RAW'], 
                                         use_astroquery=True)
        
        assert(len(_prod) == 4)
    
    