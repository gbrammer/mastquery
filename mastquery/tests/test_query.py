
from .. import query

def test_query():
    """
    Test basic query
    """
    _res = query.run_query(proposal_id=[11359], filters=['G141'])

    assert(len(_res) == 1)
    return _res

def test_products():
    """
    Get products info
    """
    _res = test_query()
    _prod = query.get_products_table(_res, extensions=['RAW'], 
                                     use_astroquery=True)
    assert(len(_prod) == 4)
    
    