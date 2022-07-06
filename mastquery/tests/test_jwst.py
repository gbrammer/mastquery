import numpy as np
from .. import jwst

def test_jwst_query():
    
    filters = jwst.CALIB_FILTERS
    filters += jwst.FULL_SUBARRAY
    filters += jwst.FINE_GUIDE
    
    filters += [{'paramName': 'expstart',
                 'values': [{'min': 59730.4, 'max': 59735.4}]}]
    
    res = jwst.query_all_jwst(filters=filters, fix=False, columns=None)
    res = jwst.query_all_jwst(filters=filters, fix=False, columns='*')
    
    assert(len(res) > 0)
    # assert(np.abs(len(res)-953) < 50)
    
    