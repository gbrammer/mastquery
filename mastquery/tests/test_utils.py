import os

from .. import utils

def test_api_download():
    """
    Test downloading with MAST API
    """
    import os
    import glob
    from astropy.table import Table
    
    _list = ['mast:HST/product/oetd2s020_sx1.fits']
    _file = os.path.basename(_list[0])
    
    _table = Table()
    _table['dataURL'] = _list
    
    _table_filename = Table()
    _table_filename['dataURL'] = _list
    _table_filename['filename'] = 'myfile.fits'
    
    _path='./'
    
    _file_list = utils.download_from_mast(_list, path=_path, min_size=0.05)
    assert(len(_file_list) == 1)
    assert(os.path.join(_path, _file) in _file_list)
    for _f in _file_list:
        os.remove(_f)
     
    _file_list = utils.download_from_mast(_table, path=_path, min_size=0.05)
    assert(len(_file_list) == 1)
    assert(os.path.join(_path, _file) in _file_list)
    for _f in _file_list:
        os.remove(_f)
    
    _file_list = utils.download_from_mast(_table_filename, path=_path,
                                          min_size=0.05)
    assert(len(_file_list) == 1)
    assert(os.path.join(_path, 'myfile.fits') in _file_list)
    for _f in _file_list:
        os.remove(_f)
    
    # Full path
    _subpath='./test_output/subdir'
    _file_list = utils.download_from_mast(_table_filename, path=_subpath,
                                          min_size=0.05)
    
    assert(os.path.exists(_subpath))
    assert(len(_file_list) == 1)
    assert(os.path.join(_subpath, 'myfile.fits') in _file_list)
    for _f in _file_list:
        os.remove(_f)
    
    os.removedirs(_subpath)