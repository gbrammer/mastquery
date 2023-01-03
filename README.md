# MASTQUERY
User-friendly tools for using the MAST Mashup API (https://mast.stsci.edu/api/v0/index.html)

## Installation:

    # From PIP
    pip install mastquery
    
    # *OR* from latest version of the respository
    git clone https://github.com/gbrammer/mastquery.git
    cd mastquery
    pip install . 
    
## Demo:

See also [demo.ipynb](https://github.com/gbrammer/mastquery/blob/master/examples/demo.ipynb).

```python
>>> from mastquery import query, fetch

### Query associations
>>> tab = query.run_query(box=None, proposal_id=[11359],
                         instruments=['WFC3/IR'], 
                         filters=['G141'],
                         base_query=query.DEFAULT_QUERY)

>>> print(tab['obs_id', 'filter', 'exptime', 'proposal_id'])
  obs_id   filter  exptime  proposal_id
========   ======  =======  ===========
ib6o23010   G141    4211.7        11359

### Data products
>>> prod = query.get_products_table(tab, extensions=['RAW'])

>>> print(prod['observation_id', 'filter', 'productFilename'])
observation_id filter  productFilename  
============== ====== ==================
     ib6o23rsq   G141 ib6o23rsq_raw.fits
     ib6o23ruq   G141 ib6o23ruq_raw.fits
     ib6o23ryq   G141 ib6o23ryq_raw.fits
     ib6o23s0q   G141 ib6o23s0q_raw.fits
     
### Fetch products
>>> s3_lines = fetch.make_curl_script(prod, script_name=None, s3_sync=True)
>>> print(s3_lines[0])
aws s3 sync --request-payer requester --exclude="*.*" --include="*raw.fits" s3://stpubdata/hst/public/ib6o/ib6o23rsq/ .//
``` 
