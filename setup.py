#from distutils.core import setup
#from distutils.extension import Extension
from setuptools import setup
from setuptools.extension import Extension

import subprocess

import os

#update version
args = 'git describe --tags'
p = subprocess.Popen(args.split(), stdout=subprocess.PIPE)
long_version = p.communicate()[0].decode("utf-8").strip()
spl = long_version.split('-')

if len(spl) == 3:
    main_version = spl[0]
    commit_number = spl[1]
    version_hash = spl[2]
    version = f'{main_version}.dev{commit_number}'
else:
    version_hash = '---'
    version = long_version
    
# version = "0.1.8"
# version = "0.1.9" # Fix file extensions
# version = "0.2.0" # Minor changes
# version = "0.3.0" # New target formatting
# version = "0.4.0" # Try to use astroquery rather than MAST CAOM
# version = "0.4.1" # Bug in shapes
#version = "1.0" # Stable release, dustmaps dust
#version = "1.0.1" # Use new get_mw_dust in overlaps
#version = "1.1" # mast products workaround
# version = "1.2" # jwst
# version = "1.3" # polystr
# version = "1.4" # force_rate in fetch
# version = "1.5" # guide star query in jwst.query_guidestar_log
version = "1.5.1" # bug fix for pip versioning

# Set this to true to add install_requires to setup
# Turned off for incremental builds as it kills "reload(mastquery.query)" 
if 0:
    install_requires=[
         'astropy>=2.0.0',
         'astroquery>=0.3.0',
         'lxml>=3.8.0',
         'numpy>=1.10.2',
         'geos>=0.2.1',
         'shapely>=1.5.16',
         'matplotlib>=2.0.2',
         'descartes>=1.0.2']
else:
    install_requires = []    
    
#lines = open('grizli/version.py').readlines()
version_str =f"""# git describe --tags
__version__ = "{version}"
__long_version__ = "{long_version}"
__version_hash__ = "{version_hash}"
"""

fp = open('mastquery/version.py','w')
fp.write(version_str)
fp.close()
print('Git version: {0}'.format(version))

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "mastquery",
    version = version,
    author = "Gabriel Brammer",
    author_email = "gbrammer@gmail.com",
    description = "Python tools for querying the MAST Archive",
    license = "MIT",
    url = "https://github.com/gbrammer/mastquery",
    download_url = "https://github.com/gbrammer/mastquery/tarball/{0}".format(version),
    packages=['mastquery'],
    classifiers=[
        "Development Status :: 1 - Planning",
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
    install_requires=install_requires,
    package_data={'mastquery': []},
)
