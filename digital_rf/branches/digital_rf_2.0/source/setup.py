"""setup file for the digital_rf_hdf5 python module
$Id: setup.py 786 2015-07-07 15:01:51Z brideout $
"""

from distutils.core import setup, Extension

setup(name="digital_rf_hdf5",
        version="1.1.3",
        description="Python tools to read and write digital rf data in Hdf5 format",
        author="Bill Rideout",
        author_email="brideout@haystack.mit.edu",
        url="http://www.haystack.mit.edu/~brideout/",
        py_modules=['digital_rf_hdf5'],
        ext_modules=[Extension("_py_rf_write_hdf5",
                              ["_py_rf_write_hdf5.c", "rf_write_hdf5.c"],
                              libraries=["hdf5"])
                    ])