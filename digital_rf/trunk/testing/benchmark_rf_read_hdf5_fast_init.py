"""benchmark_rf_read_hdf5_fast_init.py is a script to benchmark reading the Hdf5 files produced by benchmark_rf_write_hdf5.py

Differs from benchmark_rf_read_hdf5.py is that only sparse metadata used. 

$Id: benchmark_rf_read_hdf5_fast_init.py 669 2014-10-15 18:44:44Z brideout $
"""
# standard python imports
import os, os.path, sys
import time

# third party imports
import numpy

# Millstone imports
import digital_rf_hdf5

# constants
WRITE_BLOCK_SIZE = 1024
FILE_SAMPLES = 1000000
N_WRITES = int(1e9/WRITE_BLOCK_SIZE)


def test_read(channel_name,test_read_obj, read_size=FILE_SAMPLES):
    """test_read measures the speed of reading all the data back into numpy arrays for
    a given channel_name and digital_rf_hdf5.read_hdf5 object
    
    read_size must not be greater than FILES_SAMPLES (or you run out of data)
    """
    if sys.platform == 'darwin':
        os.system('purge')
    t2 = time.time()
    start_index, end_index = test_read_obj.get_bounds(channel_name)
    print('get_bounds returned %i - %i and took %f seconds' % (start_index, end_index, time.time() - t2))
    
    next_sample = start_index
    t2 = time.time()
    count = 0
    for i in range(100000):
        arr = test_read_obj.read_vector_raw(next_sample, read_size, channel_name)
        next_sample += read_size
        count += 1
        if count % 10000 == 0:
            print('%i out of 100000' % (count))
    seconds = time.time() - t2
    speedMB = (read_size*100000*4)/(1.0E6*seconds)
    print('Total read time %i seconds, speed %1.2f MB/s' % (int(seconds), speedMB))

t = time.time()
test_read_obj = digital_rf_hdf5.read_hdf5('/tmp/benchmark', load_all_metadata=False)
print('metadata analysis took %f seconds' % (time.time() - t))

print("\nTest 0.1 - read Hdf5 files with no compress, no checksum, small read size - channel name = junk0")
test_read('junk0', test_read_obj, 1000)

print("\nTest call to reload to update metadata")
t = time.time()
test_read_obj.reload()
print('reload took %f seconds' % (time.time() - t))


print("\nTest 1 -read Hdf5 files with no compress, but with level 9 checksum - channel name = junk1")
test_read('junk1', test_read_obj, 1000)


print("\nTest 2 - read Hdf5 files with compress, and with level 9 checksum - channel name = junk2")
test_read('junk2', test_read_obj, 1000)
