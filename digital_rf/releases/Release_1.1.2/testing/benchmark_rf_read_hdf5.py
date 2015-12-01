"""benchmark_rf_read_hdf5.py is a script to benchmark reading the Hdf5 files produced by benchmark_rf_write_hdf5.py

$Id: benchmark_rf_read_hdf5.py 636 2014-09-26 18:34:12Z brideout $
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
    if read_size > FILE_SAMPLES:
        raise ValueError, 'read_size %i must not be greater than FILE_SAMPLES %i' % (read_size, FILE_SAMPLES)
    if sys.platform == 'darwin':
        os.system('purge')
    t2 = time.time()
    start_index, end_index = test_read_obj.get_bounds(channel_name)
    print('get_bounds returned %i - %i and took %f seconds' % (start_index, end_index, time.time() - t2))
    t2 = time.time()
    result = test_read_obj.get_continuous_blocks(start_index, end_index, channel_name)
    print('get_continuous_blocks returned <%s> and took %f seconds' % (str(result), time.time() - t2))
    
    next_sample = start_index
    t2 = time.time()
    count = 0
    for i in range(1000):
        arr = test_read_obj.read_vector_raw(next_sample, read_size, channel_name)
        next_sample += read_size
        count += 1
        if count % 100 == 0:
            print('%i out of 1000' % (count))
    seconds = time.time() - t2
    speedMB = (read_size*1000*4)/(1.0E6*seconds)
    print('Total read time %i seconds, speed %1.2f MB/s' % (int(seconds), speedMB))

t = time.time()
test_read_obj = digital_rf_hdf5.read_hdf5('/tmp/benchmark', load_all_metadata=True)
print('metadata analysis took %f seconds' % (time.time() - t))

print("\nTest 0 - read Hdf5 files with no compress, no checksum - channel name = junk0")
test_read('junk0', test_read_obj)

print("\nTest 0.1 - read Hdf5 files with no compress, no checksum, small read size - channel name = junk0")
test_read('junk0', test_read_obj, 1000)

print("\nTest call to reload to update metadata")
t = time.time()
test_read_obj.reload()
print('reload took %f seconds' % (time.time() - t))


print("\nTest 1 -read Hdf5 files with no compress, but with level 9 checksum - channel name = junk1")
test_read('junk1', test_read_obj)


print("\nTest 2 - read Hdf5 files with compress, and with level 9 checksum - channel name = junk2")
test_read('junk2', test_read_obj)
