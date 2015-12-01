"""benchmark_rf_write_hdf5.py is a script to to run the same test as in benchmark_rf_write_hdf5.c,
except in C.  All tests single subchannel.

$Id: benchmark_rf_write_hdf5.py 411 2014-05-23 20:37:12Z brideout $
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
FILES_PER_DIR = 1000
SAMPLE_RATE = 1E9

# start 2014-03-09 12:30:30 plus one sample
start_global_index = (1394368230 * SAMPLE_RATE) + 1


# data to write
data_int16 = numpy.zeros((WRITE_BLOCK_SIZE, 2), dtype='i2')
# make random
for i in range(WRITE_BLOCK_SIZE):
    j = i*2
    k = i*2 + 1
    data_int16[i][0]=(j%32768)*(j+8192)*(j%13);
    data_int16[i][1]=(k%32768)*(k+8192)*(k%13);
    
print('creating top level dir /tmp/benchmark')
os.system("rm -rf /tmp/benchmark ; mkdir /tmp/benchmark")

print("Test 0 - simple single write to multiple files, no compress, no checksum - channel 0")
os.system("rm -rf /tmp/benchmark/junk0 ; mkdir /tmp/benchmark/junk0")
print("Start writing")
channelObj = digital_rf_hdf5.write_hdf5_channel('/tmp/benchmark/junk0', 'i2', FILE_SAMPLES, FILES_PER_DIR, start_global_index, SAMPLE_RATE, 'Fake_uuid', 0, False)
t = time.time()
for i in range(N_WRITES):
    channelObj.rf_write(data_int16)
channelObj.close()
seconds = time.time() - t
speedMB = (N_WRITES*4*WRITE_BLOCK_SIZE)/(1.0E6*seconds)
print('Total time %i seconds, speed %1.2f MB/s' % (int(seconds), speedMB))


print("Test 1 - simple single write to multiple files, no compress, with checksum - channel 1")
os.system("rm -rf /tmp/benchmark/junk1 ; mkdir /tmp/benchmark/junk1")
print("Start writing")
channelObj = digital_rf_hdf5.write_hdf5_channel('/tmp/benchmark/junk1', 'i2', FILE_SAMPLES, FILES_PER_DIR, start_global_index, SAMPLE_RATE, 'Fake_uuid', 0, True)
t = time.time()
for i in range(N_WRITES):
    channelObj.rf_write(data_int16)
channelObj.close()
seconds = time.time() - t
speedMB = (N_WRITES*4*WRITE_BLOCK_SIZE)/(1.0E6*seconds)
print('Total time %i seconds, speed %1.2f MB/s' % (int(seconds), speedMB))


print("Test 2 - simple single write to multiple files, compress to level 9, with checksum - channel 2")
os.system("rm -rf /tmp/benchmark/junk2 ; mkdir /tmp/benchmark/junk2")
print("Start writing")
channelObj = digital_rf_hdf5.write_hdf5_channel('/tmp/benchmark/junk2', 'i2', FILE_SAMPLES, FILES_PER_DIR, start_global_index, SAMPLE_RATE, 'Fake_uuid', 9, True)
t = time.time()
for i in range(N_WRITES):
    channelObj.rf_write(data_int16)
channelObj.close()
seconds = time.time() - t
speedMB = (N_WRITES*4*WRITE_BLOCK_SIZE)/(1.0E6*seconds)
print('Total time %i seconds, speed %1.2f MB/s' % (int(seconds), speedMB))
