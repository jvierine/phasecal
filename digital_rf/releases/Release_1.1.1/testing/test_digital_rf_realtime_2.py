"""test_digital_rf_realtime_2.py is a script to test the digital_rf_hdf5 module
under realtime conditions

In this realtime test, RF files are added to last subdirectory, and a subdirectory is also added.

$Id: test_digital_rf_realtime_2.py 525 2014-08-04 18:22:16Z brideout $
"""
usage = 'python test_digital_rf_realtime_2.py <0 for fast metadata, 1 for full metadata>'

# standard python imports
import os, os.path, sys
import datetime, time
import traceback
import glob
import shutil

# third party imports
import numpy
import h5py

# Millstone imports
import digital_rf_hdf5

if len(sys.argv) != 2:
    print(usage)
    print('1 or 0 required')
    sys.exit(-1)
    
try:
    update = int(sys.argv[1])
    if update not in (1,0):
        raise ValueError, ''
except:
    print(usage)
    print('1 or 0 required')
    sys.exit(-1)
    
if update == 0:
    update = False
else:
    update = True

# base data
num_subchannels = 4
base_data = []
for i in range(100):
    if i < 50:
        base_data.append([2,3]*num_subchannels)
    else:
        base_data.append([4,6]*num_subchannels)  
                                                    
        
# constants
sample_rate = 1.0E2
files_per_directory = 10
# start 2014-03-09 12:30:30 plus one sample
start_global_index = (1394368230 * sample_rate) + 1

# set up top level directory
os.system("rm -rf /tmp/hdf5 ; mkdir /tmp/hdf5");
        
print("Test 0 - simple single write to multiple files, no compress, no checksum - channel 0");
os.system("rm -rf /tmp/hdf5/junk0 ; mkdir /tmp/hdf5/junk0");
data_object = digital_rf_hdf5.write_hdf5_channel("/tmp/hdf5/junk0", 'i4', 40, files_per_directory, start_global_index,
                                                 sample_rate, "FAKE_UUID_0", 0, False, True, num_subchannels=num_subchannels);
data = numpy.array(base_data, numpy.int32)
for i in range(10):
    result = data_object.rf_write(data);
data_object.close();
print("done write");

# set up fake realtime data by copying files
os.system('rm -rf /tmp/hdf52')
os.system('mkdir /tmp/hdf52')
os.system('mkdir /tmp/hdf52/junk1')
os.system('cp -r /tmp/hdf5/junk0/2014-03-09T12-30-30 /tmp/hdf52/junk1/')
os.system('mkdir /tmp/hdf52/junk1/2014-03-09T12-30-34')
files = glob.glob('/tmp/hdf5/junk0/2014-03-09T12-30-34/*')
files.sort()
for thisFile in files[:5]:
    shutil.copy(thisFile, '/tmp/hdf52/junk1/2014-03-09T12-30-34/')

# sleep for 4 seconds to make sure system knows all files closed
time.sleep(4)

# read
testReadObj = digital_rf_hdf5.read_hdf5(['/tmp/hdf52'], load_all_metadata=False)
channels = testReadObj.get_channels()
print(channels)

print('working on junk1')
start_index, end_index = testReadObj.get_bounds('junk1')
print(('bounds are: ', start_index, end_index))
cont_data_arr = testReadObj.get_continuous_blocks(start_index, end_index, 'junk1')
print(('continuous data is ', cont_data_arr))
result = testReadObj.read_vector_raw(cont_data_arr[0][0], cont_data_arr[0][1], 'junk1')
print('got %i samples' % (len(result)))


# simulate realtime update
time.sleep(5)
for thisFile in files[5:]:
    shutil.copy(thisFile, '/tmp/hdf52/junk1/2014-03-09T12-30-34/')
    cmd = 'touch %s' % (os.path.join('/tmp/hdf52/junk1/2014-03-09T12-30-34/', os.path.basename(thisFile)))
    os.system(cmd)
os.system('cp -r /tmp/hdf5/junk0/2014-03-09T12-30-38 /tmp/hdf52/junk1/')
print('new data added')
time.sleep(5)
testReadObj.reload()
start_index, end_index = testReadObj.get_bounds('junk1')
print(('bounds are: ', start_index, end_index))
cont_data_arr = testReadObj.get_continuous_blocks(start_index, end_index, 'junk1')
print(('continuous data is ', cont_data_arr))
result = testReadObj.read_vector_raw(cont_data_arr[0][0], cont_data_arr[0][1]+200, 'junk1')
print('got %i samples' % (len(result)))


print('Overall test passed')
