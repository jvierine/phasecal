"""test_millstone_data_access.py is a test of using digital_rf_hdf5.read_hdf5 for 
reading small amounts of data from a large data set.

$Id: test_millstone_data_access.py 394 2014-05-23 17:25:57Z brideout $
"""

# constants - modify as needed
top_level_dir = '/data0/ringbuffer_h5'
test_chan = 'misa-l'

import digital_rf_hdf5
import time
t = time.time()
o = digital_rf_hdf5.read_hdf5(top_level_dir)
print('overall init took %f' % (time.time() - t))
print('The channels are <%s>' % (str(o.get_channels())))
b = o.get_bounds(test_chan)
print('Bounds of channel %s are <%s>' % (test_chan, str(b)))

print('Test getting small subsets of get_continuous_blocks')
for i in range(5):

    start_index = b[0] + long(0.2*i*(b[1]-b[0]))
    end_index = start_index + 10000000
    t = time.time()
    print('continuous blocks between %i and %i are <%s>' % (start_index, end_index, str(o.get_continuous_blocks(start_index,end_index,test_chan))))
    print('took %f' % (time.time() - t))

print('Test getting small subsets of read_vector at different places')
for i in range(5):

    start_index = b[0] + long(0.03*i*(b[1]-b[0]))
    end_index = start_index + 10000000
    t = time.time()
    vector = o.read_vector(start_index,end_index-start_index,test_chan)
    print('read_vector between %i and %i are len=%i, value=<%s>' % (start_index, end_index, len(vector), str(vector)))
    print('took %f' % (time.time() - t))

print('Test getting small subsets of contiguous read_vector (so assuming no reload needed after the first')
start_index = b[0] + long(0.032*(b[1]-b[0]))
for i in range(5):
    end_index = start_index + 10000000
    t = time.time()
    vector = o.read_vector(start_index,end_index-start_index,test_chan)
    print('read_vector between %i and %i are len=%i, value=<%s>' % (start_index, end_index, len(vector), str(vector)))
    print('took %f' % (time.time() - t))
    start_index = end_index

t = time.time()
o.reload()
print('reload took %f' % (time.time() - t))
b = o.get_bounds(test_chan)
print('Bounds of channel %s are now <%s>' % (test_chan, str(b)))

metadata_dict = o.get_rf_file_metadata(test_chan)
print('metadata rf file dict:')
keys = metadata_dict.keys()
keys.sort()
for key in keys:
    print('%s: %s' % (str(key), str(metadata_dict[key])))
    
metadata_file = o.get_metadata(test_chan)
