"""py_rf_write_hdf5.py is a module that exposes the capabilities of the C rf_write_hdf5 library

$Id: py_rf_write_hdf5.py 281 2014-03-06 16:43:29Z brideout $
"""

# standard python imports
import os, os.path, sys
import types

# third party imports
import numpy

# Millstone imports
import _py_rf_write_hdf5  # c extension

class rf_write_hdf5_channel:
    """The class rf_write_hdf5_channel is an object used to write rf data to Hdf5 files as specified
    in the http://www.haystack.mit.edu/pipermail/rapid-dev/2014-February/000273.html email thread.
    
    Usage - this class has only two public methods: an init, and a write where a user passes in numpy arrays.
    It also has the following attributes:
        total_samples_written - total number of samples written in channel (does not include gaps)
        total_gaps - total number of samples left as default in channel
        next_available_sample - the index of the next sample available for writing.  Equal to 
            (total_samples_written + total_gaps)
    """
    
    def __init__(self, directory, dtype_str, max_samples, uuid_str, compression_level=0, checksum=False):
        """__init__ creates an rf_write_hdf5_channel
        
        Inputs:
            directory - the directory where this channel is to be written.  Must already exist and be writable
            
            dtype_str - format of numpy data in string format.  String is format as passed into numpy.dtype().
                For example, numpy.dtype('>i4').  For now accepts any legal byte-order character (No character means native),
                and one of 'i1', 'u1', 'i2', 'u2', 'i4', 'u4', 'i8', 'u8', 'f', or 'd'.
                
            max_samples - number of samples in each Hdf5 file
            
            uuid_str - uuid string that will tie the data files to the Hdf5 metadata
            
            compression_level - 0 for no compression (default), 1-9 for varying levels of gzip compression (1 least compression, least CPU,
                9 most compression, most CPU)
                
            checksum - if True, use Hdf5 checksum capability, if False (default) no checksum.
        """
        if not os.access(directory, os.W_OK):
            raise IOError, 'Directory %s does not exist or is not writable' % (directory)
        self.directory = directory
        
        # use numpy to get all needed info about this datatype
        self.dtype = numpy.dtype(dtype_str)
        self.byteorder = self.dtype.byteorder
        if self.byteorder == '=':
            # simplify C code by convertion here
            if sys.byteorder == 'big':
                self.byteorder = '>'
            else:
                self.byteorder = '<'
        
        if max_samples < 1 or max_samples > 1E11:
            raise ValueError, 'max_samples must be between 1 and 1E11, not %s' % (str(max_samples))
        self.max_samples = long(max_samples)
        
        if type(uuid_str) != types.StringType:
            raise ValueError, 'uuid_str must be StringType, not %s' % (str(type(uuid_str)))
        self.uuid = str(uuid_str)
        
        if compression_level not in range(10):
            raise ValueError, 'compression_level must be 0-9, not %s' % (str(compression_level))
        self.compression_level = compression_level
        
        self.checksum = bool(checksum)
        
        # call the underlying C extension, which will call the C init method
        self._channelObj = _py_rf_write_hdf5.init(directory, self.byteorder, self.dtype.char, self.dtype.itemsize,
                                                  self.max_samples, uuid_str, compression_level, int(self.checksum))
        
        # set the next available sample to write at
        self.next_avail_sample = long(0)
        self.total_samples_written = long(0)
        self.total_gaps = long(0)
        
        
        
    def iq_write(self, arr, next_sample=None):
        """iq_write writes a numpy array to Hdf5.
        
        Inputs - arr - numpy array of data of size Nx2, where I data goes in the first column, and Q data in the second.
                    Error will be raised if its not the same data type set in init.
                next_sample - global index of next sample to write to.  Default is self.next_avail_sample.  Error raised
                    if next_sample < self.next_avail_sample
                    
        Returns: self.next_avail_sample
        """
        # verify input arguments
        if arr.shape[1] != 2:
            raise ValueError, 'array to write must have shape Nx2, not %s' % (str(arr.shape))
        
        if arr.dtype != self.dtype:
            raise ValueError, 'arr has dtype %s, but dtype set in init was %s' % (str(arr.dtype), str(self.dtype))
        
        if next_sample == None:
            next_sample = self.next_avail_sample
        else:
            next_sample = long(next_sample)
        if next_sample < self.next_avail_sample:
            raise ValueError, 'Trying to write at sample %i, but next available sample is %i' % (next_sample, self.next_avail_sample)
        
        vector_length = int(arr.shape[0])
        
        result = _py_rf_write_hdf5.iq_write(self._channelObj, arr, next_sample, vector_length)
        
        # update index attributes
        self.total_gaps += next_sample - self.next_avail_sample
        self.total_samples_written += vector_length
        self.next_avail_sample += (next_sample - self.next_avail_sample) + vector_length
        
        
    def __del__(self):
        """__del__ frees the C object before calling the normal __del__
        """
        if hasattr(self, '_channelObj'):
            _py_rf_write_hdf5.free(self._channelObj)
            
        
        