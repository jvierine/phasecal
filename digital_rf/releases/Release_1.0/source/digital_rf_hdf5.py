"""digital_rf_hdf5.py is a module that allows python to write and read Hdf5 digital rf data.

It uses h5py to read, and exposes the capabilities of the C rf_write_hdf5 library to write.

It has two classes:
    write_hdf5_channel
    read_hdf5
    
    

$Id: digital_rf_hdf5.py 396 2014-05-23 17:38:45Z brideout $
"""

# standard python imports
import os, os.path, sys
import types
import glob
import traceback
import datetime, time

# third party imports
import numpy
import h5py

# Millstone imports
import _py_rf_write_hdf5  # c extension

def get_unix_time(unix_sample_index, sample_rate):
    """get_unix_time returns a tuple of (datetime, picosecond) given an input unix_sample_index and sample rate
    
    The returned datetime will contain microsecond precision.  Picosecond is also returned for users wanting
    greater precision than is available in the datatime object, where picoseond is the number of picoseconds 
    since the last second in the returned datetime.
    
    Inputs:
        unix_sample_index - number of samples at given sample rate since UT midnight 1970-01-01 (to be precise,
        the unix time at unix second (unix_sample_index / sample_rate).
        
        sample_rate - sample rate in Hz.
        
    Returned value will be precise to the picosecond if sample rate is evenly divisible by 1.0 (no fractional part),
    because in that case integer calculations are used throughout.
    
    Returns (datetime, picosecond since last second)
    """
    year, month, day, hour, minute, second, picosecond = _py_rf_write_hdf5.get_unix_time(unix_sample_index, sample_rate)
    dt = datetime.datetime(year, month, day, hour, minute, second, microsecond=long(picosecond/1.0E6))
    return((dt, picosecond))


class write_hdf5_channel:
    """The class write_hdf5_channel is an object used to write rf data to Hdf5 files as specified
    in the http://www.haystack.mit.edu/pipermail/rapid-dev/2014-February/000273.html email thread.
    """
    
    def __init__(self, directory, dtype_str, samples_per_file, files_per_directory, start_global_index, sample_rate, uuid_str,
                 compression_level=0, checksum=False, is_complex=True, num_subchannels=1, marching_periods=True):
        """__init__ creates an write_hdf5_channel
        
        Inputs:
            directory - the directory where this channel is to be written.  Must already exist and be writable
            
            dtype_str - format of numpy data in string format.  String is format as passed into numpy.dtype().
                For example, numpy.dtype('>i4').  For now accepts any legal byte-order character (No character means native),
                and one of 'i1', 'u1', 'i2', 'u2', 'i4', 'u4', 'i8', 'u8', 'f', or 'd'.
                
            samples_per_file - number of samples in each Hdf5 file
            
            files_per_directory - number files per subdirectory in form YYYY-MM-DDTHH:MM:SS before new subdirectory created
            
            start_global_index - the start time of the first sample in units of (unix_timestamp * sample_rate)
            
            sample_rate - sample rate in Hz
            
            uuid_str - uuid string that will tie the data files to the Hdf5 metadata
            
            compression_level - 0 for no compression (default), 1-9 for varying levels of gzip compression (1 least compression, least CPU,
                9 most compression, most CPU)
                
            checksum - if True, use Hdf5 checksum capability, if False (default) no checksum.
            
            is_complex - if True (the default) data is IQ. If false, each sample has a single value.
            
            num_subchannels - number of subchannels to write simultaneously.  Default is 1.
            
            marching_periods - if True, have matching periods written to stdout when writing. False - do not.
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
        
        if samples_per_file < 1 or samples_per_file > 1E11:
            raise ValueError, 'samples_per_file must be between 1 and 1E11, not %s' % (str(samples_per_file))
        self.samples_per_file = long(samples_per_file)
        
        if files_per_directory < 1:
            raise ValueError, 'files_per_directory must be positive, not %s' % (str(files_per_directory))
        self.files_per_directory = long(files_per_directory)
        
        if start_global_index < 0:
            raise ValueError, 'start_global_index cannot be negative (%s)' % (str(start_global_index))
        self.start_global_index = long(start_global_index)
        
        if sample_rate <= 0.0:
            raise ValueError, 'sample_rate must be positive, not %s' % (str(sample_rate))
        self.sample_rate = float(sample_rate)
        
        if type(uuid_str) != types.StringType:
            raise ValueError, 'uuid_str must be StringType, not %s' % (str(type(uuid_str)))
        self.uuid = str(uuid_str)
        
        if compression_level not in range(10):
            raise ValueError, 'compression_level must be 0-9, not %s' % (str(compression_level))
        self.compression_level = compression_level
        
        self.checksum = bool(checksum)
        self.is_complex = bool(is_complex)
        
        if num_subchannels < 1:
            raise ValueError, 'Number of subchannels must be at least one, not %i' % (num_subchannels)
        self.num_subchannels = int(num_subchannels)
        
        if marching_periods:
            use_marching_periods = 1
        else:
            use_marching_periods = 0
        
        # call the underlying C extension, which will call the C init method
        self._channelObj = _py_rf_write_hdf5.init(directory, self.byteorder, self.dtype.char, self.dtype.itemsize,
                                                  self.samples_per_file, self.files_per_directory, self.start_global_index,
                                                  self.sample_rate, uuid_str, compression_level, int(self.checksum),
                                                  int(self.is_complex), self.num_subchannels, use_marching_periods)
        
        if not self._channelObj:
            raise ValueError, 'Failed to create write_hdf5_channel'
        
        # set the next available sample to write at
        self._next_avail_sample = long(0)
        self._total_samples_written = long(0)
        self._total_gap_samples = long(0)
        
        
        
    def rf_write(self, arr, next_sample=None):
        """rf_write writes a numpy array to Hdf5.   Must have the same number of subchannels as declared in init.
        For single valued data, number of columns == number of subchannels.  For complex data, there are two types
        of input arrays that are allowed:
            1. An array without column names with number of columns = 2*num_subchannels.  I/Q are assumed to be interleaved.
            2. A structured array with column names r and i, as stored in the Hdf5 file.  Then the shape will
                be N * num_subchannels, because numpy considered the r/i data as one piece of data.
                
        Here's an example of one way to create a structured numpy array with complex data with dtype int16:
        
        arr_data = numpy.ones((num_rows, num_subchannels), dtype=[('r', numpy.int16), ('i', numpy.int16)])
        for i in range(num_subchannels):
            for j in range(num_rows):
                arr_data[j,i]['r'] = 2
                arr_data[j,i]['i'] = 3
                
        The same data could be passed in via the array created as:
        
        arr_data = numpy.ones((num_rows, num_subchannels*2), dtype=numpy.int16)
        
        Inputs - arr - numpy array of data of size described above if complex, and size num_rows if not. 
                 Error will be raised if its not the same data type set in init.
                 
                next_sample - global index of next sample to write to.  Default is self._next_avail_sample.  Error raised
                    if next_sample < self._next_avail_sample
                    
        Returns: self._next_avail_sample
        """
        # verify input arr argument
        self._verify_input(arr)
        
        if next_sample == None:
            next_sample = self._next_avail_sample
        else:
            next_sample = long(next_sample)
        if next_sample < self._next_avail_sample:
            raise ValueError, 'Trying to write at sample %i, but next available sample is %i' % (next_sample, self._next_avail_sample)
        
        vector_length = int(arr.shape[0])
        
        result = _py_rf_write_hdf5.rf_write(self._channelObj, arr, next_sample)
        
        # update index attributes
        self._total_gap_samples += next_sample - self._next_avail_sample
        self._total_samples_written += vector_length
        self._next_avail_sample += (next_sample - self._next_avail_sample) + vector_length
        
        
    def rf_write_blocks(self, arr, global_sample_arr, block_sample_arr):
        """rf_write_blocks writes a data with interleaved gaps to Hdf5 files
        
        Inputs - arr - numpy array of data. See rf_write for a complete description of allowed forms.
                 
                global_sample_arr an array len < N, > 0 of type numpy.uint64 that sets the global sample index for
                     each continuous block of data in arr.  Must be increasing, and first value must be >= self._next_avail_sample
                     of ValueError raised.
                
                block_sample_arr an array len = len(global_sample_arr) of type numpy.uint64.  Values are the index into arr of each block.
                    Values must be < len(arr).  First value must be zero.  Increments between value must be > 0 and less than the 
                    corresponding increment in global_sample_arr
                    
        Returns: self._next_avail_sample
        """
        # verify input arr argument
        self._verify_input(arr)
        
        if global_sample_arr[0] < self._next_avail_sample:
            raise ValueError, 'first value in global_sample_arr must be at least %i, not %i' % (self._next_avail_sample,
                                                                                                global_sample_arr[0])
            
        if block_sample_arr.dtype != numpy.uint64:
            raise ValueError, 'block_sample_arr has dtype %s, but needs to have numpy.uint64' % (str(block_sample_arr.dtype))
        
        if block_sample_arr[0] != 0:
            raise ValueError, 'first value in block_sample_arr must be 0, not %i' % (block_sample_arr[0])
        
        if len(global_sample_arr) != len(block_sample_arr):
            raise ValueError, 'len of global_sample_arr (%i) must equal len of block_sample_arr (%i)' % (len(global_sample_arr), 
                                                                                                         len(block_sample_arr))
        
        # data passed initial tests, try to write
        result = _py_rf_write_hdf5.rf_block_write(self._channelObj, arr, global_sample_arr, block_sample_arr)
        
        # update index attributes
        self._total_gap_samples += global_sample_arr[0] - self._next_avail_sample # potential gap between writes
        self._total_gap_samples += (global_sample_arr[-1] - global_sample_arr[0]) - (block_sample_arr[-1]) # gaps within write
        self._total_samples_written += len(arr)
        self._next_avail_sample = global_sample_arr[-1] + (len(arr) - block_sample_arr[-1])
        
        
    def get_total_samples_written(self):
        """get_total_samples_written returns the total number of samples written in channel (does not include gaps)
        """
        return(self._total_samples_written)
    
    
    def get_next_available_sample(self):
        """get_next_available_sample returns the index of the next sample available for writing.  Equal to 
            (total_samples_written + total_gap_samples)
        """
        return(self.next_available_sample)
    
    
    def get_total_gap_samples(self):
        """get_total_gap_samples returns the total number of samples left as default in channel
        """
        return(self._total_gap_samples)
        
        
    def close(self):
        """close frees the C object and closes the last Hdf5 file
        """
        _py_rf_write_hdf5.free(self._channelObj)
        
        
    def _verify_input(self, arr):
        """_verify_input checks for valid and consistent input arrays being passed in to write commands.
        Throws ValueError if invalid.
        
        Input:
            arr - see rf_write method for a complete description of allowed values
        """
        if self.is_complex:
            # there are two allowed ways to pass in complex data - see which one used
            if arr.dtype.names != None:
                # this must be be r/i format:
                for name in ('r', 'i'):
                    if name not in arr.dtype.names:
                        raise ValueError, 'column names must be r and i, not %s' % (str(arr.dtype.names))
                    if arr.dtype[name] != self.dtype:
                        raise ValueError, 'column %s must have dtype %s, not %s' % (name, str(self.dtype), str(arr.dtype[name]))
                if len(arr.dtype.names) != 2:
                    raise ValueError, 'column names must be only r and i, not %s' % (str(arr.dtype.names))
                if arr.shape[1] != self.num_subchannels:
                    raise ValueError, 'complex array in r/i form must have shape N x num_subchannels, not %s' % (str(arr.shape))
            else:
                if arr.shape[1] != 2*self.num_subchannels:
                    raise ValueError, 'complex array in flat form must have shape N x 2*num_subchannels, not %s' % (str(arr.shape))
                if arr.dtype != self.dtype:
                    raise ValueError, 'arr has dtype %s, but dtype set in init was %s' % (str(arr.dtype), str(self.dtype))
        
        else: # single value checks
            if arr.dtype != self.dtype:
                raise ValueError, 'arr has dtype %s, but dtype set in init was %s' % (str(arr.dtype), str(self.dtype))
            if len(arr.shape) != 1:
                raise ValueError, 'single valued array must just have one dimension, not shape %s' % (str(arr.shape))
        
            
            
class read_hdf5:
    """The class read_hdf5 is an object used to read rf data from Hdf5 files as specified
    in the http://www.haystack.mit.edu/pipermail/rapid-dev/2014-February/000273.html email thread.
    
    This class allows random access to the rf data.
    
    """
    def __init__(self, top_level_directory_arg):
        """__init__ will verify the data in top_level_directory_arg is as expected.  It will analyze metadata about all
        Hdf5 files so that other methods can return quickly
        
        Inputs:
            directory_arg - either a single top level, directory, or a list.  A directory can be a file system path or a url,
            where the url points to a top level directory.
            
        A top level directory must contain <channel_name>/<YYYY-MM-DDTHH:MM:SS/rf@<unix_seconds>.<%03i milliseconds>.h5
        
        If more than one top level directory contains the same channel_name subdirectory, this is considered the same channel.  An error
        is raised if their sample rates differ, or if their time periods overlap.
        
        This method will create the following attributes:
        
        self._top_level_dir_dict - a dictionary with keys = top_level_directory string, value = access mode (eg, 'local', 'file', or 'http')
            This attribute is static, that is, it is not updated when self.reload() called
            
        self._channel_dict - a dictionary with keys = channel_name, and value is a _channel_metadata object.
        
        """
        
        # first, make top_level_directory_arg a list if a string
        if type(top_level_directory_arg) == types.StringType:
            top_level_arg = [top_level_directory_arg]
        else:
            top_level_arg = top_level_directory_arg
        
        # create static attribute self._top_level_dir_dict
        self._top_level_dir_dict = {}
        for top_level_directory in top_level_arg:
            if top_level_directory[0:7] == 'file://':
                self._top_level_dir_dict[top_level_directory] = 'file'
            elif top_level_directory[0:7] == 'http://':
                self._top_level_dir_dict[top_level_directory] = 'http'
            else:
                self._top_level_dir_dict[top_level_directory] = 'local'
                
        self._channel_dict = {}
        self._channel_dir_metadata = {}
        
        self.reload()
        
        
    def reload(self, start_sample=None, sample_length=None, channel=None):
        """reload updates the attributes self._channel_location_dict and self._channel_dict.  If start_sample and/or
        sample_length == None, then only get high level metadata.  If start_sample and sample_length given, then only 
        updated detailed metadata as needed.
        
            Inputs:
                start_sample - sample number if units of samples since 1970-01-01.  If None, only refresh high level metadata.
                sample_length - number of samples.  If given, only refresh low-level metadata near sample extent. If None, 
                    only refresh high level metadata.
                channel - channel to update.  If None, only refresh high level metadata.
        """
        if start_sample == None or sample_length == None or channel == None:
            high_level = True
        else:
            high_level = False
            
        if high_level:
            self._high_level_reload()
        else:
            self._low_level_reload(start_sample, sample_length, channel)
            
                    
                
    def get_channels(self):
        """get_channels returns a alphabetically sorted list of channels in this read_hdf5 object
        """
        channels = self._channel_dict.keys()
        channels.sort()
        return(channels)
    
    
    def get_bounds(self, channel_name):
        """get_bounds returns a tuple of (first_unix_sample, last_unix_sample) for a given channel name
        """
        channel_metadata = self._channel_dict[channel_name]
        return((long(channel_metadata.unix_start_sample), 
                long(channel_metadata.unix_start_sample + channel_metadata.sample_extent)))
        
        
    
    def get_rf_file_metadata(self, channel_name):
        """get_rf_file_metadata returns a dictionary of metadata found as attributes in the Hdf5 file /rf_data
        dataset for the given channel name.
        """
        return(self._channel_dict[channel_name].metadata_dict)
    
    
    
    def get_metadata(self, channel_name, timestamp=None):
        """get_metadata returns a h5py.File object pointing to the metadata*.h5 file at the top level of the 
        channel directory.  The user is responsible for closing that file when done.  If timestamp == None,
        the latest metadata*.h5 will be returned.  Otherwise it will open the earliest metadata file with timestamp
        greater than or equal to timestamp
        """
        # first, get a sorted list of all metadata*.hf files
        metadata_file_list = []
        metadata_basename_list = [] # to make sure there are no repeated basenames
        channel = self._channel_dict[channel_name]
        for top_level_dir_obj in channel.top_level_dir_meta_list:
            metadata_files = glob.glob(os.path.join(top_level_dir_obj.top_level_dir, 
                                                    top_level_dir_obj.channel_name, 'metadata@*.h5'))
            metadata_files.sort()
            for metadata_file in metadata_files:
                basename = os.path.basename(metadata_file)
                if basename in metadata_basename_list:
                    raise IOError, 'found repeated metadata file names in channel %s' % (channel_name)
                metadata_basename_list.append(basename)
                metadata_file_list.append(metadata_file)
                
        if len(metadata_file_list) == 0:
            raise IOError, 'No metadata files found in channel %s' % (channel_name)
                
        # open right metadata file
        if timestamp == None:
            return(h5py.File(metadata_file_list[-1], 'r'))
        else:
            rightFile = None
            for metadata_file in metadata_file_list:
                basename = os.path.basename(metadata_file)
                this_timestamp = long(basename[len('metadata@'):basename.find('.')])
                if this_timestamp <= timestamp:
                    rightFile = metadata_file
                else:
                    break
            if rightFile == None:
                raise IOError, 'All metadata files found in channel %s after timestamp' % (channel_name, timestamp)
            return(h5py.File(rightFile, 'r'))
            
        
        
        
        
    def get_continuous_blocks(self, start_unix_sample, stop_unix_sample, channel_name):
        """get_continuous_blocks returns a numpy array of dtype u64 and shape (N,2) where the first
        column represents the unix_sample of a continuous block of data, and the second column represents the
        number of samples in that continuous block.  Only samples between (start_unix_sample, stop_unix_sample)
        inclusive will be returned.
        
        Calls the private method _get_continuous_blocks.  If that raises a _MissingMetadata exception, calls
        reload to get missing metadata, and then retries _get_continuous_blocks.
        
        Returns IOError if no blocks found
        
        Inputs:
            start_unix_sample, stop_unix_sample - only samples between (start_unix_sample, stop_unix_sample)
                inclusive will be returned.  Value of both are samples since 1970-01-01
                
            channel_name - channel to examine
        """
        try:
            return(self._get_continuous_blocks(start_unix_sample, stop_unix_sample, channel_name))
        except _MissingMetadata:
            self.reload(start_unix_sample, stop_unix_sample-start_unix_sample, channel_name)
            # try again
            return(self._get_continuous_blocks(start_unix_sample, stop_unix_sample, channel_name))
                
        
    
    def read_vector(self, unix_sample, vector_length, channel_name):
        """read_vector returns a numpy array of dim(up to num_samples, num_subchannels) of the dtype in the Hdf5 files.
        
        If complex data, real and imag data will have names 'r' and 'i' if underlying data are integers 
        or be numpy complex data type if underlying data floats.
        
        Calls private method _read_vector.  If that method raises a _MissingMetadata exception, calls
        reload, then tries again.
        
        Inputs:
            unix_sample - the number of samples since 1970-01-01 at start of data
            
            vector_length - the number of continuous samples to include
            
            channel_name - the channel name to use
        
        This method will raise an IOError error if the returned vector would include any missing data. 
        It will also raise an IOError is any of the files needed to read the data have been deleted.  
        This is possible because metadata on which this call is based might be out of date.
        """
        try:
            return(self._read_vector(unix_sample, vector_length, channel_name))
        except _MissingMetadata:
            self.reload(unix_sample, vector_length, channel_name)
            return(self._read_vector(unix_sample, vector_length, channel_name))
        
        
    def read_vector_c81d(self, unix_sample, vector_length, channel_name, subchannel=0):
        """read_vector_c81d returns a numpy vector of complex8 type, no matter the dtype of the Hdf5 file
        or the number of channels. Error thrown if subchannel doesn't exist
        
        Calls private method _read_vector.  If that method raises a _MissingMetadata exception, calls
        reload, then tries again.
        
        Inputs:
            unix_sample - the number of samples since 1970-01-01 at start of data
            
            vector_length - the number of continuous samples to include
            
            channel_name - the channel name to use
            
            subchannel - which subchannel to use.  Default is 0 (first)
        
        This method will raise an IOError error if the returned vector would include any missing data. 
        It will also raise an IOError is any of the files needed to read the data have been deleted.  
        This is possible because metadata on which this call is based might be out of date.
        """
        try:
            z = self._read_vector(unix_sample, vector_length, channel_name)
        except _MissingMetadata:
            self.reload(unix_sample, vector_length, channel_name)
            z = self._read_vector(unix_sample, vector_length, channel_name)
            
        if z.shape[1] < subchannel + 1:
            raise ValueError, 'Returned data has only %i subchannels, does not have subchannel %i' % (z.shape[1], subchannel)
        
        if z.dtype == numpy.complex64:
            return(z[:,subchannel])
        elif z.dtype in (numpy.complex128, numpy.complex256):
            return(numpy.array(z[:,subchannel], dtype=numpy.complex64))
        
        slice = (z[:,subchannel])
        if not hasattr(slice.dtype, 'names'):
            raise ValueError, 'Single valued channels cannot be cast to complex'
        elif slice.dtype.names == None:
            raise ValueError, 'Single valued channels cannot be cast to complex'
        slice = numpy.array(slice['r'] + slice['i']*1.0j, dtype=numpy.complex64)
        return(slice)
        
        
    def _read_vector(self, unix_sample, vector_length, channel_name):
        """_read_vector is a private method returns a numpy array of dim(up to num_samples, 2) of the dtype in the Hdf5 files.
        
        Raises _MissingMetadata error is encounters not up to data metadata.
        
        Inputs:
            unix_sample - the number of samples since 1970-01-01 at start of data
            
            vector_length - the number of continuous samples to include
            
            channel_name - the channel name to use
        
        This method will raise an IOError error if the returned vector would include any missing data. 
        It will also raise an IOError is any of the files needed to read the data have been deleted.  
        This is possible because metadata on which this call is based might be out of date.
        """
        if vector_length < 1:
            raise IOError, 'Number of samples requested must be greater than 0, not %i' % (vector_length)
        
        # make sure everything is a long
        unix_sample = long(unix_sample)
        vector_length = long(vector_length)
        
        channel_metadata = self._channel_dict[channel_name]
        
        ret_array = None
        first_unix_sample = None
        
        for top_level_dir in channel_metadata.top_level_dir_meta_list:
            if top_level_dir.unix_start_sample + top_level_dir.sample_extent < unix_sample:
                # this top level dir is too early
                continue
            if unix_sample + vector_length < top_level_dir.unix_start_sample:
                # this top level dir is too late
                continue
            this_array, unix_sample = top_level_dir.get_continuous_vector(max(unix_sample, top_level_dir.unix_start_sample),
                                                                          min(unix_sample + vector_length, top_level_dir.unix_start_sample + top_level_dir.sample_extent))
            ret_array = self._combine_continuous_vectors(ret_array, this_array, first_unix_sample, unix_sample)
            if first_unix_sample == None:
                first_unix_sample = unix_sample
            
        if ret_array == None:
            raise IOError, 'No data found for channel %s between %i and %i' % (channel_name, unix_sample, 
                                                                               unix_sample + vector_length)
                
        return(ret_array)
    
    
    def _get_continuous_blocks(self, start_unix_sample, stop_unix_sample, channel_name):
        """_get_continuous_blocks is a private method that returns a numpy array of dtype u64 and shape (N,2) where the first
        column represents the unix_sample of a continuous block of data, and the second column represents the
        number of samples in that continuous block.  Only samples between (start_unix_sample, stop_unix_sample)
        inclusive will be returned.
        
        Raises _MissingMetadata exception if reload needs to be called. 
        
        Returns IOError if no blocks found
        
        Inputs:
            start_unix_sample, stop_unix_sample - only samples between (start_unix_sample, stop_unix_sample)
                inclusive will be returned.  Value of both are samples since 1970-01-01
                
            channel_name - channel to examine
        """ 
        channel_metadata = self._channel_dict[channel_name]
        
        ret_array = numpy.array([], dtype=numpy.uint64)
        for top_level_dir in channel_metadata.top_level_dir_meta_list:
            if top_level_dir.unix_start_sample + top_level_dir.sample_extent < start_unix_sample:
                # this top level dir is too early
                continue
            if stop_unix_sample < top_level_dir.unix_start_sample:
                # this top level dir is too late
                continue
            this_array = top_level_dir.get_continuous_blocks(max(start_unix_sample, top_level_dir.unix_start_sample),
                                                             min(stop_unix_sample, top_level_dir.unix_start_sample + top_level_dir.sample_extent))
            ret_array = self._combine_blocks(ret_array, this_array, top_level_dir.samples_per_file)
            
        if len(ret_array) == 0:
            raise IOError, 'No data found for channel %s between %i and %i' % (channel_name, start_unix_sample, 
                                                                               stop_unix_sample)
            
        return(ret_array)
    
    
    def _high_level_reload(self):
        """_high_level_reload updates the attributes self._channel_location_dict and self._channel_dict using only
        high level metadata
        """
        # first update the channel list
        channel_dict = {} # a temporary dict with key = channels, value = list of top level directories where found
        for top_level_dir in self._top_level_dir_dict.keys():
            channels_found = self._get_channels_in_dir(top_level_dir)
            for channel in channels_found:
                channel_name = os.path.basename(channel)
                if channel_dict.has_key(channel_name):
                    channel_dict[channel_name].append(top_level_dir)
                else:
                    channel_dict[channel_name] = [top_level_dir]
                    
        # next throw away any metadata where the entire channel not longer exists
        remove_keys = []
        for channel_name in self._channel_dict.keys():
            if channel_name not in channel_dict.keys():
                # channel no longer exists
                remove_keys.append(channel_name)
        if len(remove_keys):
            for remove_key in remove_keys:
                del self._channel_dict[remove_key] 
                
                    
        # update all channels
        for channel_name in channel_dict.keys():
            if not self._channel_dict.has_key(channel_name):
                # a new channel is found - create it
                top_level_dir_metadata_list = []
                for top_level_dir in channel_dict[channel_name]:
                    new_top_level_metaddata = _top_level_dir_metadata(top_level_dir, channel_name,
                                                                      self._top_level_dir_dict[top_level_dir])
                    top_level_dir_metadata_list.append(new_top_level_metaddata)
                top_level_dir_metadata_list.sort()
                new_channel_metadata = _channel_metadata(channel_name, top_level_dir_meta_list = top_level_dir_metadata_list)
                new_channel_metadata.update()
                self._channel_dict[channel_name] = new_channel_metadata
                
            else:
                # handle any changes to the top_level_dir list
                chan_obj = self._channel_dict[channel_name] # just to shorten the following code
                found_dirs = [] # a list of directories in _channel_metadata.top_level_dir_meta_list found
                for top_level_dir in channel_dict[channel_name]:
                    found = False
                    for chan_top_dir in chan_obj.top_level_dir_meta_list:
                        if top_level_dir == chan_top_dir.top_level_dir:
                            found_dirs.append(chan_top_dir.top_level_dir)
                            found = True
                            break
                    if not found:
                        # this is a new top level
                        new_top_level_meta = _top_level_dir_metadata(top_level_dir, channel_name,
                                                                     self._top_level_dir_dict[top_level_dir])
                        chan_obj.add_top_level(new_top_level_meta)
                        found_dirs.append(top_level_dir)
                        
                # make sure all top level dirs in chan metadata still exist
                for chan_top_dir in chan_obj.top_level_dir_meta_list:
                    if chan_top_dir.top_level_dir not in found_dirs:
                        # this top level dir no longer has data
                        chan_obj.remove_top_level_metadata(chan_top_dir.top_level_dir)
                        
                chan_obj.update()
                
                
    def _low_level_reload(self, start_sample,  sample_length, channel):
        """_low_level_reload updates the attributes self._channel_location_dict and self._channel_dict only
        between start_sample,  sample_length, but with detailed metadata for a given channel.
        
        Inputs:
            start_sample - sample number if units of samples since 1970-01-01.
            sample_length - number of samples.
            channel - channel to update
        """
        chan_obj = self._channel_dict[channel]
        chan_obj.update(start_sample,  sample_length)
    
    
    
    def _combine_continuous_vectors(self, first_array, second_array, first_start_sample, second_start_sample):
        """_combine_continuous_vectors returns the concatenation of first_array and second_array,  Raises error
        if two vectors are not continuous.
        
        Inputs:
            first_array - first array to combine.  If None, just return second_array
            second_array - second_array to merge at end of first
            first_start_sample - unix_sample of first sample in first_array.  None if first_array == None
            second_start_sample - unix_sample of first sample in second_array
        """
        
        if first_array == None:
            return(second_array)
        
        if len(first_array) != second_start_sample - first_start_sample:
            raise IOError, '_combine_continuous_vectors trying to combine two non-continuous vectors'
        
        return(numpy.concatenate((first_array, second_array)))
    
    
    def _combine_blocks(self, first_array, second_array, samples_per_file):
        """_combine_blocks combines two numpy array of dtype u64 and shape (N,2) where the first
        column represents the unix_sample of a continuous block of data, and the second column represents the
        number of samples in that continuous block. The first row of the second array may or may not be contiguous
        with the last row of the first array.  If it is contiguous, that row will not be included, and the
        number of samples in that first row will instead be added to the last row of first_array. If not contiguous,
        the two arrays are simply concatenated
        """
        if len(first_array) == 0:
            return(second_array)
        is_contiguous = False
        if first_array[-1][0] + first_array[-1][1] > second_array[0][0]:
            raise IOError, 'overlapping data found in top level directories %i %i' % \
                (first_array[-1][0] + first_array[-1][1], second_array[0][0])
        if first_array[-1][0] + first_array[-1][1] == second_array[0][0]:
            is_contiguous = True
        if is_contiguous:
            first_array[-1][1] += second_array[0][1]
            if len(second_array) == 1:
                return(first_array)
            else:
                return(numpy.concatenate([first_array, second_array[1:]]))
        else:
            return(numpy.concatenate([first_array, second_array]))
    
    
    def _get_channels_in_dir(self, top_level_dir):
        """_get_channels_in_dir returns a list of channel names found in top_level_dir
        
        Inputs:
            top_level_dir - string indicating top_level_dir
        """
        # define glob string for sub_directories in form YYYY-MM-DDTHH:MM:SS
        sub_directory_glob = '[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9]:[0-9][0-9]:[0-9][0-9]'
        
        retList = []
        access_mode = self._top_level_dir_dict[top_level_dir]
        # for now only local access
        if access_mode not in ('local'):
            raise ValueError, 'access_mode %s not yet implemented' % (access_mode)
        
        if access_mode == 'local':
            potential_channels = glob.glob(os.path.join(top_level_dir, '*', sub_directory_glob))
            for potential_channel in potential_channels:
                channel_name = os.path.dirname(potential_channel)
                if channel_name not in retList:
                    retList.append(channel_name)
                    
        return(retList)
    
    
    def _needs_updating(self, channel_name, top_level_dir):
        """_needs_updating returns True if channel_name, top_level_dir needs updating, False if not.
        
        Inputs:
            channel_name - channel_name (also directory under top_level_dir)
            top_level_dir - top level directory where channel_name directory is to be found
        """
        access_mode = self._top_level_dir_dict[top_level_dir]
        # for now only local access
        if access_mode not in ('local'):
            raise ValueError, 'access_mode %s not yet implemented' % (access_mode)
        
        recarray = self
        

    
    
class _channel_metadata:
    """The _channel_metadata is a private class to hold and access metadata about a particular digital_rf channel.
    A channel can extend over one of more top level directories.
    """
    
    def __init__(self, channel_name, unix_start_sample = 0, sample_extent = 0, top_level_dir_meta_list = []):
        """__init__ creates a new _channel_metadata object
        
        Inputs:
            channel_name - channel name (name of subdirectory defining this channel)
            unix_start_sample - unix start sample - first sample time in unix timeseconds * sample rate. If default
                0, then unknown
            sample_extent - number of samples between first and last in data
            top_level_dir_meta_list - a time ordered list of _top_level_dir_metadata objects.  Default is empty list
            
        """
        self.channel_name = channel_name
        self.unix_start_sample = long(unix_start_sample)
        self.sample_extent = long(sample_extent)
        self.top_level_dir_meta_list = top_level_dir_meta_list
        self.top_level_dir_meta_list.sort()
        self.metadata_dict = {} # stores all metadata for this _channel_metadata
        self.update()
        
        
    def update(self, start_sample=None, sample_length=None):
        """update will cause this _channel_metadata object to update itself.
        
        If start_sample and/or sample_length == None, then only get high level metadata.  If 
        start_sample and sample_length given, then only updated detailed metadata as needed.
        
        Inputs:
            start_sample - sample number if units of samples since 1970-01-01.  If None, only refresh high level metadata.
            sample_length - number of samples.  If given, only refresh low-level metadata near sample extent. If None, 
                only refresh high level metadata.
        """
        for top_level_meta in self.top_level_dir_meta_list:
            top_level_meta.update(start_sample, sample_length)
            if not self.metadata_dict.has_key('uuid_str'):
                for key in top_level_meta.metadata_dict.keys():
                    self.metadata_dict[key] = top_level_meta.metadata_dict[key]
        self.reset_indices()
    
    
    def add_top_level(self, top_level_dir_meta):
        """add_top_level will add a new _top_level_dir_metadata object to self.top_level_dir_meta_list
        
        Inputs:
            top_level_dir_meta - new _top_level_dir_metadata object to add
        """
        self.top_level_dir_meta_list.append(top_level_dir_meta)
        self.top_level_dir_meta_list.sort()
        self.update()
        
        
    def remove_top_level_metadata(self, top_level_dir):
        """remove_top_level_metadata removes all metadata associated with this channel and one top_level_directory.  Raise
        ValueError is top_level_dir not found
        """
        remove_list = []
        for i, top_level_dir_meta in enumerate(self.top_level_dir_meta_list):
            if top_level_dir == top_level_dir_meta.top_level_dir:
                remove_list.append(i)
                
        if len(remove_list) == 0:
            raise ValueError, 'No directory %s found in this channel' % (top_level_dir)
        elif len(remove_list) > 1:
            raise ValueError, 'More than one directory %s found in this channel' % (top_level_dir)
        
        self.top_level_dir_meta_list.pop(remove_list[0])
        
        self.reset_indicies()
        
        
    def reset_indices(self):
        """reset_indices recalculates self.unix_start_sample and self.sample_extent based on self.top_level_dir_meta_list
        """
        if len(self.top_level_dir_meta_list) == 0:
            # set to unknown
            self.unix_start_sample = long(0)
            self.sample_extent = long(0)
            return
        
        self._verify_non_overlapping_data()
        self.unix_start_sample = long(self.top_level_dir_meta_list[0].unix_start_sample)
        self.sample_extent = long(self.top_level_dir_meta_list[-1].unix_start_sample + self.top_level_dir_meta_list[-1].sample_extent - self.unix_start_sample)
        
        
    def _verify_non_overlapping_data(self):
        """_verify_non_overlapping_data raises an error if any overlapping top level directories found
        """
        for i, record in enumerate(self.top_level_dir_meta_list):
            if i == 0:
                last_unix_start_sample = record.unix_start_sample
                last_sample_extent = record.sample_extent
            else:
                this_unix_start_sample = record.unix_start_sample
                this_sample_extent = record.sample_extent
                if last_unix_start_sample + last_sample_extent > this_unix_start_sample:
                    raise IOError, 'Overlapping samples found in top level dir %s' % (record.top_level_dir)
                last_unix_start_sample = this_unix_start_sample
                last_sample_extent = this_sample_extent
        

class _top_level_dir_metadata:
    """The _top_level_dir_metadata is a private class to hold and access metadata about a particular digital_rf channel in
    a particular top level directory.
    """
    
    def __init__(self, top_level_dir, channel_name, access_mode, unix_start_sample = 0, sample_extent = 0, 
                 samples_per_file=0, sub_directory_recarray=None, sub_directory_dict=None):
        """__init__ creates a new _top_level_dir_metadata
        
        Inputs:
            top_level_dir - full path the top level directory that contains the parent channel_name
            channel_name - the channel_name subdirectory name
            access_mode - string giving access mode (eg, 'local', 'file', or 'http')
            unix_start_sample - unix start sample - first sample time in unix timeseconds * sample rate. If default
                0, then unknown
            sample_extent - number of samples between first and last in data. If default 0, then unknown
            samples_per_file - number of samples per file. If default 0, then unknown
            sub_directory_recarray - a ordered numpy recarray with one row describing summary information about a single
                sub_directory in that channel/top_level_dir named YYYY-MM-DDTHH:MM:SS with the following columns: 
                1) 'subdirectory' - in form YYYY-MM-DDTHH:MM:SS, 
                2) 'unix_start_sample' (for that subdirectory).  May be zero in no detailed metadata for this directory yet
                3) 'sample_extent' eg, number of samples to last sample in that subdirectory.  May be zero if no detailed
                    metadata yet for this subdirectory.
                4) 'file_count' number of Hdf5 data files in directory.  Will be zero if no detailed
                    metadata yet for this subdirectory.
                5) 'last_timestamp' UTC time stamp of latest data file in directory. .  Will be zero if no detailed
                    metadata yet for this subdirectory
                Order is by subdirectory and/or unix_start_sample
            sub_directory_dict - a dictionary with key = sub_directory, value = _sub_directory_metadata object
            
        """
        self.top_level_dir = top_level_dir
        self.channel_name = channel_name
        self.access_mode = access_mode
        self.unix_start_sample = long(unix_start_sample)
        self.sample_extent = long(sample_extent)
        self.samples_per_file = long(samples_per_file)
        self.sub_directory_recarray = sub_directory_recarray
        self.sub_directory_dict = sub_directory_dict
        self.metadata_dict = {} # to be populated by rf file metadata
        
        # data type of sub_directory_array
        self.data_t = numpy.dtype([('subdirectory', numpy.str_, 512), ('unix_start_sample', numpy.uint64, 1), ('sample_extent', numpy.uint64, 1),
                                   ('file_count', numpy.uint64, 1), ('last_timestamp', numpy.double, 1)])
        
        if self.sub_directory_recarray == None:
            # create empty array
            self.sub_directory_recarray = numpy.array([], dtype=self.data_t)
            
        # define glob strings for sub_directories in form YYYY-MM-DDTHH:MM:SS and rf files in form rf@*.*.h5
        self._sub_directory_glob = '[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9]:[0-9][0-9]:[0-9][0-9]'
        self._rf_file_glob = 'rf@[0-9]*.[0-9][0-9][0-9].h5'
        
        
        
    def update(self, start_sample=None, sample_length=None):
        """update will cause this _top_level_dir_metadata object to update itself.
        
        If start_sample and/or sample_length == None, then only get high level metadata.  If 
        start_sample and sample_length given, then only updated detailed metadata as needed.
        
        Inputs:
            start_sample - sample number if units of samples since 1970-01-01.  If None, only refresh high level metadata.
            sample_length - number of samples.  If given, only refresh low-level metadata near sample extent. If None, 
                only refresh high level metadata.
        """
        if start_sample == None or sample_length == None:
            high_level = True
        else:
            high_level = False
            
        if high_level:
            self._high_level_reload()
        else:
            self._low_level_reload(start_sample, sample_length)
           
           
            

            
            
            
    def get_continuous_blocks(self, start_unix_sample, stop_unix_sample):
        """get_continuous_blocks returns a numpy array of dtype u64 and shape (N,2) where the first
        column represents the unix_sample of a continuous block of data, and the second column represents the
        number of samples in that continuous block.  Only samples between (start_unix_sample, stop_unix_sample)
        inclusive will be returned.
        
        Returns IOError if no blocks found
        
        Inputs:
            start_unix_sample, stop_unix_sample - only samples between (start_unix_sample, stop_unix_sample)
                inclusive will be returned.
        """
        # to improve speed, do searchsorted to get first index to look into
        first_index = numpy.searchsorted(self.sub_directory_recarray['unix_start_sample'], 
                                         numpy.array([start_unix_sample]))
        first_index = first_index[0]
        if first_index > 0:
            first_index -= 1
        ret_array = numpy.array([], dtype=numpy.uint64)
        for i in range(first_index, len(self.sub_directory_recarray)):
            this_start_sample = long(self.sub_directory_recarray['unix_start_sample'][i])
            if this_start_sample > stop_unix_sample:
                break
            if stop_unix_sample < this_start_sample:
                continue
            this_extent = long(self.sub_directory_recarray['sample_extent'][i])
            if this_extent == 0:
                raise _MissingMetadata, 'this_extent == 0'
            
            # now check that subdirectories with metadata are still up to date
            base_subdirectory = self.sub_directory_recarray['subdirectory'][i]
            file_count, last_timestamp = self._get_subdirectory_file_info(base_subdirectory)
            self.sub_directory_dict[base_subdirectory].update_if_needed(file_count, last_timestamp)
            
            sub_dir_metadata = self.sub_directory_dict[self.sub_directory_recarray['subdirectory'][i]]
            this_array = sub_dir_metadata.get_continuous_blocks(max(start_unix_sample, this_start_sample),
                                                                min(stop_unix_sample, this_start_sample + this_extent))
            ret_array = self._combine_blocks(ret_array, this_array, self.samples_per_file)
            
        return(ret_array)
    
    
    
    def get_continuous_vector(self, start_unix_sample, stop_unix_sample):
        """get_continuous_vector returns a tuple of (numpy array of data, first unix_sample in returned data)
        Only samples between (start_unix_sample, stop_unix_sample) (excludes stop_unix_sample) will be returned.
        
        Throws IOError if no blocks found
        
        Inputs:
            start_unix_sample, stop_unix_sample - only samples between (start_unix_sample, stop_unix_sample)
                (excludes stop_unix_sample) will be returned.
        """
        # to improve speed, do searchsorted to get first index to look into
        first_index = numpy.searchsorted(self.sub_directory_recarray['unix_start_sample'], 
                                         numpy.array([start_unix_sample]))
        first_index = first_index[0]
        if first_index > 0:
            first_index -= 1
        ret_array = None
        first_unix_sample = None
        for i in range(first_index, len(self.sub_directory_recarray)):
            this_start_sample = long(self.sub_directory_recarray['unix_start_sample'][i])
            if this_start_sample > stop_unix_sample:
                break
            if stop_unix_sample < this_start_sample:
                continue
            this_extent = long(self.sub_directory_recarray['sample_extent'][i])
            if this_extent == 0:
                raise _MissingMetadata, 'this_extent == 0'
            if this_start_sample + this_extent <= start_unix_sample:
                continue
            
            # now check that subdirectories with metadata are still up to date
            base_subdirectory = self.sub_directory_recarray['subdirectory'][i]
            file_count, last_timestamp = self._get_subdirectory_file_info(base_subdirectory)
            self.sub_directory_dict[base_subdirectory].update_if_needed(file_count, last_timestamp)
            
            sub_dir_metadata = self.sub_directory_dict[self.sub_directory_recarray['subdirectory'][i]]
            this_array, unix_sample = sub_dir_metadata.get_continuous_vector(max(start_unix_sample, this_start_sample),
                                                                              min(stop_unix_sample, this_start_sample + this_extent))
            
            ret_array = self._combine_continuous_vectors(ret_array, this_array, first_unix_sample, unix_sample)
            if first_unix_sample == None:
                first_unix_sample = unix_sample
            
        return((ret_array, first_unix_sample))
    
    
    def _high_level_reload(self):
        """_high_level_reload updates only high level metadata
        """
        update_needed = False # will be set to True if any subdirectory updated
        base_subdirectory_list = self._get_subdirectories()
        dt1970 = datetime.datetime(1970,1,1)
        
        # first pass is to remove any subdirectories that have disappeared
        rows_to_delete_arr = []
        for i, subdirectory in enumerate(self.sub_directory_recarray['subdirectory']):
            if subdirectory not in base_subdirectory_list:
                rows_to_delete_arr.append(i)
                del self.sub_directory_dict[subdirectory]
        if len(rows_to_delete_arr) > 0:
            rows_to_delete_arr = numpy.array(rows_to_delete_arr, numpy.int64)
            self.sub_directory_recarray = numpy.delete(self.sub_directory_recarray, rows_to_delete_arr)
            
        
        # next pass creates any new rows in _sub_directory_metadata
        for base_subdirectory in base_subdirectory_list:
            try:
                file_count, last_timestamp = self._get_subdirectory_file_info(base_subdirectory)
            except IOError:
                update_needed = True
                new_sub_dir_meta = _sub_directory_metadata(self.top_level_dir, self.channel_name, 
                                                           self.access_mode, base_subdirectory)
                if self.sub_directory_dict != None:
                    self.sub_directory_dict[base_subdirectory] = new_sub_dir_meta
                else:
                    self.sub_directory_dict = {base_subdirectory: new_sub_dir_meta}
                # extend self.sub_directory_recarray by one
                self.sub_directory_recarray.resize(len(self.sub_directory_recarray) + 1)
                
                # get estimate of first sample without IO
                subDirDT = datetime.datetime.strptime(os.path.basename(base_subdirectory), '%Y-%m-%dT%H:%M:%S')
                total_secs = (subDirDT - dt1970).total_seconds()
                if len(self.metadata_dict.keys()) == 0:
                    # force it to be created
                    new_sub_dir_meta.get_last_sample()
                    self.metadata_dict = new_sub_dir_meta.metadata_dict
                sample_index = total_secs * self.metadata_dict['sample_rate']
                self.sub_directory_recarray[-1] = (base_subdirectory, sample_index, 0, 0, 0) # use default values
            
        self.unix_start_sample = self.sub_directory_dict[base_subdirectory_list[0]].get_first_sample()
        last_sample = long(self.sub_directory_dict[base_subdirectory_list[-1]].get_last_sample())
        self.sample_extent = long(last_sample - self.unix_start_sample)
            
            
            
    def _low_level_reload(self, start_sample, sample_length):
        """_low_level_reload updates detailed metadata for only those subdirectories that need it
        
        Inputs:
            start_sample - sample number if units of samples since 1970-01-01.
            sample_length - number of samples.
        """
        tolerance = 1.001
        
        # see which subdirectories need to be updated
        update_indices = []
        last_start_sample = None
        for i in range(len(self.sub_directory_recarray)):
            this_start_sample = self.sub_directory_recarray['unix_start_sample'][i]
            if last_start_sample != None:
                if last_start_sample > start_sample + long(sample_length * tolerance):
                    break
            this_extent = self.sub_directory_recarray['sample_extent'][i]
            if this_extent > 0:
                if this_start_sample + long(this_extent * tolerance) > start_sample:
                    if i not in update_indices:
                        update_indices.append(i)
                continue
            elif this_start_sample < start_sample:
                last_start_sample = this_start_sample
                continue
            elif last_start_sample != None:
                if i-1 not in update_indices:
                    update_indices.append(i-1)
            last_start_sample = this_start_sample
            
        # see if last item should be appended
        if last_start_sample <= start_sample + sample_length:
            if i not in update_indices:
                update_indices.append(i)
                
        update_indices.sort()
        for i in update_indices:
            sub_directory = self.sub_directory_recarray['subdirectory'][i]
            sub_dir_meta = self.sub_directory_dict[sub_directory]
            sub_dir_meta.update()
            first_unix_sample, sample_extent, file_count, samples_per_file, last_timestamp = \
                sub_dir_meta.get_summary_metadata()
            self.sub_directory_recarray[i] = (sub_directory, first_unix_sample, sample_extent, file_count,
                                                   last_timestamp)
            
            # handle samples_per_file
            if self.samples_per_file == 0:
                self.samples_per_file = long(samples_per_file)
            elif self.samples_per_file != long(samples_per_file):
                raise IOError, 'Samples per file changed from %i to %i with subdirectory %s' % \
                    (self.samples_per_file, samples_per_file, base_subdirectory)
    
    
    
    def _verify_non_overlapping_data(self):
        """_verify_non_overlapping_data raises an error if any overlapping subdirectories found
        """
        for i, record in enumerate(self.sub_directory_recarray):
            if i == 0:
                last_unix_start_sample = record['unix_start_sample']
                last_sample_extent = record['sample_extent']
            else:
                this_unix_start_sample = record['unix_start_sample']
                this_sample_extent = record['sample_extent']
                if last_unix_start_sample + last_sample_extent > this_unix_start_sample:
                    raise IOError, 'Overlapping samples found in subdirectory %s' % (record['subdirectory'])
                last_unix_start_sample = this_unix_start_sample
                last_sample_extent = this_sample_extent
            
                
                
    def _get_subdirectory_file_info(self, subdirectory):
        """_get_subdirectory_file_info returns a tuple ot (num_files, last_timestamp) for a given
        subdirectory using the self.sub_directory_recarray recarray.  Raises IOError if subdirectory
        not found in recarray.
        """
        result = numpy.argwhere(self.sub_directory_recarray['subdirectory'] == subdirectory)
        if len(result) == 0:
            raise IOError, 'subdirectory %s not found' % (subdirectory)
        if len(result) > 1:
            raise ValueError, 'got unexpected result %s' % (str(result))
        return((self.sub_directory_recarray['file_count'][result[0][0]], 
                self.sub_directory_recarray['last_timestamp'][result[0][0]]))
        
        
        
    def _combine_continuous_vectors(self, first_array, second_array, first_start_sample, second_start_sample):
        """_combine_continuous_vectors returns the concatenation of first_array and second_array,  Raises error
        if two vectors are not continuous.
        
        Inputs:
            first_array - first array to combine.  If None, just return second_array
            second_array - second_array to merge at end of first
            first_start_sample - unix_sample of first sample in first_array.  None if first_array == None
            second_start_sample - unix_sample of first sample in second_array
        """
        
        if first_array == None:
            return(second_array)
        
        if len(first_array) != second_start_sample - first_start_sample:
            raise IOError, '_combine_continuous_vectors trying to combine two non-continuous vectors'
        
        return(numpy.concatenate((first_array, second_array)))
                
        
        
    def _get_subdirectories(self):
        """_get_subdirectories returns a sorted list of base subdirectory names
        """
        # for now only local access
        if self.access_mode not in ('local'):
            raise ValueError, 'access_mode %s not yet implemented' % (access_mode)
        subdirectory_list = glob.glob(os.path.join(self.top_level_dir, self.channel_name, self._sub_directory_glob))
        subdirectory_list.sort()
        retList = [] # only return those with files
        for subdirectory in subdirectory_list:
            if len(glob.glob(os.path.join(subdirectory, '*.h5'))) > 0:
                retList.append(subdirectory)
        return(retList)
    
    
    
    def _combine_blocks(self, first_array, second_array, samples_per_file):
        """_combine_blocks combines two numpy array of dtype u64 and shape (N,2) where the first
        column represents the unix_sample of a continuous block of data, and the second column represents the
        number of samples in that continuous block. The first row of the second array may or may not be contiguous
        with the last row of the first array.  If it is contiguous, that row will not be included, and the
        number of samples in that first row will instead be added to the last row of first_array. If not contiguous,
        the two arrays are simply concatenated
        """
        if len(first_array) == 0:
            return(second_array)
        is_contiguous = False
        if first_array[-1][0] + first_array[-1][1] > second_array[0][0]:
            raise IOError, 'overlapping data found in top level directories %i %i' % \
                (first_array[-1][0] + first_array[-1][1], second_array[0][0])
        if first_array[-1][0] + first_array[-1][1] == second_array[0][0]:
            is_contiguous = True
        if is_contiguous:
            first_array[-1][1] += second_array[0][1]
            if len(second_array) == 1:
                return(first_array)
            else:
                return(numpy.concatenate([first_array, second_array[1:]]))
        else:
            return(numpy.concatenate([first_array, second_array]))
    
        
        
    def __cmp__(self, other):
        """__cmp__ compares two _top_level_dir_metadata objects
        """
        # only the same channel can be compared
        if self.channel_name != other.channel_name:
            raise ValueError, 'Cannot compare mismatched channel names %s and %s' % (self.channel_name, other.channel_name)
        
        if self.unix_start_sample != 0 and other.unix_start_sample != 0:
            return(cmp(self.unix_start_sample, other.unix_start_sample))
        
        # use subdirectory names instead
        # for now only local access
        if self.access_mode not in ('local'):
            raise ValueError, 'access_mode %s not yet implemented' % (access_mode)
        
        first_subdirectory_list = glob.glob(os.path.join(self.top_level_dir, self.channel_name, self._sub_directory_glob))
        first_subdirectory_list.sort()
        if len(first_subdirectory_list) == 0:
            raise ValueError, 'Cannot compare top level directory because it has no data' % (self.top_level_dir)
        first_subdirectory = os.path.basename(first_subdirectory_list[0])
        
        second_subdirectory_list = glob.glob(os.path.join(other.top_level_dir, other.channel_name, self._sub_directory_glob))
        second_subdirectory_list.sort()
        if len(second_subdirectory_list) == 0:
            raise ValueError, 'Cannot compare top level directory because it has no data' % (other.top_level_dir)
        second_subdirectory = os.path.basename(second_subdirectory_list[0])
        
        return(cmp(first_subdirectory, second_subdirectory))
    
    
class _sub_directory_metadata:
    """The _sub_directory_metadata is a private class to hold and access metadata about a particular digital_rf channel in
    a particular subdirectory.
    """
    
    def __init__(self, top_level_dir, channel_name, access_mode, subdirectory):
        """__init__ creates a new _sub_directory_metadata object
        
        Inputs:
            top_level_dir - full path the top level directory that contains the parent channel_name
            channel_name - the channel_name subdirectory name
            access_mode - string giving access mode (eg, 'local', 'file', or 'http')
            subdirectory - subdirectory name in form YYYY-MM-DDTHH:MM:SS
            
        Affects:
            Sets self.metadata to None.  When update called, self.metadata will be set either to a numpy.recarray
            with columns:  
                1. unix_sample_index - number of samples since 1970-01-01 to the start of a contiguous data block (uint64_t)
                2. file_index - where in the file this contiguous block of data begins
                3. rf_basename (25 char string)
            or to a string giving the full path to an Hdf5 file with that recarray on disk on an Hdf5 file under the dataset
            name /sub_directory_metadata
            
            Also sets self.continuous_metadata to None.  When update called, self.continuous_metadata will be set 
            either to a numpy.recarray about block of contiguous data with columns:  
                1. unix_sample_index - number of samples since 1970-01-01 to the start of a contiguous data block (uint64_t)
                2. sample_extent - number of continuous samples
            or to a string giving the full path to an Hdf5 file with that recarray on disk on an Hdf5 file under the dataset
            name /sub_directory_continuous_metadata.
            
            Also sets self.samples_per_file and self.file_count and self.last_timestamp to None.  Will be set with 
            first call to update
          
        """
        self.top_level_dir = top_level_dir
        self.channel_name = channel_name
        self.access_mode = access_mode
        self.subdirectory = subdirectory
        self.metadata = None
        self.cont_metadata = None
        self.samples_per_file = None
        self.file_count = None
        self.last_timestamp = None # timestamp of last file in UTC
        self.metadata_dict = {} # to be populated by rf file metadata
        
        self._rf_file_glob = 'rf@[0-9]*.[0-9][0-9][0-9].h5'
        
        # data type of sub_directory_array
        self.data_t = numpy.dtype([('unix_sample_index', numpy.uint64, 1), ('file_index', numpy.uint64, 1), 
                                   ('rf_basename', numpy.str_, 25)])
        self.cont_data_t = numpy.dtype([('unix_sample_index', numpy.uint64, 1), ('sample_extent', numpy.uint64, 1)])
        
        # see if there already exists a Hdf5 metadata file
        if self.access_mode == 'local':
            metadata_file = os.path.join(top_level_dir, channel_name, subdirectory, 'sub_directory_metadata.h5')
            if os.access(metadata_file, os.R_OK):
                self.metadata = metadata_file
                
        # see if there already exists a Hdf5 continuous metadata file
        if self.access_mode == 'local':
            cont_metadata_file = os.path.join(top_level_dir, channel_name, subdirectory, 'sub_directory_continuous_metadata.h5')
            if os.access(cont_metadata_file, os.R_OK):
                self.cont_metadata = cont_metadata_file
                
        # if not set to a string, set to an empty recarray
        if self.metadata == None:
            self.metadata = numpy.array([], dtype=self.data_t)
        if self.cont_metadata == None:
            self.cont_metadata = numpy.array([], dtype=self.cont_data_t)
            
            
    def get_summary_metadata(self):
        """get_summary_metadata returns a tuple of (first_unix_sample, sample_extent, file_count, samples_per_file,
        last_timestamp) for this _sub_directory_metadata object.
        
        Raises IOError if no self.metadata
        """
        if type(self.metadata) in (types.StringType,):
            self.read_into_memory()
            
        if len(self.metadata) == 0:
            raise IOError, 'Must call update before calling get_summary_metadata, or subdirectory %s empty' % (self.subdirectory)
        
        first_unix_sample = long(self.metadata['unix_sample_index'][0])
        last_unix_sample = long(self.metadata['unix_sample_index'][-1]) + ((self.samples_per_file - long(self.metadata['file_index'][-1])) - 1)
        return((first_unix_sample, 1+(last_unix_sample-first_unix_sample), self.file_count, self.samples_per_file,
                self.last_timestamp))
        
        
    def update_if_needed(self, file_count, last_timestamp):
        """update_if_needed calls update only if input file_count or last_timestamp indicate an update
        is needed. Returns True if update actually called, False otherwise
        
        Inputs:
            file_count - number of file in subdirectory when last checked
            last_timestamp - UTC timestamp of last file in subdirectory when last checked
        """
        if type(self.metadata) in (types.StringType,):
            self.read_into_memory()
            
        # for now only local access
        if self.access_mode not in ('local'):
            raise ValueError, 'access_mode %s not yet implemented' % (access_mode)
        
        rf_file_list = glob.glob(os.path.join(self.top_level_dir, self.channel_name, self.subdirectory, 
                                              self._rf_file_glob))
        if len(rf_file_list) == 0:
            raise IOError, 'subdirectory %s empty' % (self.subdirectory)
        
        rf_file_list.sort()
        if len(rf_file_list) != file_count:
            self.update()
            return(True)
        elif abs(self._get_utc_timestamp(rf_file_list[-1]) - last_timestamp) > 2.0:
            # leave margin for error
            self.update()
            return(True)
        return(False)
                
                
                
    def update(self):
        """update updates self.metadata.  If it was a file name, it reads that data into memory, and then updates it
        """
        if type(self.metadata) in (types.StringType,):
            self.read_into_memory()
            
        # for now only local access
        if self.access_mode not in ('local'):
            raise ValueError, 'access_mode %s not yet implemented' % (access_mode)
        
        rf_file_list = glob.glob(os.path.join(self.top_level_dir, self.channel_name, self.subdirectory, 
                                              self._rf_file_glob))
        rf_file_list.sort()
        rf_file_basename_list = [os.path.basename(rf_file) for rf_file in rf_file_list]
        
        # first check to see if we can update things quickly if the data is continuous
        if self._update_continuous_data(rf_file_basename_list, rf_file_list):
            return
        
        unique_rf_basenames = numpy.unique(self.metadata['rf_basename'])
        
        # first step is to delete all lines where the rf file has been deleted
        rows_to_delete_arr = numpy.array([], dtype=numpy.int64)
        for rf_basename in unique_rf_basenames:
            if rf_basename not in rf_file_basename_list:
                # get a list of all rows with that file
                result = numpy.argwhere(self.metadata['rf_basename'] == rf_basename)
                result = result.flatten()
                rows_to_delete_arr = numpy.concatenate((rows_to_delete_arr, result))
        if len(rows_to_delete_arr) > 0:
            self.metadata = numpy.delete(self.metadata, rows_to_delete_arr)
            unique_rf_basenames = numpy.unique(self.metadata['rf_basename'])
            
        # the next step is to add rows from each file where it does not yet exist in self.metadata
        if len(self.metadata['rf_basename']) > 0:
            first_file_index = rf_file_basename_list.index(self.metadata['rf_basename'][-1]) + 1
        else:
            first_file_index = 0
        # we are only looping over files not already in self.metadata, and all data will be appended
        for i, rf_file_basename in enumerate(rf_file_basename_list[first_file_index:]):
            # verify the last file is not still being written
            if rf_file_basename == rf_file_basename_list[-1]:
                if self._file_is_open(rf_file_list[-1]):
                    self.file_count = len(rf_file_basename_list) - 1 # last file not counted
                    if self.file_count > 0:
                        self.last_timestamp = self._get_utc_timestamp(rf_file_list[-2])
                    continue
                else:
                    self.file_count = len(rf_file_basename_list)
                    self.last_timestamp = self._get_utc_timestamp(rf_file_list[-1])
            added_rows = self._get_new_rows(rf_file_basename)
            if added_rows == None:
                continue
            if len(self.metadata_dict.keys()) == 0:
                self.metadata_dict = self._get_rf_metadata(rf_file_basename)
            self.metadata = numpy.concatenate((self.metadata, added_rows))
            
        self._update_cont_metadata()
            
            
    def get_continuous_blocks(self, start_unix_sample, stop_unix_sample):
        """get_continuous_blocks returns a numpy array of dtype u64 and shape (N,2) where the first
        column represents the unix_sample of a continuous block of data, and the second column represents the
        number of samples in that continuous block.  Only samples between (start_unix_sample, stop_unix_sample)
        inclusive will be returned.
        
        Returns IOError if no blocks found
        
        Inputs:
            start_unix_sample, stop_unix_sample - only samples between (start_unix_sample, stop_unix_sample)
                inclusive will be returned.
        """
        # to improve speed, do searchsorted to get first index to look into
        first_index = numpy.searchsorted(self.cont_metadata['unix_sample_index'], 
                                         numpy.array([start_unix_sample]))
        first_index = first_index[0]
        if first_index > 0:
            first_index -= 1
            
        # for now, deal with two edges later
        bool_arr = self.cont_metadata['unix_sample_index'] >= start_unix_sample
        bool_arr1 = self.cont_metadata['unix_sample_index'] <= stop_unix_sample
        bool_arr = numpy.logical_and(bool_arr, bool_arr1)
        metadata_slice = self.cont_metadata[bool_arr]
        
        # handle front edge
        ones = numpy.ones((len(bool_arr),))
        zeros = numpy.zeros((len(bool_arr),))
        indices = numpy.nonzero(numpy.where(bool_arr, ones, zeros))
        if len(indices[0]):
            first_index = indices[0][0]
            if first_index != 0:
                # in this case, we may need to add one line before, but modify it
                previous_sample_index = self.cont_metadata[first_index-1]['unix_sample_index']
                previous_sample_extent = self.cont_metadata[first_index-1]['sample_extent']
                if previous_sample_index + previous_sample_extent >= start_unix_sample:
                    metadata_slice = numpy.concatenate((self.cont_metadata[first_index-1:first_index],
                                                       metadata_slice))
                    metadata_slice[0] = (start_unix_sample, previous_sample_extent - \
                        (start_unix_sample - previous_sample_index))
            # else there is no need to fix front edge
        else:
            previous_sample_index = self.cont_metadata[0]['unix_sample_index']
            previous_sample_extent = self.cont_metadata[0]['sample_extent']
            if previous_sample_index + previous_sample_extent >= start_unix_sample:
                metadata_slice = numpy.zeros((1,), dtype=self.cont_data_t)
                metadata_slice[0] = (start_unix_sample, previous_sample_extent - \
                    (start_unix_sample - previous_sample_index))
                
        # fix end if need
        last_sample_index = metadata_slice[-1]['unix_sample_index']
        last_sample_extent = metadata_slice[-1]['sample_extent']
        real_last_sample_extent = 1 + (stop_unix_sample - last_sample_index)
        if real_last_sample_extent < last_sample_extent:
            metadata_slice[-1]['sample_extent'] = real_last_sample_extent
        
        ret_arr = numpy.zeros((len(metadata_slice), 2), dtype=numpy.uint64)
        ret_arr[:,0] = metadata_slice['unix_sample_index']
        ret_arr[:,1] = metadata_slice['sample_extent']
        
        return(ret_arr)
    
    
    def get_continuous_vector(self, start_unix_sample, stop_unix_sample):
        """get_continuous_vector returns a tuple of (numpy array of data, first unix_sample in returned data)
        Only samples between (start_unix_sample, stop_unix_sample) (excludes stop_unix_sample) will be returned.
        
        Returns IOError if no blocks found
        
        Inputs:
            start_unix_sample, stop_unix_sample - only samples between (start_unix_sample, stop_unix_sample)
                (excludes stop_unix_sample) will be returned.
        """
        # first verify no gaps in this range
        self._verify_no_gaps(start_unix_sample, stop_unix_sample)
        
        # to improve speed, do searchsorted to get first index to look into
        first_index = numpy.searchsorted(self.metadata['unix_sample_index'], 
                                         numpy.array([start_unix_sample]))
        first_index = long(first_index[0])
        if first_index == len(self.metadata):
            first_index -= 1
        elif self.metadata['unix_sample_index'][first_index] > start_unix_sample:
            first_index -= 1
        ret_array = None
        first_unix_sample = None
        this_hdf5_file = None
        samples_read = 0
        samples_to_read = stop_unix_sample - start_unix_sample
        for i in range(first_index, len(self.metadata)):
            if self.metadata['unix_sample_index'][i] >= stop_unix_sample:
                raise IOError, 'Did not get expected read - debug'
            
            this_hdf5_file = self.metadata['rf_basename'][i]
            full_hdf5_file = os.path.join(self.top_level_dir, self.channel_name, self.subdirectory, this_hdf5_file)
            
            # get max possible length of this read as block_len
            if i == len(self.metadata) - 1:
                # last index
                block_len = self.samples_per_file - long(self.metadata['file_index'][i])
            elif self.metadata['rf_basename'][i+1] == this_hdf5_file:
                block_len = long(self.metadata['file_index'][i+1]) - long(self.metadata['file_index'][i])
            else:
                block_len = self.samples_per_file - long(self.metadata['file_index'][i])
                
            # next get file start index
            if i == first_index:
                offset = long(start_unix_sample) - long(self.metadata['unix_sample_index'][first_index])
                block_len -= offset # we will not get the full read predicted above
                start_file_index = long(self.metadata['file_index'][i]) + offset
            else:
                start_file_index = long(self.metadata['file_index'][i])
                
            # next get read len
            if samples_read + block_len > samples_to_read:
                read_len = samples_to_read - samples_read
            else:
                read_len = block_len
                
            # finally - read it!!!
            f = h5py.File(full_hdf5_file, 'r')
            rf_data = f['/rf_data'][start_file_index:start_file_index + read_len]
            f.close()
            
            if ret_array == None:
                ret_array = rf_data
            else:
                ret_array = numpy.concatenate((ret_array, rf_data))
            samples_read += read_len
            if samples_read == samples_to_read:
                break
                
        return((ret_array, start_unix_sample))
        
        
    def read_into_memory(self):
        """read_into_memory reads data into recarray form into self.metadata
        """
        # temp only
        pass
    
    
    
    def write_to_memory(self):
        """write_to_memory writes the numpy.recarray in self.metadata to an Hdf5 file, releases recarray from memory
        """
        # temp only
        pass
    
    
    def get_first_sample(self):
        """get_first_sample returns the first sample index in this subdirectory.  May be exact (if self.metadata
        not == None) or an estimate based one subdirectory naming convention.
        """
        if len(self.metadata) > 0:
            return(self.metadata['unix_sample_index'][0])
        
        rf_file_list = glob.glob(os.path.join(self.top_level_dir, self.channel_name, self.subdirectory, 
                                              self._rf_file_glob))
        rf_file_list.sort()
        
        new_rows = self._get_new_rows(os.path.basename(rf_file_list[0]))
        return(new_rows['unix_sample_index'][0])
    
    
    def get_last_sample(self):
        """get_last_sample returns the last sample in a subdirectory.  
        """
        if len(self.metadata) > 0 and len(self.metadata_dict.keys()):
            return(self.metadata['unix_sample_index'][-1] + self.metadata_dict['samples_per_file'] - self.metadata['file_index'][-1])
        
        rf_file_list = glob.glob(os.path.join(self.top_level_dir, self.channel_name, self.subdirectory, 
                                              self._rf_file_glob))
        rf_file_list.sort()
        
        if len(self.metadata_dict.keys()) == 0:
            self.metadata_dict = self._get_rf_metadata(os.path.basename(rf_file_list[0]))
            
        if not self._file_is_open(rf_file_list[-1]):
            index = -1
        else:
            index = -2
            
        new_rows = self._get_new_rows(os.path.basename(rf_file_list[index]))
        return(new_rows['unix_sample_index'][-1] + self.metadata_dict['samples_per_file'] - new_rows['file_index'][-1])
        
        
    
    
    
    def _update_continuous_data(self, rf_file_basename_list, rf_file_list):
        """_update_continuous_data updates all metadata if data in subdirectory is continuous, then return True.  Does nothing
        and returns False if not continuous data.  
        
        Determines if continuous by looking only at first and last file in rf_file_basename_list.
        
        Inputs:
            rf_file_basename_list - sorted list of basenames in subdirectory
            rf_file_list - sorted list of full names  in subdirectory
        """
        first_index = None
        last_index = None
        
        # get info from first file
        for i, rf_file_basename in enumerate(rf_file_basename_list):
            if first_index == None:
                first_row = self._get_new_rows(rf_file_basename)
                if first_row != None:
                    if len(first_row) > 1:
                        return(False)
                    first_index = i
                    first_sample = first_row['unix_sample_index']
                    break
                else:
                    continue
                
        if first_index == None:
            return(False)
        
        # get info from last file
        if not self._file_is_open(rf_file_list[-1]):
            last_row = self._get_new_rows(rf_file_basename_list[-1])
            last_index = len(rf_file_basename_list)
            self.last_timestamp = self._get_utc_timestamp(rf_file_list[-1])
        else:
            last_row = self._get_new_rows(rf_file_basename_list[-2])
            last_index = len(rf_file_basename_list) - 1
            self.last_timestamp = self._get_utc_timestamp(rf_file_list[-2])
        if last_row == None:
            return(False)
        if len(last_row) > 1:
            return(False)
        last_sample = last_row['unix_sample_index']
        self.file_count = last_index - first_index
        if (self.file_count - 1) * self.samples_per_file < last_sample - first_sample:
            # data gaps detected
            return(False)
        
        # create all metadata
        self.metadata = numpy.recarray((self.file_count,), dtype = self.data_t)
        sample_data = numpy.arange(0,self.file_count*self.samples_per_file, self.samples_per_file,dtype=numpy.int64)
        sample_data += first_sample
        self.metadata['unix_sample_index'] = sample_data
        self.metadata['file_index'][:] = 0
        self.metadata['rf_basename'] = rf_file_basename_list[first_index:last_index]
         
        self._update_cont_metadata()
        return(True)
    
    
    
    
    def _verify_no_gaps(self, start_unix_sample, stop_unix_sample):
        """_verify_no_gaps raises an IOError if there is a gap between start_unix_sample, stop_unix_sample
        """
         # to improve speed, do searchsorted to get first index to look into
        first_index = numpy.searchsorted(self.cont_metadata['unix_sample_index'], 
                                         numpy.array([start_unix_sample]))
        first_index = first_index[0]
        if first_index == len(self.cont_metadata):
            first_index -= 1
        elif self.cont_metadata['unix_sample_index'][first_index] > start_unix_sample:
            first_index -= 1
        offset = start_unix_sample - long(self.cont_metadata['unix_sample_index'][first_index])
        if self.cont_metadata['sample_extent'][first_index] - offset < stop_unix_sample - start_unix_sample:
            raise IOError, 'gap found between samples %i and %i' % (start_unix_sample, stop_unix_sample)
        
    
    
    def _get_new_rows(self, rf_file_basename):
        """_get_new_rows is a private method that returns all needed rows for self.metadata in the correct recarray
        format for rf_file_basename, or None if that file has disappeared
        
        Inputs:
            rf_file_basename - rf file to examine

        Throws IOError if global indices overlap with previous metadata
        """
        # read data from /rf_data_index
        fullname = os.path.join(self.top_level_dir, self.channel_name, self.subdirectory, rf_file_basename)
        try:
            f = h5py.File(fullname, 'r')
        except IOError:
            # presumably file deleted
            return(None)
        rf_data_index = f['/rf_data_index']
        samples_per_file = f['rf_data'].attrs['samples_per_file'][0]
        if self.samples_per_file == None:
            self.samples_per_file = int(samples_per_file)
        elif self.samples_per_file != int(samples_per_file):
            raise IOError, 'Illegal change in samples_per_file from %i to %i in file %s' % (self.samples_per_file, int(samples_per_file),
                                                                                            fullname)
            
        # create recarray
        new_rows = numpy.zeros((len(rf_data_index),),dtype=self.data_t)
        new_rows['unix_sample_index'] = rf_data_index[:,0]
        new_rows['file_index'] = rf_data_index[:,1]
        new_rows['rf_basename'] = rf_file_basename
        
        f.close()
        
        return(new_rows)
    
    
    def _get_rf_metadata(self, rf_file_basename):
        """_get_rf_metadata is a private method that returns a dictionary of all metadata stored in each rf file,
        or empty dict if that file has disappeared
        
        Inputs:
            rf_file_basename - rf file to examine

        Returns dictionary with string keys:
            sample_rate
            samples_per_file
            uuid_str
        """
        ret_dict = {}
        fullname = os.path.join(self.top_level_dir, self.channel_name, self.subdirectory, rf_file_basename)
        try:
            f = h5py.File(fullname, 'r')
        except IOError:
            return({})
        dataset = f['/rf_data']
        for attr in dataset.attrs:
            ret_dict[str(attr)] = dataset.attrs[attr]
        
        f.close()
        return(ret_dict)
        
    
    
    def _combine_blocks(self, first_array, second_array, samples_per_file):
        """_combine_blocks combines two numpy array of dtype u64 and shape (N,2) where the first
        column represents the unix_sample of a continuous block of data, and the second column represents the
        number of samples in that continuous block. The first row of the second array may or may not be contiguous
        with the last row of the first array.  If it is contiguous, that row will not be included, and the
        number of samples in that first row will instead be added to the last row of first_array. If not contiguous,
        the two arrays are simply concatenated
        """
        if len(first_array) == 0:
            return(second_array)
        is_contiguous = False
        if first_array[-1][0] + first_array[-1][1] > second_array[0][0]:
            raise IOError, 'overlapping data found in top level directories %i, %i' % \
                (first_array[-1][0] + first_array[-1][1], second_array[0][0])
        if first_array[-1][0] + first_array[-1][1] == second_array[0][0]:
            is_contiguous = True
        if is_contiguous:
            first_array[-1][1] += second_array[0][1]
            if len(second_array) == 1:
                return(first_array)
            else:
                return(numpy.concatenate([first_array, second_array[1:]]))
        else:
            return(numpy.concatenate([first_array, second_array]))
        
        
    def _update_cont_metadata(self):
        """_update_cont_metadata completely rebuilds self.cont_metadata
        """
        cont_meta = []
        # handle empty dir case
        if len(self.metadata) == 0:
            self.cont_metadata = numpy.zeros((len(cont_meta),), dtype=self.cont_data_t)
            return
        for i in range(len(self.metadata)):
            if i == 0:
                cont_meta.append([self.metadata['unix_sample_index'][0], 0])
                last_sample = self.metadata['unix_sample_index'][0]
                last_index = 0
                continue
            this_sample = self.metadata['unix_sample_index'][i]
            this_index = self.metadata['file_index'][i]
            if this_index == 0:
                num_samples = self.samples_per_file - last_index
            else:
                num_samples = this_index - last_index
                if num_samples < 1:
                    raise ValueError, 'bug in self.metadata'
            cont_meta[-1][1] += num_samples
            if this_sample - last_sample == num_samples:
                if this_index != 0:
                    raise ValueError, 'bug 2 in self.metadata'
            else:
                cont_meta.append([this_sample, 0])
            last_sample = this_sample
            last_index = this_index
                
        # handle end of last file
        edge_samples = self.samples_per_file - last_index

        cont_meta[-1][1] += edge_samples
        
        cont_meta = numpy.array(cont_meta)
        
        # create self.cont_metadata
        self.cont_metadata = numpy.zeros((len(cont_meta),), dtype=self.cont_data_t)
        self.cont_metadata['unix_sample_index'] = cont_meta[:,0]
        self.cont_metadata['sample_extent'] = cont_meta[:,1]
        
    
    
    def _file_is_open(self, rf_file):
        """_file_is_open returns True if rf_file might be open (or corrupt), False otherwise
        """
        if time.time() - os.path.getmtime(rf_file) < 3:
            return(True)
        else:
            try:
                f = h5py.File(rf_file)
                f.close()
                return(False)
            except:
                return(True)
        
    def _get_utc_timestamp(self, fullfile):
        """_get_utc_timestamp returns the last modification timestamp of fullfile in UTC
        """
        # for now only local access
        if self.access_mode not in ('local'):
            raise ValueError, 'access_mode %s not yet implemented' % (access_mode)
        
        return(os.path.getmtime(fullfile) - time.timezone)
    
    
class _MissingMetadata(Exception):
    """_MissingMetadata is a Exception that will be raised when metadata needs to be updated
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
        
    
    