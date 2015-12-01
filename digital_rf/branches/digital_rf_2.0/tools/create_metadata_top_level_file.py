"""create_metadata_top_level_file.py is a temporary piece of code to create top level
metadata files (metadata.h5) for metadata channels.

$Id: create_metadata_top_level_file.py 827 2015-11-06 15:01:44Z brideout $
"""
usage = """python create_metadata_top_level_file.py <destDir> <subdirectory_secs> <file_secs> <file_name> <samples_per_second> <has_idx>
where <destDir> is where metadata.h5 will be written, subdirectory_secs is the integer number of seconds per 
top level directory in form YYYY-MM-DDTHH-MM-SS, file_secs is the integer number of seconds per file
in the form <file_name>@<second>.h5, sample rate is samples per second, and has_idx is 1 if file has an idx column indicating sample
metadata recorded at, 0 if only a single measurment without index.
"""

# standard python imports
import os, os.path, sys

# third party imports
import h5py
import numpy

def create_metadata_top_level_file(destDir, subdirectory_secs, file_secs, file_name, samples_per_second, has_idx):
    """create_metadata_top_level_file creates metadata.h5 in destDir
    
    Inputs:
        destDir - where metadata.h5 will be written
        subdirectory_secs - the integer number of seconds per top level directory in form YYYY-MM-DDTHH-MM-SS
        file_secs - the integer number of seconds per file
        file_name - first part of file name in form <file_name>@<second>.h5
        samples_per_second - samples per second
        has_idx - 1 if file has an idx column indicating sample metadata recorded at, 
            0 if only a single measurment without index
    """
    if has_idx not in (0, 1):
        raise ValueError, 'has_idx must be 0 or 1, not <%s>' % (str(has_idx))
    
    # hack - too complex to pass in from command line!!!
    if has_idx:
        recarray = numpy.array([('idx', 'int', False), 
                                ('modeid', 'int', False),
                                ('sample_rate', 'float', True),
                                ('sweepid', 'int', False),
                                ('sweepnum', 'int', False)], 
                               dtype=[('column', '|S256'), ('type', '|S5'), ('isSingleValued', bool)])
    else:
        print('no idx')
        recarray = numpy.array([('cal_temp', 'float', True), 
                                ('cal_temp_misa', 'float', True), 
                                ('cal_temp_zenith', 'float', True), 
                                ('misa_off', 'float', True), 
                                ('misa_off_100_khz_median', 'float', True), 
                                ('misa_on', 'float', True), 
                                ('misa_on_100_khz_median', 'float', True), 
                                ('misa_sr', 'float', True), 
                                ('misa_tx_power', 'float', True), 
                                ('n_off_samples', 'float', True), 
                                ('n_on_samples', 'float', True), 
                                ('res_misa_off_dc', 'complex', True), 
                                ('res_misa_on_cd', 'complex', True), 
                                ('res_zenith_off_dc', 'complex', True), 
                                ('res_zenith_on_cd', 'complex', True), 
                                ('t0', 'float', True), 
                                ('t1', 'float', True), 
                                ('temp_misa', 'float', True), 
                                ('temp_misa_median', 'float', True), 
                                ('temp_zenith', 'float', True), 
                                ('temp_zenith_median', 'float', True), 
                                ('zenith_off', 'float', True), 
                                ('zenith_off_100_khz_median', 'float', True), 
                                ('zenith_on', 'float', True), 
                                ('zenith_on_100_khz_median', 'float', True), 
                                ('zenith_sr', 'float', True), 
                                ('zenith_tx_power', 'float', True)], 
                               dtype=[('column', '|S256'), ('type', '|S8'), ('isSingleValued', bool)])
    with h5py.File(os.path.join(destDir, 'metadata.h5')) as f:
        f.attrs['subdirectory_cadence_seconds'] = int(subdirectory_secs)
        f.attrs['file_cadence_seconds'] = int(file_secs)
        f.attrs['samples_per_second'] = samples_per_second
        f.attrs['file_name'] = str(file_name)
        f.attrs['has_idx'] = bool(has_idx)
        f.create_dataset("fields", data=recarray)
        
    
    
if len(sys.argv) != 7:
    print(usage)
    sys.exit(-1)
    
destDir = sys.argv[1]
if not os.access(destDir, os.W_OK):
    raise IOError, 'destDir %s not writable' % (destDir)
subdirectory_secs = int(sys.argv[2])
if subdirectory_secs < 2:
    raise ValueError, 'subdirectory_secs must be greater than 1, not %i' % (subdirectory_secs)
file_secs = int(sys.argv[3])
if file_secs >= subdirectory_secs or file_secs < 1:
    raise ValueError, 'file_secs must be less than subdirectory_secs but greater than zero, not %i' % (file_secs)
file_name = sys.argv[4]
samples_per_second = float(sys.argv[5])
if samples_per_second <= 0.0:
    raise ValueError, 'samples_per_second must be positive, not %f' % (samples_per_second)
has_idx = int(sys.argv[6])
if has_idx not in (0, 1):
        raise ValueError, 'has_idx must be 0 or 1, not <%s>' % (str(has_idx))
create_metadata_top_level_file(destDir, subdirectory_secs, file_secs, file_name, samples_per_second, has_idx)
print('complete')
        
        