/*******************************************************
 * digital_rf.h is the header file for the rf_write_hdf5 C library
 *
 * The rf_write_hdf5 library supports the writing of rf data
 * into the Hdf5 file as specified in the accompanying documentation.  This
 * format is design to support RF archiving with excellent random access capability.
 *
 * Note that there are four different indices used in this API, and the names of these indices
 * are named with the prefixes:
 * 	global_ - this is the overall index, where 0 represents the first sample recorded, and
 * 		the index always refers to the number of sample periods since that first sample. Note now that
 * 		what is stored in the Hdf5 file is global_index + index of first sample, which is the number of samples
 * 		between midnight UT 1970-01-01 and the start of the experiment at the given sample_rate.  This mean this
 * 		index is an absolute UTC time (leap seconds are ignored, unless they occur during data taking).
 * 	data_ - this is the index into the block of data as passed in by the user, and zero is always
 * 		the first sample in that data array passed in.
 * 	dataset_ - this index always refers to a position in the Hdf5 rf_data dataset in a particular Hdf5 file
 * 	block_ - this index always refers to a position in the Hdf5 rf_data_index dataset that stores indices
 * 		into /rf_data
 *
 * 	Written by Bill Rideout (brideout@haystack.mit.edu), in collaboration with Juha Vierinan (x@haystack.mit.edu)
 * 	and Frank Lind (flind@haystack.mit.edu)
 *
 * $Id: digital_rf.h 782 2015-07-07 14:36:49Z brideout $
 */

#ifndef _RF_WRITE_HDF5_
#define _RF_WRITE_HDF5_

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <inttypes.h>

#include "hdf5.h"
# include "H5Tpublic.h"

/* constants */
#define DIGITAL_RF_VERSION "1.0"

/* string sizes */
#define SMALL_HDF5_STR 265
#define MED_HDF5_STR 512
#define BIG_HDF5_STR 1024

/* chunk size for rf_data_index */
#define CHUNK_SIZE_RF_DATA_INDEX 100

#define DIGITAL_RF_EPOCH "1970-01-01T00:00:00Z"
#define DIGITAL_RF_TIME_DESCRIPTION "All times in this format are in number of samples since the epoch in the epoch attribute.  The first sample time will be sample_rate * UTC time at first sample.  Attribute init_utc_timestamp records this init UTC time so that a conversion to any other time is possible given the number of leapseconds difference at init_utc_timestamp.  Leapseconds that occur during data recording are included in the data."


typedef struct digital_rf_write_object {

    /* this structure encapsulates all information needed to write to a series of Hdf5 files in a directory */
	char *     directory;       		/* Channel directory name where all data is stored - will always end with a "/" */
	char *     sub_directory;               /* Present sub-directory in form YYYY-MM-DDTHH:MM:SS - will always end with a "/" */
	int        is_complex;              /* 1 if complex (IQ) data, 0 if single valued */
	int        num_subchannels;         /* number of subchannels in the data stream.  Must be at least 1. */
	int        rank;            		/* 2 if complex (IQ) data or num_subchannels > 1, 1 otherwise */
	char *     uuid_str;        		/* UUID in str form */
	uint64_t   samples_per_file;		/* number of IQ or single valued samples in any Hdf5 file */
	uint64_t   files_per_directory; 	/* number of Hdf5 files before creating a new directory. If 0, break on hour */
	uint64_t   directory_last_hour;     /* the hour of the last directory created. Used only if files_per_directory == 0  */
	uint64_t   global_start_sample;     /* time of first sample in number of samples since UT midnight 1970-01-01 */
	double     sample_rate;             /* sample rate in Hz */
	int        needs_chunking;  		/* 1 if /rf_data needs chunking (either compression or checksums used) */
	int        chunk_size;      		/* chunk size used - left at 0 if no chunking. 0 also indicates chunking property not yet set */
	hid_t      dtype_id;        		/* individual field data type as defined by hdf5.h */
	hid_t      complex_dtype_id;        /* complex compound data type if is_complex, with fields r and i */
	uint64_t   global_index;    		/* index into the next sample that could be written (global) */
	int        present_seq;     		/* The present Hdf5 file sequence. Init value is -1 */
	uint64_t   dataset_index;   		/* the next available index in open Hdf5 to write to */
	uint64_t   dataset_avail;   		/* the number of samples in the dataset available for writing to */
	uint64_t   block_index;     		/* the next available row in the open Hdf5 file/rf_data_index dataset to write to */
	hid_t      dataset;         		/* Dataset presently opened            */
	hid_t      dataspace;       		/* Dataspace used (rf_data)            */
	hid_t      filespace;       		/* filespace object used               */
	hid_t      memspace;        		/* memspace object used                */
	hid_t      hdf5_file;       		/* Hdf5 file presently opened          */
	hid_t      dataset_prop;    		/* Hdf5 dataset property               */
	hid_t      index_dataset;   		/* Hdf5 rf_data_index dataset          */
	hid_t      index_prop;      		/* Hdf5 rf_data_index property         */
	int        next_index_avail;		/* the next available row in /rf_data_index */
	int        marching_dots;           /* non-zero if marching dots desired when writing, 0 if not */
	uint64_t   init_utc_timestamp;      /* unix time when channel init called - stored as attribute in each file */

} Digital_rf_write_object;

/* Public method declarations */


int digital_rf_write_blocks_hdf5(Digital_rf_write_object *hdf5_data_object, uint64_t * global_index_arr, uint64_t * data_index_arr,
		                         uint64_t index_len, void * vector, uint64_t vector_length);


#ifdef __cplusplus
extern "C" int digital_rf_get_unix_time(uint64_t, double, int*, int*, int*,
		                     int*, int*, int*, uint64_t*); 
extern "C" Digital_rf_write_object * digital_rf_create_write_hdf5(char*, hid_t, uint64_t,
					               uint64_t, uint64_t,
				                       double, char *,
		                                       int, int, int,
		                                       int, int);
extern "C" int digital_rf_write_hdf5(Digital_rf_write_object*, uint64_t, void*,uint64_t);
extern "C" int digital_rf_close_write_hdf5(Digital_rf_write_object*);

#else
int digital_rf_get_unix_time(uint64_t global_sample, double sample_rate, int * year, int * month, int *day,
		                     int * hour, int * minute, int * second, uint64_t * picosecond);
Digital_rf_write_object * digital_rf_create_write_hdf5(char * directory, hid_t dtype_id, uint64_t samples_per_file,
					               uint64_t files_per_directory, uint64_t global_start_sample,
				                       double sample_rate, char * uuid_str,
		                                       int compression_level, int checksum, int is_complex,
		                                       int num_subchannels, int marching_dots);
int digital_rf_write_hdf5(Digital_rf_write_object *hdf5_data_object, uint64_t global_leading_edge_index, void * vector,
						  uint64_t vector_length);
int digital_rf_close_write_hdf5(Digital_rf_write_object *hdf5_data_object);
#endif

/* Private method declarations */
int digital_rf_free_hdf5_data_object(Digital_rf_write_object *hdf5_data_object);
int digital_rf_check_hdf5_directory(char * directory);
uint64_t digital_rf_write_samples_to_file(Digital_rf_write_object *hdf5_data_object, uint64_t samples_written, uint64_t * global_index_arr,
		uint64_t * data_index_arr, uint64_t index_len, void * vector, uint64_t vector_length);
int digital_rf_create_hdf5_file(Digital_rf_write_object *hdf5_data_object, uint64_t next_global_sample);
int digital_rf_create_new_directory(Digital_rf_write_object *hdf5_data_object, uint64_t next_global_sample);
int digital_rf_set_fill_value(Digital_rf_write_object *hdf5_data_object);
void digital_rf_write_metadata(Digital_rf_write_object *hdf5_data_object);
uint64_t * digital_rf_create_rf_data_index(Digital_rf_write_object *hdf5_data_object, uint64_t samples_written, uint64_t * global_index_arr,
			uint64_t * data_index_arr, uint64_t index_len, int * rows_to_write);
int digital_rf_write_rf_data_index(Digital_rf_write_object * hdf5_data_object, uint64_t * rf_data_index_arr, int block_index_len);
uint64_t digital_rf_get_global_sample(uint64_t samples_written, uint64_t * global_index_arr, uint64_t * data_index_arr,
		                              uint64_t index_len);
int digital_rf_is_little_endian(void);


#endif
