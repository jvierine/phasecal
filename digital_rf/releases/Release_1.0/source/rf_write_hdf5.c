/* Implementation of rf_hdf5 library
 *
  See rf_hdf5.h for overview of this module.

  Written 2/2014 by B. Rideout

  $Id: rf_write_hdf5.c 387 2014-05-23 14:52:44Z brideout $
*/

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#include "digital_rf.h"



/* Public method implementations */
Digital_rf_write_object * digital_rf_create_write_hdf5(char * directory, hid_t dtype_id, uint64_t samples_per_file,
													  uint64_t files_per_directory, uint64_t global_start_sample,
													  double sample_rate, char * uuid_str,
		                                              int compression_level, int checksum, int is_complex,
		                                              int num_subchannels, int marching_dots)
/*  digital_rf_create_write_hdf5 returns an Digital_rf_write_object used to write a single channel of RF data to
 * a directory, or NULL with error to standard error if failure.
 *
 * Inputs:
 * 		char * directory - a directory under which to write the resultant Hdf5 files.  Must already exist.
 * 			Hdf5 files will be stored as YYYY-MM-DDTHH:MM:SS/rf@<unix_second>.<3 digit millsecond>.h5
 * 		hid_t dtype_id - data type id as defined by hdf5.h
 * 		uint64_t samples_per_file - the number of samples in a single Hdf5 file. Exactly that number of samples
 * 			will be in that file. samples_per_file will be saved in each Hdf5 file's metadata
 * 		uint64_t files_per_directory - the number of Hdf5 files in each directory of the form YYYY-MM-DDTHH:MM:SS.
 * 			A new directory will be created when that number is reached.
 * 		uint64_t global_start_sample - The start time of the first sample in units of samples since UT midnight 1970-01-01.
 * 		double sample_rate - sample rate in Hz
 * 		char * uuid_str - a string containing a UUID generated for that channel.  uuid_str saved in
 * 			each resultant Hdf5 file's metadata.
 * 		int compression_level - if 0, no compression used.  If 1-9, level of gzip compression.  Higher compression
 * 			means smaller file size and more time used.
 * 		int checksum - if non-zero, HDF5 checksum used.  If 0, no checksum used.
 * 		int is_complex - 1 if complex (IQ) data, 0 if single-valued
 * 		int num_subchannels - the number of subchannels of complex or single valued data recorded at once.
 * 			Must be 1 or greater. Note: A single stream of complex values is one subchannel, not two.
 * 		int marching_dots - non-zero if marching dots desired when writing; 0 if not
 *
 * 	Hdf5 format
 *
 * 	/rf_data - dataset of size (samples_per_file, 1 + is_complex), datatype = dtype_id
 * 	/rf_data_index - dataset of size (num of separate block of data, 2), datatype - uint_64, length at least 1
 *  /rf_data has 7 attributes: sequence_num (int), samples_per_file (int), uuid_str (string), sample_rate (double),
 *  	is_complex (0 or 1 - int), num_subchannels (int), and computer_time (unix timestamp from computer when
 *  	file created) uint64_t
 */
{
	/* local variables */
	Digital_rf_write_object * hdf5_data_object;
	hsize_t  chunk_dims[2]  = {CHUNK_SIZE_RF_DATA_INDEX, 2};
	time_t     computer_time;

	/* allocate overall object */
	if ((hdf5_data_object = (Digital_rf_write_object *)malloc(sizeof(Digital_rf_write_object)))==0)
	{
		fprintf(stderr, "malloc failure - unrecoverable\n");
		exit(-1);
	}
	/* init everything to NULL that will be malloced */
	hdf5_data_object->directory = NULL;
	hdf5_data_object->sub_directory = NULL;
	hdf5_data_object->uuid_str = NULL;
	hdf5_data_object->dataset = 0;
	hdf5_data_object->dataset_prop = 0;
	hdf5_data_object->dataspace = 0;
	hdf5_data_object->filespace = 0;
	hdf5_data_object->memspace = 0;
	hdf5_data_object->hdf5_file = 0; /* indicates no Hdf5 file presently opened */
	hdf5_data_object->index_dataset = 0;
	hdf5_data_object->index_prop = 0;
	hdf5_data_object->next_index_avail = 0;

	/* this value not set until digital_rf_write_hdf5 called */
	hdf5_data_object->chunk_size = 0;

	/* set directory name */
	if ((hdf5_data_object->directory = (char *)malloc(sizeof(char) * (strlen(directory)+2)))==0)
	{
		fprintf(stderr, "malloc failure - unrecoverable\n");
		exit(-1);
	}
	strcpy(hdf5_data_object->directory, directory);
	if (hdf5_data_object->directory[strlen(hdf5_data_object->directory)] != '/')
		strcat(hdf5_data_object->directory, "/");

	if (digital_rf_check_hdf5_directory(hdf5_data_object->directory))
	{
		digital_rf_close_write_hdf5(hdf5_data_object);
		return(NULL);
	}

	/* set UUID */
	if ((hdf5_data_object->uuid_str = (char *)malloc(sizeof(char) * (strlen(uuid_str)+1)))==0)
	{
		fprintf(stderr, "malloc failure - unrecoverable\n");
		exit(-1);
	}
	strcpy(hdf5_data_object->uuid_str, uuid_str);

	if (compression_level < 0 || compression_level > 9)
	{
		fprintf(stderr, "Illegal compression level, must be 0-9\n");
		digital_rf_close_write_hdf5(hdf5_data_object);
		return(NULL);
	}

	if (num_subchannels < 1)
	{
		fprintf(stderr, "Illegal num_subchannels %i, must be greater than 0\n", num_subchannels);
		digital_rf_close_write_hdf5(hdf5_data_object);
		return(NULL);
	}
	hdf5_data_object->num_subchannels = num_subchannels;

	/* init other values in hdf5_data_object */
	hdf5_data_object->samples_per_file = samples_per_file;
	hdf5_data_object->files_per_directory = files_per_directory;
	hdf5_data_object->global_start_sample = global_start_sample;
	hdf5_data_object->sample_rate = sample_rate;
	hdf5_data_object->dtype_id = dtype_id;
	hdf5_data_object->global_index = 0; /* index of next sample to write */
	hdf5_data_object->present_seq = -1; /* present file seq number, -1 because no file open */
	hdf5_data_object->dataset_index = 0; /* where in the dataset the next write should occur */
	hdf5_data_object->dataset_avail = 0; /* how many samples are free to write to it the open Hdf5 file */
	hdf5_data_object->block_index = 0;   /* the next available row in the open Hdf5 file/rf_data_index dataset to write to */
	hdf5_data_object->marching_dots = marching_dots;

	/* init_utc_timestamp - stored as attribute to allow conversion to astronomical times */
	computer_time = time(NULL);
	hdf5_data_object->init_utc_timestamp = (uint64_t)computer_time;

	if (is_complex)
	{
		hdf5_data_object->is_complex = 1;
		hdf5_data_object->rank = 2;
		/* create complex compound datatype */
		hdf5_data_object->complex_dtype_id = H5Tcreate(H5T_COMPOUND, 2*H5Tget_size(hdf5_data_object->dtype_id));
		/* create r column */
		H5Tinsert (hdf5_data_object->complex_dtype_id, "r", 0, hdf5_data_object->dtype_id);
		/* create i column */
		H5Tinsert (hdf5_data_object->complex_dtype_id, "i", H5Tget_size(hdf5_data_object->dtype_id),
				   hdf5_data_object->dtype_id);
	}
	else
	{
		hdf5_data_object->is_complex = 0;
		if (hdf5_data_object->num_subchannels == 1)
			hdf5_data_object->rank = 1;
		else
			hdf5_data_object->rank = 2;
		hdf5_data_object->complex_dtype_id = (hid_t)0; /* make sure its not used by accident */
	}

	/* check for illegal values */
	if (samples_per_file == 0 || files_per_directory == 0)
	{
		fprintf(stderr, "Illegal samples_per_file or files_per_directory, must not be zero\n");
		digital_rf_close_write_hdf5(hdf5_data_object);
		return(NULL);
	}
	if (global_start_sample == 0)
	{
		fprintf(stderr, "Illegal global_start_sample, must not be zero\n");
		digital_rf_close_write_hdf5(hdf5_data_object);
		return(NULL);
	}
	if (sample_rate <= 0.0)
	{
		fprintf(stderr, "Illegal sample_rate, must be positive\n");
		digital_rf_close_write_hdf5(hdf5_data_object);
		return(NULL);
	}


	/* dataset_prop is constant except for chunk size, so we can start to set this up in init */
	hdf5_data_object->dataset_prop = H5Pcreate (H5P_DATASET_CREATE);
	if (compression_level != 0)
		H5Pset_deflate (hdf5_data_object->dataset_prop, compression_level);
	if (checksum)
		H5Pset_filter (hdf5_data_object->dataset_prop, H5Z_FILTER_FLETCHER32, 0, 0, NULL);
	if (checksum || compression_level != 0)
		hdf5_data_object->needs_chunking = 1;
	else
		hdf5_data_object->needs_chunking = 0;
	/* set fill value for data gaps according to input dtype_id */
	if (digital_rf_set_fill_value(hdf5_data_object))
	{
		digital_rf_close_write_hdf5(hdf5_data_object);
		return(NULL);
	}

	/* index_prop is constant so we can start to set this up in init */
	hdf5_data_object->index_prop = H5Pcreate (H5P_DATASET_CREATE);
	H5Pset_chunk (hdf5_data_object->index_prop, 2, chunk_dims);

	/* done - return object */
	return(hdf5_data_object);
}



int digital_rf_write_hdf5(Digital_rf_write_object *hdf5_data_object, uint64_t global_leading_edge_index, void * vector,
						  uint64_t vector_length)
/*
 * digital_rf_write_hdf5 writes a continuous block of data from vector into one or more Hdf5 files
 *
 * Inputs:
 * 		Digital_rf_write_object *hdf5_data_object - C struct created by digital_rf_create_write_hdf5
 * 		uint64_t global_leading_edge_index - index to write data to.  This is a global index with zero
 * 			representing the sample taken at the time global_start_sample specified in the init method.
 * 			Note that all values stored in Hdf5 file will have global_start_sample added, and this offset
 * 			should NOT be added by the user. Error raised and -1 returned if before end of last write.
 * 		void * vector - pointer into data vector to write
 * 		uint64_t vector_length - number of samples to write to Hdf5
 *
 * 	Affects - Writes data to existing open Hdf5 file.  May close that file and write some or all of remaining data to
 * 		new Hdf5 file.
 *
 * 	Returns 0 if success, non-zero and error written if failure.
 *
 */
{
	uint64_t data_index_arr[1] = {0};
	uint64_t index_len = 1;
	int result;

	/* just call full method digital_rf_write_blocks_hdf5 under the covers */
	result = digital_rf_write_blocks_hdf5(hdf5_data_object, &global_leading_edge_index, data_index_arr, index_len, vector, vector_length);
	return(result);

}


int digital_rf_write_blocks_hdf5(Digital_rf_write_object *hdf5_data_object, uint64_t * global_index_arr, uint64_t * data_index_arr,
		                         uint64_t index_len, void * vector, uint64_t vector_length)
/*
 * digital_rf_write_blocks_hdf5 writes blocks of data from vector into one or more Hdf5 files
 *
 * Inputs:
 * 		Digital_rf_write_object *hdf5_data_object - C struct created by digital_rf_create_write_hdf5
 * 		uint64_t * global_index_arr - an array of global indices into the samples being written.  The global
 * 			index is the total number of sample periods since data taking began, including gaps.  Note that
 * 			all values stored in Hdf5 file will have global_start_sample added, and this offset
 * 			should NOT be added by the user.  Error is raised if any value is before its expected value (meaning repeated data).
 * 		uint64_t * data_index_arr - an array of len = len(global_index_arr), where the indices are related to which
 * 			sample in the vector being passed in is being referred to in global_index_arr.  First value must be 0
 * 			or error raised.  Values must be increasing, and cannot be equal or greater than index_len or error raised.
 * 		uint_64 index_len - the len of both global_index_arr and data_index_arr.  Must be greater than 0.
 * 		void * vector - pointer into data vector to write
 * 		uint64_t vector_length - number of samples to write to Hdf5
 *
 * 	Affects - Writes data to existing open Hdf5 file.  May close that file and write some or all of remaining data to
 * 		new Hdf5 file.
 *
 * 	Returns 0 if success, non-zero and error written if failure.
 *
 */
{
	char error_str[SMALL_HDF5_STR] = "";
	uint64_t samples_written = 0; /* total samples written so far to all Hdf5 files during this write call */
	uint64_t dataset_samples_written = 0; /* number of samples written to the present file */
	hsize_t      chunk_dims[2] = {0, hdf5_data_object->num_subchannels};
	int chunk_size = 0;

	/* verify data exists */
	if (!vector)
	{
		sprintf(error_str, "Null data passed in\n");
		fprintf(stderr, "%s", error_str);
		return(-2);
	}

	/* verify not writing in the past */
	if (global_index_arr[0] < hdf5_data_object->global_index)
	{
		sprintf(error_str, "Request index %" PRIu64 " before first expected index %" PRIu64 " in digital_rf_write_hdf5\n",
				global_index_arr[0], hdf5_data_object->global_index);
		fprintf(stderr, "%s", error_str);
		return(-3);
	}

	/* set chunking if needed */
	if (hdf5_data_object->needs_chunking && !hdf5_data_object->chunk_size)
	{
		if (vector_length < hdf5_data_object->samples_per_file)
			chunk_size = vector_length;
		else
			chunk_size = hdf5_data_object->samples_per_file;
		hdf5_data_object->chunk_size = chunk_size;
		chunk_dims[0] = chunk_size;
		H5Pset_chunk (hdf5_data_object->dataset_prop, hdf5_data_object->rank, chunk_dims);
	}

	/* loop until all data written - this loop breaks multiple file writes into a series single file writes*/
	while (samples_written < vector_length)
	{
		dataset_samples_written = digital_rf_write_samples_to_file(hdf5_data_object, samples_written,
				global_index_arr, data_index_arr, index_len, vector, vector_length);
		if (dataset_samples_written == 0)
					return(-6);
		samples_written += dataset_samples_written;
	}


	return(0);
}


int digital_rf_close_write_hdf5(Digital_rf_write_object *hdf5_data_object)
/* digital_rf_close_write_hdf5 closes open Hdf5 file if needed and releases all memory associated with hdf5_data_object
 *
 * Inputs:
 * 		Digital_rf_write_object *hdf5_data_object - C struct created by digital_rf_create_write_hdf5
 */
{
	digital_rf_free_hdf5_data_object(hdf5_data_object);
	return(0);
}

int digital_rf_get_unix_time(uint64_t global_sample, double sample_rate, int * year, int * month, int *day,
		                     int * hour, int * minute, int * second, uint64_t * picosecond)
/* get_unix_time converts a global_sample and a sample rate into year, month, day
 * 	hour, minute, second, picosecond
 *
 * 	Returns 0 if success, -1 if failure.
 *
 */
{
	struct tm *gm;
	time_t unix_second;
	double unix_remainder;

	/* set values down to second using gmtime */
	unix_second = (time_t)(global_sample / sample_rate);
	gm = gmtime(&unix_second);
	if (gm == NULL)
		return(-1);
	*year = gm->tm_year + 1900;
	*month = gm->tm_mon + 1;
	*day = gm->tm_mday;
	*hour = gm->tm_hour;
	*minute = gm->tm_min;
	*second = gm->tm_sec;

	/* set picoseconds */
	if (fmod(sample_rate, 1.0) == 0.0) /* use integer logic when sample rate can be converted to an integer */
		unix_remainder = global_sample - (unix_second * (uint64_t)sample_rate);
	else
		unix_remainder = fmod((double)global_sample, sample_rate);
	*picosecond = (uint64_t)round((unix_remainder/sample_rate)*1.0E12);
	return(0);
}


/* Private Method implementations */

int digital_rf_free_hdf5_data_object(Digital_rf_write_object *hdf5_data_object)
/* digital_rf_free_hdf5_data_object frees all resources in hdf5_data_object */
{
	if (hdf5_data_object->directory != NULL)
		free(hdf5_data_object->directory);
	if (hdf5_data_object->sub_directory != NULL)
		free(hdf5_data_object->sub_directory);
	if (hdf5_data_object->uuid_str != NULL)
		free(hdf5_data_object->uuid_str);

	/* free all Hdf5 resources */
	if (hdf5_data_object->dataset)
		H5Dclose (hdf5_data_object->dataset);
	if (hdf5_data_object->dataset_prop)
		H5Pclose (hdf5_data_object->dataset_prop);
	if (hdf5_data_object->dataspace)
		H5Sclose (hdf5_data_object->dataspace);
	if (hdf5_data_object->filespace)
		H5Sclose (hdf5_data_object->filespace);
	if (hdf5_data_object->memspace)
		H5Sclose (hdf5_data_object->memspace);
	if (hdf5_data_object->index_dataset)
		H5Dclose (hdf5_data_object->index_dataset);
	if (hdf5_data_object->index_prop)
		H5Pclose (hdf5_data_object->index_prop);
	if (hdf5_data_object->hdf5_file)
		H5Fclose (hdf5_data_object->hdf5_file);
	free(hdf5_data_object);

	return(0);
}

int digital_rf_check_hdf5_directory(char * directory)
/* digital_rf_check_hdf5_directory checks if directory exists.  If it does not
 * exist, returns -1.  If okay, returns 0.
 *
 */
{
	/* local variable */
	struct stat stat_obj = {0}; /* used to run stat to determine if directory exists */
	char error_str[BIG_HDF5_STR] = "";

	/* see if directory needs to be created */
	if (stat(directory, &stat_obj))
	{
		strncat(error_str, directory, 200);
		strcat(error_str, " does not exist.\n");
		fprintf(stderr, "%s", error_str);
		return(-1);
	} else {
		/* verify its a directory */
		if(!S_ISDIR(stat_obj.st_mode))
		{
			sprintf(error_str, "The following is not a directory: %s\n", directory);
			fprintf(stderr, "%s", error_str);
			return(-1);
		}
	}
	return(0);
}


uint64_t digital_rf_write_samples_to_file(Digital_rf_write_object *hdf5_data_object, uint64_t samples_written, uint64_t * global_index_arr,
		uint64_t * data_index_arr, uint64_t index_len, void * vector, uint64_t vector_length)
/* digital_rf_write_samples_to_file writes some or all of the data to a single Hdf5 file.
 *
 * Returns the number of samples written. Returns 0 if error.
 *
 * Inputs:
 * 	Digital_rf_write_object *hdf5_data_object - the Digital_rf_write_object created by digital_rf_create_write_hdf5
 * 	uint64_t samples_written - the number of samples written to previous files during this particular user write call
 * 	uint64_t * global_index_arr - an array of global indices into the samples being written.  The global
 * 		index is the total number of sample periods since data taking began, including gaps.  Error
 * 		is raised if any value is before its expected value (meaning repeated data).
 * 	uint64_t * data_index_arr - an array of len = len(global_index_arr), where the indices are related to which
 * 		sample in the vector being passed in by the user, so that the first value is always 0,
 * 		or error raised.  Values must be increasing, and cannot be equal or greater than index_len or error raised.
 * 	uint_64 index_len - the len of both global_index_arr and data_index_arr.  Must be greater than 0.
 * 	void * vector - the vector (either single valued or complex) containing the data, or NULL pointer if no data to be written.
 * 	uint64_t vector_length - the total number of samples to write from vector.
 *
 */
{
	/* local variables */
	uint64_t samples_left_to_write, samples_to_write, last_global_index, samples_after_last_global_index;
	uint64_t next_global_index;
	int block_index_len; /* len of /rf_data_index dataset needed for this particular write */
	uint64_t * rf_data_index_arr; /* will be malloced and filled out with all data needed for rf_data_index table */
	int result;
	hsize_t size[2] = {0, hdf5_data_object->num_subchannels}; /* will be set to size of full dataset in file */
	hsize_t offset[2] = {0,0};        /* will be set to the index in file writing to */
	herr_t  status;                   /* Hdf5 error status */

	/* verify inputs are sensible */
	if (index_len < 1)
	{
		fprintf(stderr, "Illegal index_len %" PRIu64 " in digital_rf_write_samples_to_file\n", index_len);
		return(0);
	}
	if (data_index_arr[0] != 0)
	{
		fprintf(stderr, "Illegal first value %" PRIu64 " in data_index_arr, must be 0\n", data_index_arr[0]);
				return(0);
	}

	/* get all the info needed to create (or expand) /rf_data_index and fill it out */
	rf_data_index_arr = digital_rf_create_rf_data_index(hdf5_data_object, samples_written, global_index_arr,
				data_index_arr, index_len, &block_index_len);
	if (!rf_data_index_arr && block_index_len == -1)
		return(0);

    /* First we need to see if we need to open a new file */
	if (hdf5_data_object->hdf5_file == 0)
	{
		next_global_index = digital_rf_get_global_sample(samples_written, global_index_arr, data_index_arr, index_len);
		result = digital_rf_create_hdf5_file(hdf5_data_object, next_global_index);
		if (result)
		{
			free(rf_data_index_arr);
			return(0);
		}
	}

	samples_left_to_write = vector_length - samples_written;
	if (samples_left_to_write > hdf5_data_object->dataset_avail)
		samples_to_write = hdf5_data_object->dataset_avail;
	else
		samples_to_write = samples_left_to_write;

	/* create dataspace hyperslab to write to */
	if (hdf5_data_object->filespace)
		H5Sclose (hdf5_data_object->filespace);
	hdf5_data_object->filespace = H5Dget_space(hdf5_data_object->dataset);
	offset[0] = hdf5_data_object->dataset_index;
	size[0] = samples_to_write;
	H5Sselect_hyperslab(hdf5_data_object->filespace, H5S_SELECT_SET,
						offset, NULL, size, NULL);

	/* set up memspace to control write */
	if (hdf5_data_object->memspace)
		H5Sclose (hdf5_data_object->memspace);
	hdf5_data_object->memspace = H5Screate_simple(hdf5_data_object->rank, size, NULL);

	/* write rf_data */
	if (hdf5_data_object->is_complex == 0)
		status = H5Dwrite(hdf5_data_object->dataset, hdf5_data_object->dtype_id, hdf5_data_object->memspace,
						  hdf5_data_object->filespace, H5P_DEFAULT,
						  (char *)vector + (samples_written * H5Tget_size(hdf5_data_object->dtype_id) * hdf5_data_object->num_subchannels));
	else /* complex */
		status = H5Dwrite(hdf5_data_object->dataset, hdf5_data_object->complex_dtype_id, hdf5_data_object->memspace,
						  hdf5_data_object->filespace, H5P_DEFAULT,
						  (char *)vector + (samples_written * H5Tget_size(hdf5_data_object->dtype_id) * 2* hdf5_data_object->num_subchannels));

	if (status < 0)
	{
		fprintf(stderr, "Failure at H5DWrite\n");
		free(rf_data_index_arr);
		return(0);
	}

	/* write rf_data_index dataset */
	if (block_index_len > 0)
	{
		if (digital_rf_write_rf_data_index(hdf5_data_object, rf_data_index_arr, block_index_len))
		{
			free(rf_data_index_arr);
			return(0);
		}
	}

	/* advance state */
	hdf5_data_object->dataset_index += samples_to_write;
	hdf5_data_object->dataset_avail -= samples_to_write;
	/* update global_index - see meaning of rf_data_index_arr and block_index_len for logic.
	 * Basically using last line and counting extra samples after last global index */
	if (block_index_len > 0)
	{
		last_global_index = rf_data_index_arr[block_index_len*2 - 2] - hdf5_data_object->global_start_sample;
		samples_after_last_global_index = samples_to_write - (rf_data_index_arr[block_index_len*2 - 1] - rf_data_index_arr[1]);
		hdf5_data_object->global_index = last_global_index + samples_after_last_global_index;
		free(rf_data_index_arr);
	}
	else
	{
		/* if no indices were used, then all data written must have been continuous */
		hdf5_data_object->global_index += samples_to_write;
	}

	assert(hdf5_data_object->dataset_index <= hdf5_data_object->samples_per_file); /* debug */

	if (hdf5_data_object->dataset_index == hdf5_data_object->samples_per_file)
	{
		/* hdf5 file full - close it */
		H5Dclose (hdf5_data_object->dataset);
		hdf5_data_object->dataset = 0;
		H5Dclose (hdf5_data_object->index_dataset);
		hdf5_data_object->index_dataset = 0;
		H5Sclose (hdf5_data_object->dataspace);
		hdf5_data_object->dataspace = 0;
		if (hdf5_data_object->filespace)
		{
			H5Sclose (hdf5_data_object->filespace);
			hdf5_data_object->filespace = 0;
		}
		if (hdf5_data_object->memspace)
		{
			H5Sclose (hdf5_data_object->memspace);
			hdf5_data_object->memspace = 0;
		}
		H5Fclose (hdf5_data_object->hdf5_file);
		hdf5_data_object->hdf5_file = 0;
		hdf5_data_object->dataset_index = 0;
	}

	return(samples_to_write);
}


int digital_rf_create_hdf5_file(Digital_rf_write_object *hdf5_data_object, uint64_t next_global_sample)
/* digital_rf_create_hdf5_file opens a new Hdf5 file
 *
 * Inputs:
 * 	Digital_rf_write_object *hdf5_data_object - the Digital_rf_write_object created by digital_rf_create_write_hdf5
 * 	uint64_t next_global_sample - global index of next sample to write (used to create file and directory names)
 *
 * 	Creates a file with /rf_data dataset of size (samples_per_file, 2)
 *
 * 	Returns 0 if success, -1 if failure
 *
 */
{
	/* local variables */
	char datasetname[] = "rf_data";
	char fullname[BIG_HDF5_STR] = "";
	char basename[SMALL_HDF5_STR] = "";
	char error_str[BIG_HDF5_STR] = "";
	double unix_timestamp;
	uint64_t num_rows = 0;
	hsize_t  dims[2]  = {0, hdf5_data_object->num_subchannels};
	hsize_t  maxdims[2] = {0, hdf5_data_object->num_subchannels};

    if (hdf5_data_object->marching_dots)
    {
		printf(".");
		fflush(stdout);
    }

	unix_timestamp = (next_global_sample + hdf5_data_object->global_start_sample)/hdf5_data_object->sample_rate;

	/* get fullname of new Hdf5 file */
	hdf5_data_object->present_seq++; /* indicates the creation of a new file */

	/* create new directory if needed */
	if (hdf5_data_object->present_seq % hdf5_data_object->files_per_directory == 0)
		if (digital_rf_create_new_directory(hdf5_data_object, next_global_sample))
			return(-1);

	strcpy(fullname, hdf5_data_object->directory); /* previous check ensures these three commands succeed */
	strcat(fullname, hdf5_data_object->sub_directory);
	sprintf(basename, "rf@%011.3f.h5", unix_timestamp);
	strcat(fullname, basename);

	/* Create a new file. If file exists will fail. */
	hdf5_data_object->hdf5_file = H5Fcreate (fullname, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
	if (hdf5_data_object->hdf5_file < 0)
	{
		sprintf(error_str, "The following Hdf5 file could not be created, or already exists: %s\n", fullname);
		fprintf(stderr, "%s", error_str);
		hdf5_data_object->hdf5_file = 0;
		return(-1);
	}

	/* now we add the dataset to create */
	num_rows = hdf5_data_object->samples_per_file;

	dims[0] = num_rows;
	maxdims[0] = hdf5_data_object->samples_per_file;
	/* Create the data space with set dimensions. */
	if (hdf5_data_object->dataspace)
			H5Sclose (hdf5_data_object->dataspace);
	hdf5_data_object->dataspace = H5Screate_simple (hdf5_data_object->rank, dims, maxdims);
	/* create the dataset */
	if (hdf5_data_object->dataset)
			H5Dclose (hdf5_data_object->dataset);

	if (hdf5_data_object->is_complex == 0)
		hdf5_data_object->dataset = H5Dcreate2 (hdf5_data_object->hdf5_file, datasetname,
												hdf5_data_object->dtype_id,
												hdf5_data_object->dataspace, H5P_DEFAULT,
												hdf5_data_object->dataset_prop, H5P_DEFAULT);
	else
		hdf5_data_object->dataset = H5Dcreate2 (hdf5_data_object->hdf5_file, datasetname,
												hdf5_data_object->complex_dtype_id,
												hdf5_data_object->dataspace, H5P_DEFAULT,
												hdf5_data_object->dataset_prop, H5P_DEFAULT);

	hdf5_data_object->dataset_index = 0;        /* next write will be to first row */
	hdf5_data_object->dataset_avail = num_rows; /* size available to next write */

	/* last we add metadata */
	digital_rf_write_metadata(hdf5_data_object);
	return(0);
}


int digital_rf_create_new_directory(Digital_rf_write_object *hdf5_data_object, uint64_t next_global_sample)
/* digital_rf_create_new_directory creates a new subdirectory to store Hdf5 files in
 *
 * Subdirectory will be named <hdf5_data_object->directory>/YYYY-MM-DDTHH:MM:SS,
 * where time is set by next_global_sample to write
 *
 * Affects: sets hdf5_data_object->sub_directory.
 *
 * Returns 0 if success, -1 if failure. Fails if this directory already exists or can't be written.
 *
 */
{
	/* local variables */
	int year, month, day, hour, minute, second;
	uint64_t picosecond;
	char sub_directory[SMALL_HDF5_STR] = "";
	char full_directory[BIG_HDF5_STR] = "";

	if (digital_rf_get_unix_time(next_global_sample + hdf5_data_object->global_start_sample, hdf5_data_object->sample_rate,
			                     &year, &month, &day, &hour, &minute, &second, &picosecond))
		return(-1);

	sprintf(sub_directory, "%04i-%02i-%02iT%02i:%02i:%02i", year, month, day, hour, minute, second);
	strcpy(full_directory, hdf5_data_object->directory); /* directory ends with "/" */
	strcat(full_directory, sub_directory);

	if (mkdir(full_directory, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH))
	{
		fprintf(stderr, "Unable to create directory %s\n", sub_directory);
		return(-1);
	}

	if (hdf5_data_object->sub_directory != NULL)
		free(hdf5_data_object->sub_directory);

	if ((hdf5_data_object->sub_directory = (char *)malloc(sizeof(char) * (strlen(sub_directory)+2)))==0)
	{
		fprintf(stderr, "malloc failure - unrecoverable\n");
		exit(-1);
	}
	strcpy(hdf5_data_object->sub_directory, sub_directory);
	strcat(hdf5_data_object->sub_directory, "/"); /* will always end with "/" */
	return(0);
}



int digital_rf_set_fill_value(Digital_rf_write_object *hdf5_data_object)
/* digital_rf_set_fill_value sets the fill value property in hdf5_data_object->dataset_prop according to dtype_id.
 *
 * Returns 0 if success, -1 if failure
 */
{
	/* local variables */

	/* integer minimum values, first value is matching byte order, second is reversed byte order */

	int minUnsignedInt = 0;

	/* char */
	char minChar = -128;
	struct complex_char_fill_type {char r,i;} complex_char_fill;
	complex_char_fill.r = minChar;
	complex_char_fill.i = minChar;
	struct complex_uchar_fill_type {unsigned char r,i;} complex_uchar_fill;
	complex_uchar_fill.r = 0;
	complex_uchar_fill.i = 0;

	/* short */
	short minShort[2] = {INT16_MIN, 128};
	struct complex_short_fill_type {short r,i;} complex_short_fill[2];
	complex_short_fill[0].r = minShort[0];
	complex_short_fill[0].i = minShort[0];
	complex_short_fill[1].r = minShort[1];
	complex_short_fill[1].i = minShort[1];
	struct complex_ushort_fill_type {unsigned short r,i;} complex_ushort_fill;
	complex_ushort_fill.r = 0;
	complex_ushort_fill.i = 0;

	/* int */
	int minInt[2] = {INT32_MIN, 128};
	struct complex_int_fill_type {int r,i;} complex_int_fill[2];
	complex_int_fill[0].r = minInt[0];
	complex_int_fill[0].i = minInt[0];
	complex_int_fill[1].r = minInt[1];
	complex_int_fill[1].i = minInt[1];
	struct complex_uint_fill_type {unsigned int r,i;} complex_uint_fill;
	complex_uint_fill.r = 0;
	complex_uint_fill.i = 0;

	/* int64 */
	int64_t minLLong[2] = {INT64_MIN,128};
	struct complex_long_fill_type {int64_t r,i;} complex_long_fill[2];
	complex_long_fill[0].r = minLLong[0];
	complex_long_fill[0].i = minLLong[0];
	complex_long_fill[1].r = minLLong[1];
	complex_long_fill[1].i = minLLong[1];
	struct complex_ulong_fill_type {uint64_t r,i;} complex_ulong_fill;
	complex_ulong_fill.r = 0;
	complex_ulong_fill.i = 0;

	// float
	float float_fill = NAN;
	struct complex_float_fill_type {float r,i;} complex_float_fill;
	complex_float_fill.r = float_fill;
	complex_float_fill.i = float_fill;

	// double
	double double_fill = NAN;
	struct complex_double_fill_type {double r,i;} complex_double_fill;
	complex_double_fill.r = double_fill;
	complex_double_fill.i = double_fill;

	char error_str[SMALL_HDF5_STR] = "";
	H5T_class_t classType;
	H5T_sign_t signType;
	int numBytes;
	int endian_flip, write_endian;

	/* found out if the output endian differs from the host endian */
	endian_flip = 0;
	write_endian = H5Tget_order(hdf5_data_object->dtype_id);
	if (digital_rf_is_little_endian() && (write_endian == H5T_ORDER_BE))
		endian_flip = 1;
	else if ((!digital_rf_is_little_endian()) && (write_endian == H5T_ORDER_LE))
		endian_flip = 1;

	/* find out whether we are using integers or floats */
	classType = H5Tget_class( hdf5_data_object->dtype_id );
	/* find out if integer is signed, and the number of bytes */
	signType = H5Tget_sign(hdf5_data_object->dtype_id);
	numBytes = H5Tget_size(hdf5_data_object->dtype_id);

	if (classType == H5T_FLOAT && hdf5_data_object->is_complex == 0)
	{
		H5Pset_fill_value(hdf5_data_object->dataset_prop, hdf5_data_object->dtype_id, &double_fill);
	}
	else if (classType == H5T_FLOAT && hdf5_data_object->is_complex != 0 && numBytes == 4)
	{
		H5Pset_fill_value(hdf5_data_object->dataset_prop, hdf5_data_object->complex_dtype_id, &complex_float_fill);
	}
	else if (classType == H5T_FLOAT && hdf5_data_object->is_complex != 0 && numBytes == 8)
	{
		H5Pset_fill_value(hdf5_data_object->dataset_prop, hdf5_data_object->complex_dtype_id, &complex_double_fill);
	}
	else if (classType == H5T_INTEGER)
	{
		if (hdf5_data_object->is_complex == 0)
		{
			if (signType == H5T_SGN_NONE)
				H5Pset_fill_value(hdf5_data_object->dataset_prop, hdf5_data_object->dtype_id, &minUnsignedInt);
			else
			{
				switch(numBytes)
				{
					case 1:
						H5Pset_fill_value(hdf5_data_object->dataset_prop, hdf5_data_object->dtype_id, &minChar);
						break;
					case 2:
						H5Pset_fill_value(hdf5_data_object->dataset_prop, hdf5_data_object->dtype_id, &minShort[endian_flip]);
						break;
					case 4:
						H5Pset_fill_value(hdf5_data_object->dataset_prop, hdf5_data_object->dtype_id, &minInt[endian_flip]);
						break;
					case 8:
						H5Pset_fill_value(hdf5_data_object->dataset_prop, hdf5_data_object->dtype_id, &minLLong[endian_flip]);
						break;
					default:
						sprintf(error_str, "Integer type has unexpected number of bytes: %i\n", numBytes);
						fprintf(stderr, "%s", error_str);
						return(-1);
				}
			}
		}
		else /* complex data */
		{
			switch(numBytes)
			{
				case 1:
					if (signType == H5T_SGN_NONE)
						H5Pset_fill_value(hdf5_data_object->dataset_prop, hdf5_data_object->complex_dtype_id, &complex_uchar_fill);
					else
						H5Pset_fill_value(hdf5_data_object->dataset_prop, hdf5_data_object->complex_dtype_id, &complex_char_fill);
					break;
				case 2:
					if (signType == H5T_SGN_NONE)
						H5Pset_fill_value(hdf5_data_object->dataset_prop, hdf5_data_object->complex_dtype_id,
								&complex_ushort_fill);
					else
						H5Pset_fill_value(hdf5_data_object->dataset_prop, hdf5_data_object->complex_dtype_id,
								&complex_short_fill[endian_flip]);
					break;
				case 4:
					if (signType == H5T_SGN_NONE)
						H5Pset_fill_value(hdf5_data_object->dataset_prop, hdf5_data_object->complex_dtype_id,
								&complex_uint_fill);
					else
						H5Pset_fill_value(hdf5_data_object->dataset_prop, hdf5_data_object->complex_dtype_id,
								&complex_int_fill[endian_flip]);
					break;
				case 8:
					if (signType == H5T_SGN_NONE)
						H5Pset_fill_value(hdf5_data_object->dataset_prop, hdf5_data_object->complex_dtype_id,
								&complex_ulong_fill);
					else
						H5Pset_fill_value(hdf5_data_object->dataset_prop, hdf5_data_object->complex_dtype_id,
								&complex_long_fill[endian_flip]);
					break;
				default:
					sprintf(error_str, "Integer type has unexpected number of bytes: %i\n", numBytes);
					fprintf(stderr, "%s", error_str);
					return(-1);
			}
		}

	}
	else
	{
		fprintf(stderr, "Hdf5 datatype passed into dtype_id is neither integer nor float - aborting\n");
		return(-1);
	}
	return(0);

}


void digital_rf_write_metadata(Digital_rf_write_object *hdf5_data_object)
/* digital_rf_write_metadata writes the following metadata to the open file:
 * sequence, samples_per_file, uuid_str
 *
 * Inputs:
 * 	Digital_rf_write_object *hdf5_data_object - the Digital_rf_write_object created by digital_rf_create_write_hdf5
 *
 */
{
	/* local variables */
	hid_t       attribute_id, dataspace_id;  /* identifiers */
	hid_t       str_dataspace, str_type, str_attribute;
	hsize_t     dims = 1;
	time_t      computer_time;
	int64_t    u_computer_time;

	dataspace_id = H5Screate_simple(1, &dims, NULL);

	/* sequence_num */
	attribute_id = H5Acreate2 (hdf5_data_object->dataset, "sequence_num", H5T_NATIVE_INT, dataspace_id,
								 H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attribute_id, H5T_NATIVE_INT, &(hdf5_data_object->present_seq));
	H5Aclose(attribute_id);

	/* num_subchannels */
	attribute_id = H5Acreate2 (hdf5_data_object->dataset, "num_subchannels", H5T_NATIVE_INT, dataspace_id,
								 H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attribute_id, H5T_NATIVE_INT, &(hdf5_data_object->num_subchannels));
	H5Aclose(attribute_id);

	/* is_complex */
	attribute_id = H5Acreate2 (hdf5_data_object->dataset, "is_complex", H5T_NATIVE_INT, dataspace_id,
								 H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attribute_id, H5T_NATIVE_INT, &(hdf5_data_object->is_complex));
	H5Aclose(attribute_id);

	/* samples_per_file */
	attribute_id = H5Acreate2 (hdf5_data_object->dataset, "samples_per_file", H5T_NATIVE_ULLONG, dataspace_id,
								 H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attribute_id, H5T_NATIVE_ULLONG, &(hdf5_data_object->samples_per_file));
	H5Aclose(attribute_id);

	/* sample_rate */
	attribute_id = H5Acreate2 (hdf5_data_object->dataset, "sample_rate", H5T_NATIVE_DOUBLE, dataspace_id,
								 H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &(hdf5_data_object->sample_rate));
	H5Aclose(attribute_id);

	/* init_utc_timestamp */
	attribute_id = H5Acreate2 (hdf5_data_object->dataset, "init_utc_timestamp", H5T_NATIVE_ULLONG, dataspace_id,
								 H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attribute_id, H5T_NATIVE_ULLONG, &(hdf5_data_object->init_utc_timestamp));
	H5Aclose(attribute_id);


	/* computer time */
	attribute_id = H5Acreate2 (hdf5_data_object->dataset, "computer_time", H5T_NATIVE_ULLONG, dataspace_id,
								 H5P_DEFAULT, H5P_DEFAULT);
	computer_time = time(NULL);
	u_computer_time = (uint64_t)computer_time;
	H5Awrite(attribute_id, H5T_NATIVE_ULLONG, &(u_computer_time));
	H5Aclose(attribute_id);

	/* uuid_str */
	str_dataspace  = H5Screate(H5S_SCALAR);
    str_type = H5Tcopy(H5T_C_S1);
	H5Tset_size(str_type, strlen(hdf5_data_object->uuid_str)+1);
    str_attribute = H5Acreate2(hdf5_data_object->dataset, "uuid_str", str_type, str_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(str_attribute, str_type, hdf5_data_object->uuid_str);
    H5Aclose(str_attribute);

    /* epoch */
	H5Tset_size(str_type, strlen(DIGITAL_RF_EPOCH)+1);
	str_attribute = H5Acreate2(hdf5_data_object->dataset, "epoch", str_type, str_dataspace, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(str_attribute, str_type, DIGITAL_RF_EPOCH);
	H5Aclose(str_attribute);

	/* digital_rf_time_description */
	H5Tset_size(str_type, strlen(DIGITAL_RF_TIME_DESCRIPTION)+1);
	str_attribute = H5Acreate2(hdf5_data_object->dataset, "digital_rf_time_description", str_type, str_dataspace, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(str_attribute, str_type, DIGITAL_RF_TIME_DESCRIPTION);
	H5Aclose(str_attribute);

	/* digital_rf_version */
	H5Tset_size(str_type, strlen(DIGITAL_RF_VERSION)+1);
	str_attribute = H5Acreate2(hdf5_data_object->dataset, "digital_rf_version", str_type, str_dataspace, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(str_attribute, str_type, DIGITAL_RF_VERSION);
	H5Aclose(str_attribute);

    /* free resources used */
	H5Tclose(str_type);
    H5Sclose(dataspace_id);
    H5Sclose(str_dataspace);


}


uint64_t * digital_rf_create_rf_data_index(Digital_rf_write_object *hdf5_data_object, uint64_t samples_written, uint64_t * global_index_arr,
			uint64_t * data_index_arr, uint64_t index_len, int * rows_to_write)
/* digital_rf_create_rf_data_index returns a malloced block of rf_data_index index data to write into the existing Hdf5 file
 * also sets the number of rows to be written.  Number of rows to be written may be zero, in which case returns NULL.
 *
 *  Digital_rf_write_object *hdf5_data_object - the Digital_rf_write_object created by digital_rf_create_write_hdf5
 * 	uint64_t samples_written - the number of samples written to previous files during this particular user write call
 * 	uint64_t * global_index_arr - an array of global indices into the samples being written.  The global
 * 		index is the total number of sample periods since data taking began, including gaps.  Error
 * 		is raised if any value is before its expected value (meaning repeated data).
 * 	uint64_t * data_index_arr - an array of len = len(global_index_arr), where the indices are related to which
 * 		sample in the vector being passed in by the user, so that the first value is always 0,
 * 		or error raised.  Values must be increasing, and cannot be equal or greater than index_len or error raised.
 * 	uint_64 index_len - the len of data_index_arr.  Must be greater than 0.
 * 	int * rows_to_write - this int will be set to the number of rows malloced in returned data
 *
 * 	Returns a malloced uint64_t array of size (rows_to_write, 2) with columns
 * 		dataset for this particular file, or -1 if an error detected. Used to allocate or increase size of /rf_data_index.
 * 		Returned array adds hdf5_data_object->global_start_sample to all global indices so that all are zero at 0UT 1970-01-01
 * 		May be NULL if rows_to_write is 0. Caller is responsible for freeing this.
 * 		Returns NULL and rows_to_write = -1 and error printed to stderr if error detected
 */
{
	int i; /* block loop variable */
	uint64_t this_index;
	uint64_t first_index, end_index; /* index of data samples to look over */
	int first_index_found = 0; /* bool that says whether the first possible index for this file was found */
	int row_count = 0;
	uint64_t last_index = 0; /* make sure indices are increasing at least as much as global_index_arr */
	char error_str[BIG_HDF5_STR] = "";
	uint64_t * ret_arr; /* will hold malloced data to be returned */
	int rows_written = 0; /* keeps tracks of rows written */

	/* figure out the data indices that could possibly be written in the present file */
	first_index = samples_written;
	end_index = (first_index + hdf5_data_object->samples_per_file) - hdf5_data_object->dataset_index;

	/* this first pass is just to count the number of rows needed before the malloc, and to valid data is reasonable */
	if (samples_written == 0 && global_index_arr[0] < hdf5_data_object->global_index)
	{
		sprintf(error_str, "global_index_arr passed in %" PRIu64 " before minimum value of %" PRIu64 "\n",
				global_index_arr[0], hdf5_data_object->global_index);
		fprintf(stderr, "%s", error_str);
		*rows_to_write = -1;
		return(NULL);
	}
	for (i=0; i<index_len; i++)
	{
		this_index = data_index_arr[i];

		/* more input data sanity checks */
		if (i>0 && last_index >= this_index)
		{
			sprintf(error_str, "indices in data_index_arr out of order - index %i and %i\n", i-1,i);
			fprintf(stderr, "%s", error_str);
			*rows_to_write = -1;
			return(NULL);
		}
		if (i>0 && ((this_index - last_index) > (global_index_arr[i] - global_index_arr[i-1])))
		{
			sprintf(error_str, "error - indices advancing faster than global index at index %i, illegal\n", i);
			fprintf(stderr, "%s", error_str);
			*rows_to_write = -1;
			return(NULL);
		}

		/* these indices okay - see if they are needed for this file */
		if (this_index == first_index)
			first_index_found = 1;
		if (this_index >= first_index && this_index < end_index)
			/* ignore first_index if in the middle of a file and global_index indicates no gap */
			if (!(this_index == first_index && hdf5_data_object->dataset_index > 0
					&& hdf5_data_object->global_index == global_index_arr[i]))
				row_count++; /* this_index needs to be written in this file */
		last_index = this_index;
	}
	if (!first_index_found)
		row_count++; /* we will need to insert an extra index at the beginning of the write */

	/* if no indices need to be malloced, return now */
	if (row_count == 0)
	{
		*rows_to_write = 0;
		return(NULL);
	}

	/* now that we know how many rows to malloc, malloc them */
	/* allocate overall object */
	if ((ret_arr = (uint64_t *)malloc(sizeof(uint64_t)*row_count*2))==0)
	{
		fprintf(stderr, "malloc failure - unrecoverable\n");
		exit(-1);
	}

	/* next pass is to fill out ret_arr */
	for (i=0; i<index_len; i++)
	{
		this_index = data_index_arr[i];
		/* ignore first_index if in the middle of a file and global_index indicates no gap */
		if (this_index == first_index && hdf5_data_object->dataset_index > 0
				&& hdf5_data_object->global_index == global_index_arr[i])
			continue;
		if (this_index >= first_index && this_index < end_index)
		{
			if (this_index == first_index)
			{
				ret_arr[0] = global_index_arr[i] + hdf5_data_object->global_start_sample; /* global index */
				ret_arr[1] = data_index_arr[i] + hdf5_data_object->dataset_index - samples_written; /* convert from data -> dataset index */
				rows_written++;
				continue;
			}
			else if (rows_written == 0)
			{
				/* we need to write first row since at a file boundary */
				assert(i!=0); /* or else there's a bug in my logic */
				ret_arr[0] =  hdf5_data_object->global_index + hdf5_data_object->global_start_sample;
				ret_arr[1] = 0;
				rows_written++;
			}
			ret_arr[rows_written*2] = global_index_arr[i] + hdf5_data_object->global_start_sample; /* global index */
			ret_arr[rows_written*2 + 1] = data_index_arr[i] + hdf5_data_object->dataset_index - samples_written; /* convert from data -> dataset index */
			rows_written++;
		}
	}
	if (rows_written == 0)
	{
		/* we need to write first row since at a file boundary */
		ret_arr[0] =  hdf5_data_object->global_index + hdf5_data_object->global_start_sample;
		ret_arr[1] = 0;
		rows_written++;
	}
	assert(rows_written==row_count);  /* or else there's a bug in my logic */
	*rows_to_write = rows_written;

	return(ret_arr);

}


int digital_rf_write_rf_data_index(Digital_rf_write_object * hdf5_data_object, uint64_t * rf_data_index_arr, int block_index_len)
/* digital_rf_write_rf_data_index writes rf_data_index to the open Hdf5 file
 *
 * Inputs:
 *  Digital_rf_write_object *hdf5_data_object - the Digital_rf_write_object created by digital_rf_create_write_hdf5
 *  uint64_t * rf_data_index_arr - uint64_t array of size (block_index_len * 2) to write - all values are already set
 *
 *  Returns 0 if success, non-zero if error
 */
{
	/* variables for /rf_data_index */
	char index_datasetname[] = "rf_data_index";
	hsize_t  index_dims[2]  = {0, 2};
	hsize_t  dimsext[2] = {block_index_len, 2};
	hsize_t  index_maxdims[2] = {H5S_UNLIMITED, 2};
	hsize_t  offset[2] = {0, 0};
	hid_t    index_dataspace, filespace, memspace;
	herr_t      status;

	/* find out if we need to create a new dataset, or expand and existing one */
	if (hdf5_data_object->index_dataset == 0)
	{
		/* create new dataset */
		index_dims[0] = block_index_len;
		index_dataspace = H5Screate_simple (2, index_dims, index_maxdims);
		hdf5_data_object->index_dataset = H5Dcreate2 (hdf5_data_object->hdf5_file, index_datasetname,
												H5T_NATIVE_ULLONG,
												index_dataspace, H5P_DEFAULT,
												hdf5_data_object->index_prop, H5P_DEFAULT);
		status = H5Dwrite(hdf5_data_object->index_dataset, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT,
					rf_data_index_arr);
		if (status < 0)
			return(status);

		H5Sclose(index_dataspace);
		hdf5_data_object->next_index_avail = block_index_len;
	}
	else
	{
		/* expand dataset */
		index_dims[0] = hdf5_data_object->next_index_avail + block_index_len;
		status = H5Dset_extent (hdf5_data_object->index_dataset, index_dims);
		filespace = H5Dget_space (hdf5_data_object->index_dataset);
		offset[0] = hdf5_data_object->next_index_avail;
		status = H5Sselect_hyperslab (filespace, H5S_SELECT_SET, offset, NULL,
				dimsext, NULL);
		if (status < 0)
			return(status);
		memspace = H5Screate_simple (2, dimsext, NULL);
		status = H5Dwrite (hdf5_data_object->index_dataset, H5T_NATIVE_ULLONG, memspace, filespace,
		                       H5P_DEFAULT, rf_data_index_arr);
		if (status < 0)
			return(status);

		H5Sclose (memspace);
		H5Sclose (filespace);
		hdf5_data_object->next_index_avail += block_index_len;
	}
	return(0);
}


uint64_t digital_rf_get_global_sample(uint64_t samples_written, uint64_t * global_index_arr, uint64_t * data_index_arr,
		                              uint64_t index_len)
/* digital_rf_get_global_sample calculates the global_sample given samples_written using global_index_arr and data_index_arr
 *
 *  uint64_t samples_written - the number of samples written to previous files during this particular user write call
 * 	uint64_t * global_index_arr - an array of global indices into the samples being written.  The global
 * 		index is the total number of sample periods since data taking began, including gaps.
 * 	uint64_t * data_index_arr - an array of len = len(global_index_arr), where the indices are related to which
 * 		sample in the vector being passed in by the user, so that the first value is always 0
 * 	uint_64 index_len - the len of both global_index_arr and data_index_arr.  Must be greater than 0.
 *
 * 	This method tells the global index when some samples have already been written.  If samples_written == 0,
 * 	returns global_index_arr[0].  Else returns global_index_arr[i] + (samples_written-data_index_arr[i]) where
 * 	i is the highest row where data_index_arr[i] <= samples_written
 *
 */
{
	/* local variables */
	int i;
	uint64_t ret_value;

	ret_value = global_index_arr[0] + (samples_written - data_index_arr[0]);
	for (i=1; i<index_len; i++)
	{
		if (samples_written < data_index_arr[i])
			break;
		ret_value = global_index_arr[i] + (samples_written - data_index_arr[i]);
	}
	return(ret_value);
}


int digital_rf_is_little_endian(void)
/* digital_rf_is_little_endian returns 1 if local machine little-endian, 0 if big-endian
 *
 */
{
    volatile uint32_t i=0x01234567;
    // return 0 for big endian, 1 for little endian.
    return (*((uint8_t*)(&i))) == 0x67;
}
