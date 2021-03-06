/*
 * Simple example of writing Digital RF data with C API
 *
 * This simple example writes continuous complex data of short ints
 *
 * $Id: example_rf_write_hdf5.c 812 2015-09-09 19:31:35Z brideout $
 */

#include "digital_rf.h"


int main (int argc, char *argv[])
{

	/* local variables */
	Digital_rf_write_object * data_object = NULL; /* main object created by init */
	uint64_t vector_leading_edge_index = 0; /* index of the sample being written starting at zero with the first sample recorded */
	uint64_t global_start_index; /* start sample (unix time * sample_rate) of first measurement - set below */
	int i, result;

	/* dummy dataset to write */
	short data_short[100][2];

	/* writing parameters */
	double sample_rate = 100.0; /* 100 Hz sample rate - typically MUCH faster */
	uint64_t samples_per_file = 40; /* Number of samples per file - typically MUCH larger */
	uint64_t files_per_directory = 10; /* Each subdirectory will have 40*10=400 samples = 4 seconds */
	int compression_level = 1; /* low level of compression */
	int checksum = 0; /* no checksum */
	int is_complex = 1; /* complex values */
	int num_subchannels = 1; /* only one subchannel */
	int marching_periods = 0; /* no marching periods when writing */
	char uuid[100] = "Fake UUID - use a better one!";
	uint64_t vector_length = 100; /* number of samples written for each call - typically MUCH longer */

	/* init dataset */
	for (i=0; i<100; i++)
	{
		data_short[i][0] = 2*i;
		data_short[i][1] = 3*i;
	}

	/* start recording at global_start_sample */
	global_start_index = (uint64_t)(1394368230 * sample_rate) + 1; /* should represent 2014-03-09 12:30:30  and 10 milliseconds*/


	printf("Writing complex short to multiple files and subdirectores in /tmp/hdf5 channel junk0\n");
	system("rm -rf /tmp/hdf5 ; mkdir /tmp/hdf5 ; mkdir /tmp/hdf5/junk0");

	/* init */
	data_object = digital_rf_create_write_hdf5("/tmp/hdf5/junk0", H5T_NATIVE_INT, samples_per_file, files_per_directory,
			global_start_index, sample_rate, uuid, compression_level, checksum, is_complex, num_subchannels, marching_periods);

	/* write continuous data */
	for (i=0; i<7; i++) /* writing 700 samples, so should create two subdirectories (each holds 400 samples) */
	{
		result = digital_rf_write_hdf5(data_object, vector_leading_edge_index + i*100, data_short, vector_length);
		if (result)
			exit(-1);
	}

	/* close */
	digital_rf_close_write_hdf5(data_object);

	printf("example done - examine /tmp/hdf5 for data\n");
	return(0);
}
