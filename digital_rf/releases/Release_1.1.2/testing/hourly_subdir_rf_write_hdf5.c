/*
 * Test of crossing subdirectories when files_per_subdirectory == 0
 *
 * $Id: hourly_subdir_rf_write_hdf5.c 667 2014-10-15 17:50:42Z brideout $
 */
#include <time.h>
#include <stdio.h>
#include "digital_rf.h"

void digital_rf_randomize_int16(int16_t *data, int len)
{
  int i;
  for(i=0 ; i<len ; i++)
  {
    data[i]=(i%32768)*(i+8192)*(i%13);
  }
}


// length of random number buffer 
#define NUM_SUBCHANNELS 4
#define RANDOM_BLOCK_SIZE 4194304 * NUM_SUBCHANNELS
#define WRITE_BLOCK_SIZE 1000000
// set first time to be March 9, 2014
#define START_TIMESTAMP 1394368230
#define SAMPLE_RATE 1.0E4


int main (int argc, char *argv[])
{
  /* 16-bit integer data */
  int16_t *data_int16;
  uint64_t i, result;
  uint64_t vector_length;
  int n_writes;
  clock_t begin, end;
  double time_spent;
  data_int16 = (int16_t *)malloc(RANDOM_BLOCK_SIZE*sizeof(int16_t));
  vector_length=WRITE_BLOCK_SIZE;
  n_writes = 1000;
  printf("randomize data vector\n");
  digital_rf_randomize_int16(data_int16,RANDOM_BLOCK_SIZE);
  
  /* local variables */
  Digital_rf_write_object *data_object = NULL;
  uint64_t vector_leading_edge_index = 0;
  uint64_t global_start_sample = (uint64_t)(START_TIMESTAMP * SAMPLE_RATE);
  


  printf("Test 0 - simple single write to multiple files, no compress, files_per_subdirectory=0 no checksum - channel 0\n");
  system("rm -rf /tmp/hdf5/junk0 ; mkdir /tmp/hdf5/junk0");
  printf("Start writing\n");
  vector_leading_edge_index=0;
  data_object = digital_rf_create_write_hdf5("/tmp/hdf5/junk0", H5T_NATIVE_SHORT, WRITE_BLOCK_SIZE, 0, global_start_sample, SAMPLE_RATE, "FAKE_UUID_0", 0, 0, 1, NUM_SUBCHANNELS, 1);
  begin = clock();

  if (!data_object)
    exit(-1);
  for(i=0 ; i<n_writes ; i++)
  {
	result = digital_rf_write_hdf5(data_object, vector_leading_edge_index, data_int16, vector_length);
    vector_leading_edge_index+=vector_length;
    
    printf("i is %" PRIu64 "\n", i);

    if (result)
      exit(-1);
  }
  digital_rf_close_write_hdf5(data_object);
  
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("done test %1.2f MB/s\n",((double)n_writes*4*NUM_SUBCHANNELS*vector_length)/time_spent/1e6);

  /*system("rm -rf /tmp/hdf5/junk0");*/
  free(data_int16);
  return(0);
}
