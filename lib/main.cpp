// Added by Nassim.O
#include "AstroAccelerate/headers_mains.h"
#include "AstroAccelerate/device_bin.h"
#include "AstroAccelerate/device_init.h"
#include "AstroAccelerate/device_dedisperse.h"
#include "AstroAccelerate/device_dedispersion_kernel.h"

#include "AstroAccelerate/device_single_pulse_search_kernel.h" //Added by KA
#include "AstroAccelerate/device_single_pulse_search.h" //Added by KA
#include "AstroAccelerate/device_MSD_grid.h" //Added by KA
#include "AstroAccelerate/device_MSD_plane.h" //Added by KA
#include "AstroAccelerate/device_MSD_limited.h" //Added by KA
#include "AstroAccelerate/device_SNR_limited.h" //Added by KA
#include "AstroAccelerate/device_threshold.h" //Added by KA
#include "AstroAccelerate/device_single_FIR.h" //Added by KA

#include "AstroAccelerate/device_load_data.h"
#include "AstroAccelerate/device_corner_turn.h"
#include "AstroAccelerate/device_save_data.h"
#include "AstroAccelerate/host_acceleration.h"
#include "AstroAccelerate/host_allocate_memory.h"
#include "AstroAccelerate/host_analysis.h"
#include "AstroAccelerate/host_periods.h"
#include "AstroAccelerate/host_debug.h"
#include "AstroAccelerate/host_get_file_data.h"
#include "AstroAccelerate/host_get_recorded_data.h"
#include "AstroAccelerate/host_get_user_input.h"
#include "AstroAccelerate/host_help.h"
#include "AstroAccelerate/host_rfi.h"
#include "AstroAccelerate/host_stratagy.h"
#include "AstroAccelerate/host_write_file.h"

#include "AstroAccelerate/host_main_function.h" // Added by Nassim.O

#include "AstroAccelerate/params.h"

int main(int argc, char* argv[])
{		
	// Internal code variables

	// File pointers
	FILE *fp=NULL;

	// Counters and flags
	int i, t, dm_range;
	int range = 0;
	int enable_debug = 0;	
	int enable_analysis = 0;	
	int enable_acceleration = 0;	
	int enable_periodicity = 0;	
	int output_dmt = 0;	
	int *inBin = NULL;
	int *outBin = NULL;
	int *ndms = NULL;
	int maxshift = 0;
	int max_ndms = 0;
	int max_samps = 0;
	int num_tchunks = 0;
	int total_ndms = 0;
	int multi_file = 1;
	
	float max_dm = 0.0f;

	// Memory sizes and pointers
    size_t inputsize = 0;
    size_t outputsize = 0;
	size_t gpu_inputsize = 0;
	size_t gpu_outputsize = 0;
	size_t gpu_memory = 0;

    unsigned short  *input_buffer = NULL;
	float ***output_buffer = NULL;
/***************************************/
	// float **output_sps = NULL;
	// this pointer will be pass to the library, and then
	// to the SPS analysis call
	// analysis will allocate mem and send it back to host


	unsigned short  *d_input = NULL;
	float *d_output = NULL;

	float *dmshifts = NULL;

	float *user_dm_low = NULL;
	float *user_dm_high = NULL;
	float *user_dm_step = NULL;
	float *dm_low = NULL;
	float *dm_high = NULL;
	float *dm_step = NULL;

	// Telescope parameters
	int nchans = 0;
	int nsamp = 0;
	int nbits = 0;
	int nsamples = 0;
	int nifs = 0;
	int **t_processed;

	int nboots = -1;
	int ntrial_bins;
	int navdms = 1;
	int nsearch = 3;
	
	float aggression = 2.5;
	float narrow = 0.001f;
	float wide = 0.1f;

	int	maxshift_original;
	double	tsamp_original;

	long int inc = 0;

	float tstart = 0.0f;
	float tstart_local = 0.0f;
	float tsamp = 0.0f;
	float fch1 = 0.0f;
	float foff = 0.0f;
	
	// Analysis variables
	float power = 2.0f;
	float sigma_cutoff = 6.0f;

	// Timing parameters
	double start_time = omp_get_wtime();

	// Users desired de-dispersion strategy. Pick up user defined values from the CLI.
	get_user_input(&fp, argc, argv, &multi_file, &enable_debug, &enable_analysis, &enable_periodicity, &enable_acceleration, &output_dmt, &nboots, &ntrial_bins, &navdms, &narrow, &wide, &aggression, &nsearch, &inBin, &outBin, &power,	     &sigma_cutoff, &range, &user_dm_low, &user_dm_high, &user_dm_step);
	if(enable_debug == 1) debug(1, start_time, range, outBin, enable_debug, enable_analysis, output_dmt, multi_file, sigma_cutoff, power, max_ndms, user_dm_low, user_dm_high, 
	user_dm_step, dm_low, dm_high, dm_step, ndms, nchans, nsamples, nifs, nbits, tsamp, tstart, fch1, foff, maxshift, max_dm, nsamp, gpu_inputsize, gpu_outputsize, inputsize, outputsize);
	
	// Reads telescope parameters from the header of the input file and then counts the number of samples in the input data file.
	get_file_data(&fp, &nchans, &nsamples, &nsamp, &nifs, &nbits, &tsamp, &tstart, &fch1, &foff);
	if(enable_debug == 1) debug(3, start_time, range, outBin, enable_debug, enable_analysis, output_dmt, multi_file, sigma_cutoff, power, max_ndms, user_dm_low, user_dm_high,
	user_dm_step, dm_low, dm_high, dm_step, ndms, nchans, nsamples, nifs, nbits, tsamp, tstart, fch1, foff, maxshift, max_dm, nsamp, gpu_inputsize, gpu_outputsize, inputsize, outputsize);

	allocate_memory_cpu_input(&fp, gpu_memory, maxshift, num_tchunks, max_ndms, total_ndms, nsamp, nchans, nbits, range, ndms, t_processed, &input_buffer, &output_buffer, &d_input, &d_output, 
                        &gpu_inputsize, &gpu_outputsize, &inputsize, &outputsize);
	if(enable_debug == 1) debug(5, start_time, range, outBin, enable_debug, enable_analysis, output_dmt, multi_file, sigma_cutoff, power, max_ndms, user_dm_low, user_dm_high, 
	user_dm_step, dm_low, dm_high, dm_step, ndms, nchans, nsamples, nifs, nbits, tsamp, tstart, fch1, foff, maxshift, max_dm, nsamp, gpu_inputsize, gpu_outputsize, inputsize, outputsize);
 
	// needs allocate_memory, which needs stratagy, 
	// which needs init_gpu. Split all of them ? 
	// -> need to put get_recorded_data(...) out of this function
	// Store the recorded telescope data contained in the input filterbank file in the allocated memory.
	get_recorded_data(&fp, nsamp, nchans, nbits, &input_buffer, &inputsize);
     	if(enable_debug == 1) debug(7, start_time, range, outBin, enable_debug, enable_analysis, output_dmt, multi_file, sigma_cutoff, power, max_ndms, user_dm_low, user_dm_high, 
	user_dm_step, dm_low, dm_high, dm_step, ndms, nchans, nsamples, nifs, nbits, tsamp, tstart, fch1, foff, maxshift, max_dm, nsamp, gpu_inputsize, gpu_outputsize, inputsize, outputsize);

	// call to library 
	main_function(argc,
				  argv,
				  // File pointers
				  fp, 	
				  // Counters and flags
				  i, 	
				  t,
				  dm_range,
				  range,
				  enable_debug,
				  enable_analysis,
				  enable_acceleration,
				  enable_periodicity,
				  output_dmt, 
				  inBin, 
				  outBin, 
				  ndms, 
				  maxshift,
				  max_ndms,
				  max_samps,
				  num_tchunks,
				  total_ndms,
				  multi_file,
				  max_dm,
				  // Memory sizes and pointers
				  inputsize, 
				  outputsize,
				  gpu_inputsize,
				  gpu_outputsize,
				  gpu_memory,
				  input_buffer,
				  output_buffer, // pointer *output_buffer
				  d_input,
				  d_output,
				  dmshifts,
				  user_dm_low,
				  user_dm_high,
				  user_dm_step,
				  dm_low,
				  dm_high,
				  dm_step,
				  nchans, // Telescope parameters
				  nsamp,
				  nbits,
				  nsamples,
				  nifs,
				  t_processed,
				  nboots,
				  ntrial_bins,
				  navdms,
				  nsearch,
				  aggression,
				  narrow,
				  wide,
				  maxshift_original,
				  tsamp_original,
				  inc,
				  tstart,
				  tstart_local,
				  tsamp,
				  fch1,
				  foff,
				  // Analysis variables
				  power, 
				  sigma_cutoff,
				  start_time
				  );

	// exports should be done here, not in the library

	return 0;

}