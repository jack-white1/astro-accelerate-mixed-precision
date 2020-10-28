#include  <stdio.h>
#include <stdlib.h>
#include <cmath>

#include "aa_params.hpp"
#include "aa_host_rfi.hpp"
#include "aa_host_export.hpp"

namespace astroaccelerate {

  void rfi(int nsamp, int nchans, std::vector<unsigned short> &input_buffer) {
    int 	file_reducer 	= 1;
    float 	sigma_cut 	= 2.0f;
	
    float 	*stage = (float*)malloc((size_t)nsamp*(size_t)nchans*sizeof(float));

    for(int c = 0; c < nchans; c++) {
      for(int t = 0; t < (nsamp); t++) {
	stage[c * (size_t)nsamp + t] = (float)input_buffer[c  + (size_t)nchans * t];
      }
    }

    // ~~~ RFI Correct ~~~ //	
    double orig_mean = 0.0;
    double orig_var=0.0;

    // Find the mean and SD of the input data (we'll use this to rescale the data at the end of the process.

    for(int c = 0; c < nchans; c++) {
      for(int t = 0; t < (nsamp); t++) orig_mean+=stage[c * (size_t)nsamp + t];
    }
    orig_mean/=(nsamp*nchans);

    for(int c = 0; c < nchans; c++) {
      for(int t = 0; t < (nsamp); t++) {
	orig_var+=(stage[c * (size_t)nsamp + t]-orig_mean)*(stage[c * (size_t)nsamp + t]-orig_mean);
      }
    }
    orig_var/=(nsamp*nchans);
    orig_var=sqrt(orig_var);

    //printf("\n%lf\t%lf", orig_mean, orig_var);
	
    // Random Vectors
	
    float *random_chan_one = (float*)malloc(nsamp*sizeof(float));
    float *random_chan_two = (float*)malloc(nsamp*sizeof(float));

    for(int t = 0; t < nsamp; t++) {
	
      float x1, x2, w, y1, y2;

      do {
	x1 = 2.0 * ((float)rand()/(float)RAND_MAX) -1.0;
	x2 = 2.0 * ((float)rand()/(float)RAND_MAX) -1.0;
	w = x1*x1 + x2*x2;
      } while( w >= 1.0 );

      w = sqrt((-2.0 * log (w))/ w);
      y1 = x1*w;
      y2 = x2*w;


      random_chan_one[t] = y1;
      random_chan_two[t] = y2;
    }

    float *random_spectra_one = (float*)malloc(nchans*sizeof(float));
    float *random_spectra_two = (float*)malloc(nchans*sizeof(float));

    for(int c = 0; c < nchans; c++) {
	
      float x1, x2, w, y1, y2;

      do {
	x1 = 2.0 * ((float)rand()/(float)RAND_MAX) -1.0;
	x2 = 2.0 * ((float)rand()/(float)RAND_MAX) -1.0;
	w = x1*x1 + x2*x2;
      } while( w >= 1.0 );

      w = sqrt((-2.0 * log (w))/ w);
      y1 = x1*w;
      y2 = x2*w;


      random_spectra_one[c] = y1;
      random_spectra_two[c] = y2;
    }
	
    // Allocate working arrays

    int *chan_mask = (int*)malloc(nchans*sizeof(int));
    for(int c = 0; c < nchans; c++) chan_mask[c]=1;

    int *spectra_mask = (int*)malloc(nsamp*sizeof(int));
    for(int t = 0; t < nsamp; t++) spectra_mask[t]=1;
 
    double *chan_mean = (double*)malloc(nchans*sizeof(double));
    for(int c = 0; c < nchans; c++) chan_mean[c] = 0.0;

    double *chan_var = (double*)malloc(nsamp*sizeof(double));
    for(int c = 0; c < nchans; c++) chan_var[c] = 0.0;

    double *spectra_mean = (double*)malloc(nsamp*sizeof(double));
    for(int t = 0; t < nsamp; t++) spectra_mean[t] = 0.0;

    double *spectra_var = (double*)malloc(nsamp*sizeof(double));
    for(int t = 0; t < nsamp; t++) spectra_var[t] = 0.0;

    // Find the BLN and try to flatten the input data per channel (remove non-stationary component).

    for(int c = 0; c < nchans; c++) {
  
      int counter = 0;

      for(int t = 0; t < nsamp; t++) spectra_mask[t]=1;

      int finish = 0;
      int rounds = 1;

      double old_mean = 0.0;
      double old_var = 0.0;

      while(finish == 0) {
		
	counter = 0;
	chan_mean[c] = 0.0;
	for(int t = 0; t < (nsamp); t++) {
	  if(spectra_mask[t] == 1) {
	    chan_mean[c] += stage[c * (size_t)nsamp + t];
	    counter++;
	  }
	}
	if(counter == 0) {
	  printf("\nCounter zero, Channel %d", c);
	  chan_mask[c] = 0;				
	  finish = 1;
	  break;
	}
	chan_mean[c]/=(counter);

	counter = 0;
	chan_var[c] = 0.0;
	for(int t = 0; t < (nsamp); t++) {
	  if(spectra_mask[t] == 1) {
	    chan_var[c] += (stage[c * (size_t)nsamp + t]-chan_mean[c])*(stage[c * (size_t)nsamp + t]-chan_mean[c]);
	    counter++;
	  }
	}
	chan_var[c] /= (counter);
	chan_var[c] = sqrt(chan_var[c]);

	if((chan_var[c])*1000000.0 < 0.1) {
	  printf("\nVarience zero, Channel %d %d %lf %.16lf", c, rounds, chan_mean[c], chan_var[c] );
	  chan_mask[c] = 0;
	  finish = 1;
	  break;
	}

	for(int t = 0; t < (nsamp); t++) {
	  if(((stage[c * (size_t)nsamp + t]-chan_mean[c])/chan_var[c]) > sigma_cut || ((stage[c * (size_t)nsamp + t]-chan_mean[c])/chan_var[c]) < -sigma_cut) {
	    spectra_mask[t]=0;
	  } else {
	    spectra_mask[t]=1;
	  }
	}

	if(fabs(chan_mean[c] - old_mean) < 0.001 && fabs(chan_var[c] - old_var) < 0.0001 && rounds > 1) {
	  //printf("\n%d\t%d\t%.16lf\t%.16lf\t%.16lf\t%.16lf", c, rounds, (chan_mean[c]-old_mean), (chan_var[c]-old_var), chan_mean[c], chan_var[c]); 
	  finish = 1;
	}

	old_mean = chan_mean[c];
	old_var = chan_var[c];
	rounds++;
      }

      if(chan_mask[c] != 0) {
	for(int t = 0; t < (nsamp); t++) {
	  stage[c * (size_t)nsamp + t]=(stage[c * (size_t)nsamp + t]-(float)chan_mean[c])/(float)chan_var[c];
	}
      } else {
	int perm_one = (int)(((float)rand()/(float)RAND_MAX)*nsamp);
	for(int t = 0; t < nsamp; t++) {
	  stage[c * (size_t)nsamp + t] = random_chan_one[(t+perm_one)%nsamp];
	}
	chan_mean[c] = 0.0;
	chan_var[c]  = 1.0;
	chan_mask[c] = 1;
      }
    }

    // Find the BLN and try to flatten the input data per spectra (remove non-stationary component).

    for(int t = 0; t < (nsamp); t++) {

      int counter = 0;

      for(int c = 0; c < nchans; c++) chan_mask[c]=1;

      int finish = 0;
      int rounds = 1;

      double old_mean = 0.0;
      double old_var = 0.0;

      while(finish == 0) {
		
	counter = 0;
	spectra_mean[t] = 0.0;
	for(int c = 0; c < nchans; c++) {
	  if(chan_mask[c] == 1) {
	    spectra_mean[t]+=stage[c * (size_t)nsamp + t];
	    counter++;
	  }
	}
	if(counter == 0) {
	  printf("\nCounter zero, Spectra %d", t);
	  spectra_mask[t] = 0;
	  finish = 1;
	  break;
	}
	spectra_mean[t] /= (counter);

	counter = 0;
	spectra_var[t] = 0.0;
	for(int c = 0; c < nchans; c++) {
	  if(chan_mask[c] == 1) {
	    spectra_var[t] += (stage[c * (size_t)nsamp + t]-spectra_mean[t])*(stage[c * (size_t)nsamp + t]-spectra_mean[t]);
	    counter++;
	  }
	}
	spectra_var[t] /= (counter);
	spectra_var[t] = sqrt(spectra_var[t]);

	if((spectra_var[t])*1000000.0 < 0.1) {
	  printf("\nVarience zero, Spectra %d %d %lf %.16lf", t, rounds, spectra_mean[t], spectra_var[t] );
	  spectra_mask[t] = 0;
	  finish = 1;
	  break;
	}

	if(spectra_mask[t] != 0) {
	  for(int c = 0; c < nchans; c++) {
	    if(((stage[c * (size_t)nsamp + t]-spectra_mean[t])/spectra_var[t]) > sigma_cut || ((stage[c * (size_t)nsamp + t]-spectra_mean[t])/spectra_var[t]) < -sigma_cut) {
	      chan_mask[c]=0;
	    } else {
	      chan_mask[c]=1;
	    }
	  }
	}
			
	if(fabs(spectra_mean[t] - old_mean) < 0.001 && fabs(spectra_var[t] - old_var) < 0.0001 && rounds > 1) {
	  //printf("\n%d\t%d\t%.16lf\t%.16lf\t%.16lf\t%.16lf", t, rounds, (spectra_mean[t] - old_mean), (spectra_var[t] - old_var), spectra_mean[t], spectra_var[t]); 
	  finish = 1;
	}

	old_mean = spectra_mean[t];
	old_var = spectra_var[t];
	rounds++;
      }
		
      if(spectra_mask[t] != 0) {
	for(int c = 0; c < nchans; c++) {
	  stage[c * (size_t)nsamp + t]=(stage[c * (size_t)nsamp + t]-(float)spectra_mean[t])/(float)spectra_var[t];
	}
      } else {
	int perm_one = (int)(((float)rand()/(float)RAND_MAX)*nchans);
	for(int c = 0; c < nchans; c++) {
	  stage[c * (size_t)nsamp + t] = random_spectra_one[(c+perm_one)%nchans];
	}
	spectra_mean[t] = 0.0;
	spectra_var[t]  = 1.0;
	spectra_mask[t] = 1;
      }
    }

    double mean_rescale = 0.0;
    double var_rescale  = 0.0;

    // Find the mean and SD of the mean and SD...
    int finish = 0;
    int rounds = 1;
    int counter = 0;

    double mean_of_mean = 0.0;
    double var_of_mean  = 0.0;
    double mean_of_var  = 0.0;
    double var_of_var   = 0.0;

    double old_mean_of_mean = 0.0;
    double old_var_of_mean  = 0.0;
    double old_mean_of_var  = 0.0;
    double old_var_of_var   = 0.0;

    for(int c = 0; c < nchans; c++) chan_mask[c]=1;

    while(finish == 0) {

      mean_of_mean = 0.0;
      counter = 0;
      for(int c = 0; c < nchans; c++) {
	if(chan_mask[c] == 1) {			
	  mean_of_mean+=chan_mean[c];
	  counter++;
	}
      }
      mean_of_mean/=counter;

      var_of_mean = 0.0;
      counter = 0;
      for(int c = 0; c < nchans; c++) {
	if(chan_mask[c] == 1) {			
	  var_of_mean+=(chan_mean[c] - mean_of_mean)*(chan_mean[c] - mean_of_mean);
	  counter++;
	}
      }
      var_of_mean/=(counter);
      var_of_mean=sqrt(var_of_mean);

      mean_of_var = 0.0;
      counter = 0;
      for(int c = 0; c < nchans; c++) {
	if(chan_mask[c] == 1) {			
	  mean_of_var+=chan_var[c];
	  counter++;
	}
      }
      mean_of_var/=counter;

      var_of_var = 0.0;
      counter = 0;
      for(int c = 0; c < nchans; c++) {
	if(chan_mask[c] == 1) {			
	  var_of_var+=(chan_var[c] - mean_of_var)*(chan_var[c] - mean_of_var);
	  counter++;
	}
      }
      var_of_var/=(counter);
      var_of_var=sqrt(var_of_var);

      for(int c = 0; c < nchans; c++) if(fabs(chan_mean[c] - mean_of_mean)/var_of_mean > sigma_cut || fabs(chan_var[c] - mean_of_var)/var_of_var > sigma_cut) chan_mask[c] = 0;

      if(fabs(mean_of_mean - old_mean_of_mean)   < 0.001 &&
	 fabs(var_of_mean  - old_var_of_mean )   < 0.001 &&
	 fabs(mean_of_var  - old_mean_of_var )   < 0.001 &&
	 fabs(var_of_var   - old_var_of_var  )   < 0.001)  {
			
	finish = 1;

      }
		
      old_mean_of_mean = mean_of_mean;
      old_var_of_mean  = var_of_mean;
      old_mean_of_var  = mean_of_var;
      old_var_of_var   = var_of_var;
      rounds++;
    }
	
    printf("\n0 %lf %lf", mean_of_mean, var_of_mean);
    printf("\n0 %lf %lf", mean_of_var,  var_of_var);

    mean_rescale = mean_of_mean;
    var_rescale  = mean_of_var;
	
    float clipping_constant = 0.0;
    for(int c = 0; c < nchans; c++) clipping_constant += chan_mask[c];
    clipping_constant = (nchans - clipping_constant)/nchans;
    clipping_constant = sqrt(-2.0 * log(clipping_constant * 2.506628275));

    // Perform channel replacement
    for(int c = 0; c < nchans; c++) {
      if(fabs((chan_mean[c]-mean_of_mean)/var_of_mean) > clipping_constant && fabs((chan_var[c]-mean_of_var)/var_of_var) > clipping_constant) {
	//printf("\nReplacing Channel %d %lf %lf", c, chan_mean[c], chan_var[c]);
	int perm_one = (int)(((float)rand()/(float)RAND_MAX)*nsamp);
	for(int t = 0; t < (nsamp); t++) {
	  stage[(c) * (size_t)nsamp + t] = random_chan_two[(t+perm_one)%nsamp];
	}
      }
    }
	
    finish = 0;
    rounds = 1;
    counter = 0;

    mean_of_mean = 0.0;
    var_of_mean  = 0.0;
    mean_of_var  = 0.0;
    var_of_var   = 0.0;

    old_mean_of_mean = 0.0;
    old_var_of_mean  = 0.0;
    old_mean_of_var  = 0.0;
    old_var_of_var   = 0.0;

    for(int t = 0; t < (nsamp); t++) spectra_mask[t] = 1;	

    while(finish == 0) {

      mean_of_mean = 0.0;
      counter = 0;
      for(int t = 0; t < (nsamp); t++) {
	if(spectra_mask[t] == 1) {			
	  mean_of_mean+=spectra_mean[t];
	  counter++;
	}
      }
      mean_of_mean/=counter;

      var_of_mean = 0.0;
      counter = 0;
      for(int t = 0; t < (nsamp); t++) {
	if(spectra_mask[t] == 1) {			
	  var_of_mean+=(spectra_mean[t] - mean_of_mean)*(spectra_mean[t] - mean_of_mean);
	  counter++;
	}
      }
      var_of_mean/=(counter);
      var_of_mean=sqrt(var_of_mean);

      mean_of_var = 0.0;
      counter = 0;
      for(int t = 0; t < (nsamp); t++) {
	if(spectra_mask[t] == 1) {			
	  mean_of_var+=spectra_var[t];
	  counter++;
	}
      }
      mean_of_var/=counter;

      var_of_var = 0.0;
      counter = 0;
      for(int t = 0; t < (nsamp); t++) {
	if(spectra_mask[t] == 1) {			
	  var_of_var+=(spectra_var[t] - mean_of_var)*(spectra_var[t] - mean_of_var);
	  counter++;
	}
      }
      var_of_var/=(counter);
      var_of_var=sqrt(var_of_var);

      for(int t = 0; t < (nsamp); t++) if(fabs(spectra_mean[t] - mean_of_mean)/var_of_mean > sigma_cut || fabs(spectra_var[t] - mean_of_var)/var_of_var > sigma_cut) spectra_mask[t] = 0;

      if(fabs(mean_of_mean - old_mean_of_mean)   < 0.001 &&
	 fabs(var_of_mean  - old_var_of_mean )   < 0.001 &&
	 fabs(mean_of_var  - old_mean_of_var )   < 0.001 &&
	 fabs(var_of_var   - old_var_of_var  )   < 0.001)  {
			
	finish = 1;

      }
		
      old_mean_of_mean = mean_of_mean;
      old_var_of_mean  = var_of_mean;
      old_mean_of_var  = mean_of_var;
      old_var_of_var   = var_of_var;
      rounds++;
    }

    printf("\n0 %lf %lf", mean_of_mean, var_of_mean);
    printf("\n0 %lf %lf", mean_of_var,  var_of_var);

    clipping_constant = 0.0;
    for(int t = 0; t < nsamp; t++) clipping_constant += spectra_mask[t];
    clipping_constant = (nsamp - clipping_constant)/nsamp;
    clipping_constant = sqrt(-2.0 * log(clipping_constant * 2.506628275));

    // Perform spectral replacement
    for(int t = 0; t < (nsamp); t++){
      if(fabs((spectra_mean[t]-mean_of_mean)/var_of_mean) > clipping_constant && fabs((spectra_var[t]-mean_of_var)/var_of_var) > clipping_constant) {
	//printf("\nReplacing Spectral %d %lf %lf", t, spectra_mean[t], spectra_var[t]);		       
	int perm_one = (int)(((float)rand()/(float)RAND_MAX)*nchans);
	for(int c = 0; c < nchans; c++) {
	  stage[(c) * (size_t)nsamp + t] = random_spectra_two[(c+perm_one)%nchans];
	}
      }
    }


    for(int c = 0; c < nchans; c++) {
      for(int t = 0; t < (nsamp); t++) {
	//(*input_buffer)[c  + (size_t)nchans * t] = (unsigned char) ((stage[c * (size_t)nsamp + t]*orig_var)+orig_mean);
	input_buffer[c  + (size_t)nchans * t] = (unsigned char) ((stage[c * (size_t)nsamp + t]*var_rescale)+mean_rescale);
      }
    }

    FILE *fp_mask = fopen ("masked_chans.txt", "w+");
    for(int c = 0; c < nchans; c++) {
      for(int t = 0; t < (nsamp)/file_reducer; t++) {
	//fprintf(fp_mask, "%d ", (unsigned char)((stage[c * (size_t)nsamp + t]*orig_var)+orig_mean));
	fprintf(fp_mask, "%d ", (unsigned char)((stage[c * (size_t)nsamp + t]*var_rescale)+mean_rescale));
      }
      fprintf(fp_mask, "\n");
    }
    fclose(fp_mask);

    printf("\n%lf %lf", mean_rescale/orig_mean, var_rescale/orig_var);


    free(chan_mask);
    free(spectra_mask);
    free(chan_mean);
    free(chan_var);
    free(spectra_mean);
    free(spectra_var);
    free(stage);
  }


	//---------------------------------------------------------------------------------
	//-------> Kahan MSD
	void d_kahan_summation(float *signal, size_t nDMs, size_t nTimesamples, size_t offset, double *result, double *error, bool outlier_rejection, double old_mean, double old_stdev, double sigma){
		double sum;
		double sum_error;
		double a,b;
		size_t nElements = 0;
		
		double low  = old_mean - sigma*old_stdev;
		double high = old_mean + sigma*old_stdev;
		
		sum=0;
		sum_error=0;
		for(size_t d=0;d<nDMs; d++){
			for(size_t s=0; s<(nTimesamples-offset); s++){
				double sample = signal[(size_t) (d*nTimesamples + s)];
				if(outlier_rejection && (sample<low || sample>high)){
				}
				else{
					a = sample - sum_error;
					b = sum + a;
					sum_error = (b - sum);
					sum_error = sum_error - a;
					sum = b;
					nElements++;
				}
			}
		}
		*result = sum/nElements;
		*error = sum_error;
	}

	void d_kahan_sd(float *signal, size_t nDMs, size_t nTimesamples, size_t offset, double mean, double *result, double *error, bool outlier_rejection, double old_mean, double old_stdev, double sigma){
		double sum;
		double sum_error;
		double a,b,dtemp;
		size_t nElements = 0;
		
		double low  = old_mean - sigma*old_stdev;
		double high = old_mean + sigma*old_stdev;
		
		sum=0;
		sum_error=0;
		for(size_t d=0;d<nDMs; d++){
			for(size_t s=0; s<(nTimesamples-offset); s++){
				double sample = signal[(size_t) (d*nTimesamples + s)];
				if(outlier_rejection && (sample<low || sample>high)){
				}
				else{
					dtemp=(sample - sum_error - mean);
					a=dtemp*dtemp;
					b=sum+a;
					sum_error=(b-sum);
					sum_error=sum_error-a;
					sum=b;
					nElements++;
				}
			}
		}
		*result=sqrt(sum/nElements);
		*error=sum_error;
	}


	void MSD_Kahan(float *h_input, size_t nDMs, size_t nTimesamples, size_t offset, double *mean, double *sd, bool outlier_rejection, double sigma){
		double error, signal_mean, signal_sd;
		double old_mean, old_stdev;
		size_t nElements=nDMs*(nTimesamples-offset);
		
		d_kahan_summation(h_input, nDMs, nTimesamples, offset, &signal_mean, &error, false, 0, 1.0, 1.0);
		d_kahan_sd(h_input, nDMs, nTimesamples, offset, signal_mean, &signal_sd, &error, false, 0, 1.0, 1.0);
		
		old_mean = signal_mean; old_stdev = signal_sd;
		//printf("    Before outlier rejection: %e - %e\n", old_mean, old_stdev);
		if(outlier_rejection){
			for(int f=0; f<5; f++){
				d_kahan_summation(h_input, nDMs, nTimesamples, offset, &signal_mean, &error, true, old_mean, old_stdev, sigma);
				d_kahan_sd(h_input, nDMs, nTimesamples, offset, signal_mean, &signal_sd, &error, true, old_mean, old_stdev, sigma);
				
				old_mean = signal_mean; old_stdev = signal_sd;
				//printf("      Iteration %d of outlier rejection: %e - %e\n", f, old_mean, old_stdev);
			}
		}

		*mean=signal_mean;
		*sd=signal_sd;
	}
	//-------> Kahan MSD
	//---------------------------------------------------------------------------------

	void export_input_data(size_t nsamp, size_t nchans, std::vector<unsigned short> &input_buffer, const char *filename){
		float *temp;
		temp = new float[nsamp*nchans];
		
		for(size_t t=0; t<nsamp; t++){
			// transfer data into temporary array
			for(size_t c=0; c<nchans; c++){
				temp[c + nchans*t] = (float) input_buffer[c  + nchans*t];
			}
		}
		
		Export_data_as_list(temp, nchans, 1.0, 0, nsamp/10, 1.0, 0, filename, 5000);

		delete[] temp;
	}

	void input_data_renormalization(size_t nsamp, size_t nchans, std::vector<unsigned short> &input_buffer, bool enable_outlier_rejection, float sigma){
		int per;
		
		export_input_data(nsamp, nchans, input_buffer, "before_data");
		
		float *temp_spectrum;
		temp_spectrum = new float[nchans];
		per = 0;
		for(size_t t=0; t<nsamp; t++){
			// transfer data into temporary array
			for(size_t c=0; c<nchans; c++){
				temp_spectrum[c] = (float) input_buffer[c  + nchans*t];
			}
			// calculation of MSD
			double mean, stdev;
			MSD_Kahan(temp_spectrum, 1, nchans, 0, &mean, &stdev, enable_outlier_rejection, sigma);
			// renormalization
			mean = mean - 127.5;
			for(size_t c=0; c<nchans; c++){
				unsigned short value;
				float result = (float)input[c  + nchans*t] - mean;
				if(result<0) value = 0;
				else if(result>255) value = 255;
				else value = (unsigned short) result;
				input[c  + nchans*t] = value;
			}
			if(t%(nsamp/100)==0) {
				if(per==0 || per==25 || per==50 || per==75) printf("%d%%", per);
				else printf(".");
				fflush(stdout);
				per++;
			}
		}
		printf("100%%\n");
		delete [] temp_spectrum;
		
		per = 0;
		float *temp_time;
		temp_time = new float[nsamp];
		for(size_t c=0; c<nchans; c++){
			// transfer data into temporary array
			for(size_t t=0; t<nsamp; t++){
				temp_time[t] = (float) input_buffer[c  + nchans*t];
			}
			// calculation of MSD
			double mean, stdev;
			MSD_Kahan(temp_time, 1, nchans, 0, &mean, &stdev, enable_outlier_rejection, sigma);
			// renormalization
			mean = mean - 127.5;
			for(size_t t=0; t<nsamp; t++){
				unsigned short value;
				float result = (float)input[c  + nchans*t] - mean;
				if(result<0) value = 0;
				else if(result>255) value = 255;
				else value = (unsigned short) result;
				input[c  + nchans*t] = value;
			}
			if(c%(nchans/100)==0) {
				if(per==0 || per==25 || per==50 || per==75) printf("%d%%", per);
				else printf(".");
				fflush(stdout);
				per++;
			}
		}
		printf("100%%\n");		
		delete [] temp_time;
		
		export_input_data(nsamp, nchans, input_buffer, "after_data");
	}
	
	void CPU_corner_turn(float *data, float *CT_data, size_t primary_size, size_t secondary_size){
		for(int s=0; s<secondary_size; s++){
			for(int p=0; p<primary_size; p++){
				CT_data[p*secondary_size + s]=data[s*primary_size + p];
			}
		}
	}
	
	void dedisperse_DM_renormalization(size_t nDMs, size_t nTimesamples, float *input_buffer, bool enable_outlier_rejection, float sigma){
		int per;
		
		float *CTtemp;
		CTtemp = new float[nTimesamples*nDMs];
		printf("Transpose..."); fflush(stdout);
		CPU_corner_turn(input_buffer, CTtemp, nTimesamples, nDMs);
		printf(" done.\n");
		
		per = 0;
		float *temp_time;
		temp_time = new float[nDMs];
		for(size_t t=0; t<nTimesamples; t++){
			// transfer data into temporary array
			for(size_t d=0; d<nDMs; d++){
				temp_time[d] = CTtemp[d + t*nDMs];
			}
			// calculation of MSD
			double mean, stdev;
			MSD_Kahan(temp_time, 1, nDMs, 0, &mean, &stdev, enable_outlier_rejection, sigma);
			// renormalization
			mean = mean - 127.5;
			for(size_t d=0; d<nDMs; d++){
				CTtemp[d + t*nDMs] = CTtemp[d + t*nDMs] - mean;
			}
			if(t%(nTimesamples/100)==0) {
				if(per==0 || per==25 || per==50 || per==75) printf("%d%%", per);
				else printf(".");
				fflush(stdout);
				per++;
			}
		}
		printf("100%%\n");	

		printf("Transpose..."); fflush(stdout);
		CPU_corner_turn(CTtemp, input_buffer, nDMs, nTimesamples);
		printf(" done.\n");
		
		delete [] temp_time;
		delete [] CTtemp;
	}

	void dedisperse_time_renormalization(size_t nDMs, size_t nTimesamples, float *input_buffer, bool enable_outlier_rejection, float sigma){
		int per;
		
		per = 0;
		float *temp_time;
		temp_time = new float[nTimesamples];
		for(size_t d=0; d<nDMs; d++){
			// transfer data into temporary array
			for(size_t t=0; t<nTimesamples; t++){
				temp_time[t] = input_buffer[d*nTimesamples + t];
			}
			// calculation of MSD
			double mean, stdev;
			MSD_Kahan(temp_time, 1, nTimesamples, 0, &mean, &stdev, enable_outlier_rejection, sigma);
			// renormalization
			mean = mean - 127.5;
			for(size_t t=0; t<nTimesamples; t++){
				input_buffer[d*nTimesamples + t] = input_buffer[d*nTimesamples + t] - mean;
			}
			if(d%(nDMs/100)==0) {
				if(per==0 || per==25 || per==50 || per==75) printf("%d%%", per);
				else printf(".");
				fflush(stdout);
				per++;
			}
		}
		printf("100%%\n");		
		delete [] temp_time;
	}

} //namespace astroaccelerate
