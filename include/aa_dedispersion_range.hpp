#ifndef ASTRO_ACCELERATE_AA_DEDISPERSION_RANGE_HPP
#define ASTRO_ACCELERATE_AA_DEDISPERSION_RANGE_HPP

#include <stdio.h>

namespace astroaccelerate {

class aa_dedispersion_range {
private:
	double c_dm_low;         /** Lowest value of the dispersion measure in this dedispersion range [pc*cm-3]*/
	double c_dm_high;        /** Highest value of the dispersion measure in this dedispersion range [pc*cm-3]*/
	double c_dm_step;        /** Step in dispersion measure [pc*cm-3]*/
	double c_sampling_time;  /** Sampling time before binning [s]*/
	int c_inbin;             /** Binning factor*/
	unsigned long int c_nTimesamples; /** Number of time-samples before binning */
	int c_nDMs;              /** Number of DM-trials in this dedispersion range*/

public:
	/** \brief Constructor for Dedispersion_Range. */
	aa_dedispersion_range() {
		c_dm_step = 0;
		c_dm_low = 0;
		c_dm_high = 0;
		c_inbin = 0;
		c_nTimesamples = 0;
		c_nDMs = 0;
		c_sampling_time = 0;
	}

	/** \brief Constructor for Dedispersion_Range. */
	aa_dedispersion_range(
		double dm_low, 
		double dm_high, 
		double dm_step, 
		int inBin, 
		unsigned long int range_timesamples_before_binning, 
		int ndms, 
		double sampling_time_before_binning
	) {
		c_dm_step = dm_step;
		c_dm_low = dm_low;
		c_dm_high = dm_high;
		c_inBin = inBin;
		c_nTimesamples = range_timesamples_before_binning;
		c_nDMs = ndms;
		c_sampling_time = sampling_time_before_binning;
	}

	/** \brief Method to assign or change values of an instance after construction. */
	void Assign(
		double dm_low, 
		double dm_high, 
		double dm_step, 
		int inBin, 
		unsigned long int range_timesamples_before_binning, 
		int ndms, 
		double sampling_time_before_binning
	) {
		c_dm_step = dm_step;
		c_dm_low = dm_low;
		c_dm_high = dm_high;
		c_inBin = inBin;
		c_nTimesamples = range_timesamples_before_binning;
		c_nDMs = ndms;
		c_sampling_time = sampling_time_before_binning;
	}

	double dm_low() return(c_dm_low);
	double dm_high() return(c_dm_high);
	double dm_step() return(c_dm_step);
	double sampling_time() return(c_sampling_time);
	int inBin() return(c_inBin);
	int nDMs() return(c_nDMs);
	unsigned long int nTimesamples() return(c_nTimesamples);
};

} // namespace astroaccelerate
#endif // ASTRO_ACCELERATE_AA_DEDISPERSION_RANGE_HPP