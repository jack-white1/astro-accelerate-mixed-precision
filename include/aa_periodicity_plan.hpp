#ifndef ASTRO_ACCELERATE_AA_PERIODICITY_PLAN_HPP
#define ASTRO_ACCELERATE_AA_PERIODICITY_PLAN_HPP

#include <stdio.h>
#include "aa_dedispersion_range.hpp"

namespace astroaccelerate {

/**
 * \class aa_periodicity_plan aa_periodicity_plan.hpp "include/aa_periodicity_plan.hpp"
 * \brief Class to set a periodicity plan.
 * \details A periodicity plan is required to create a periodicity strategy.
 * \author AstroAccelerate team.
 */

class aa_periodicity_plan {
public:
	float sigma_cutoff; /** User selected value of sigma cutoff. Any event with SNR (sigma) below this value will not be selected as candidate */
	float sigma_outlier_rejection_threshold; /** User selected sigma for outlier rejection. Any value with sigma greater than this will be rejected from calculation of mean and standard deviation. */
	int   nHarmonics; /** User selected number of harmonics performed by the periodicity search. */
	int   candidate_algorithm; /** User selected flag to select reduction algorithm to select candidates. */
	bool  enable_outlier_rejection; /** User selected flag to enable or disable outlier rejection when calculating mean and standard deviation. */
	bool  enable_interpolation; /** Enable or disable interpolation for periodicity search */
	bool  enable_spectrum_whitening; /** Enable or disable spectrum whitening (removal of the red noise) for periodicity search */
	bool  pad_to_nearest_higher_pow2; /** Whether the periodicity will pad data to the nearest power of two for the Fourier transform. True by default */

	/** \brief Trivial constructor for aa_periodicity_plan. */
	aa_periodicity_plan() {
		sigma_cutoff = 0;
		sigma_outlier_rejection_threshold = 0;
		nHarmonics = 0;
		candidate_algorithm = 0;
		enable_outlier_rejection = false;
		enable_interpolation = false;
		enable_spectrum_whitening = false;
		pad_to_nearest_higher_pow2 = true;
	}

	/*
	// This constructor is here in case periodicity search would be independent not allowed at the moment
	aa_periodicity_plan(
		std::vector<aa_dedispersion_range> &t_ddtr_ranges,
		const float &t_sigma_cutoff,
		const float &t_sigma_outlier_rejection_threshold,
		const int &t_nHarmonics,
		const int &t_candidate_algorithm,
		const bool &t_enable_outlier_rejection,
		const bool &t_enable_interpolation,
		const bool &t_enable_spectrum_whitening,
		const bool &t_pad_to_nearest_higher_pow2
	) {
		sigma_cutoff = t_sigma_cutoff;
		sigma_outlier_rejection_threshold = t_sigma_outlier_rejection_threshold;
		nHarmonics = t_nHarmonics;
		candidate_algorithm = t_candidate_algorithm;
		enable_outlier_rejection = t_enable_outlier_rejection;
		enable_interpolation = t_enable_interpolation;
		enable_spectrum_whitening = t_enable_spectrum_whitening;
		pad_to_nearest_higher_pow2 = t_pad_to_nearest_higher_pow2;
		
		ddtr_ranges = t_ddtr_ranges;
	}
	*/
	
	/** \brief Constructor for aa_periodicity_plan that sets all member data on construction. */
	aa_periodicity_plan(
		aa_ddtr_strategy &ddtr_strategy;
		const float &t_sigma_cutoff,
		const float &t_sigma_outlier_rejection_threshold,
		const int &t_nHarmonics,
		const int &t_candidate_algorithm,
		const bool &t_enable_outlier_rejection,
		const bool &t_enable_interpolation,
		const bool &t_enable_spectrum_whitening,
		const bool &t_pad_to_nearest_higher_pow2
	) {
		sigma_cutoff = t_sigma_cutoff;
		sigma_outlier_rejection_threshold = t_sigma_outlier_rejection_threshold;
		nHarmonics = t_nHarmonics;
		candidate_algorithm = t_candidate_algorithm;
		enable_outlier_rejection = t_enable_outlier_rejection;
		enable_interpolation = t_enable_interpolation;
		enable_spectrum_whitening = t_enable_spectrum_whitening;
		pad_to_nearest_higher_pow2 = t_pad_to_nearest_higher_pow2;

		unsigned long int nTimesamples = 0;
		int nTimechunks = ddtr_strategy.m_num_tchunks;
		for(int f = 0; f < nTimechunks; f++){
			nTimesamples = nTimesamples + ddtr_strategy.m_t_processed[0][f];
		}
		
		double sampling_time = ddtr_strategy.m_metadata.tsamp();
		
		int nRanges = ddtr_strategy.m_ndms.size();
		for (int f = 0; f < nRanges; f++) {
			double dm_low = (double) ddtr_strategy[f].str_dm.low;
			double dm_high = (double) ddtr_strategy[f].str_dm.high;
			double dm_step = (double) ddtr_strategy[f].str_dm.step;
			int inbin = ddtr_strategy[f].str_dm.inBin;
			int nDMs = ddtr_strategy.m_ndms[f];
			aa_dedispersion_range newrange(dm_low, dm_high, dm_step, inbin, nTimesamples, nDMs, sampling_time);
			printf("dm_low=%f; dm_high=%f; dm_step=%f; inbin=%d; nTimesamples=%zu; nDMs=%d; sampling_time=%e;\n", dm_low, dm_high, dm_step, inbin, nTimesamples, nDMs, sampling_time); 
			ddtr_ranges.push_back(newrange);
		}
	}


	//----------------------------- Getters ---------------------------
	int nRanges() const {
		return((int) ddtr_ranges.size());
	}
	
	aa_dedispersion_range get_range(int id) const {
		if (id < (int) ddtr_ranges.size()){
			return(ddtr_ranges[id]);
		}
		else return(NULL);
	}

	float sigma_cutoff() const {
		return sigma_cutoff;
	}

	float sigma_outlier_rejection_threshold() const {
		return sigma_outlier_rejection_threshold;
	}

	/** \returns A number of harmonics used in harmonic summing algorithm. */
	int nHarmonics() const {
		return nHarmonics;
	}

	/** \returns type of the reduction algorithm for candidate selection. */
	int candidate_algorithm() const {
		return candidate_algorithm;
	}

	/** \returns A boolean indicating whether the outlier rejection will be enabled, for an instance of aa_periodicity_strategy. */
	bool outlier_rejection() const {
		return enable_outlier_rejection;
	}

	/** \returns A boolean indicating whether the interpolation will be enabled, for an instance of aa_periodicity_strategy. */
	bool interpolation() const {
		return enable_interpolation;
	}

	/** \returns A boolean indicating whether the spectrum whitening will be enabled, for an instance of aa_periodicity_strategy. */
	bool spectrum_whitening() const {
		return enable_spectrum_whitening;
	}

	/** \returns A boolean indicating whether the time-series will be padded to nearest power of two. */
	bool pad_to_nearest_higher_pow2() const {
		return(pad_to_nearest_higher_pow2);
	}
	
private:
	std::vector<aa_dedispersion_range> ddtr_ranges; /** DDTR - plan as calculated by ddtr strategy. */

};
} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_PERIODICITY_PLAN_HPP
