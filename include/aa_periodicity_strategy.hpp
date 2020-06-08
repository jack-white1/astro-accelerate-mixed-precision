#ifndef ASTRO_ACCELERATE_AA_PERIODICITY_STRATEGY_HPP
#define ASTRO_ACCELERATE_AA_PERIODICITY_STRATEGY_HPP

#include <stdio.h>

#include "aa_strategy.hpp"
#include "aa_periodicity_plan.hpp"
#include "aa_dedispersion_range.hpp"
#include "aa_periodicity_processing_chunks.hpp"
#include "aa_log.hpp"
#include "aa_params.hpp"
#include "aa_device_MSD_plane_profile.hpp"

namespace astroaccelerate {

/**
 * \class aa_periodicity_strategy aa_periodicity_strategy.hpp "include/aa_periodicity_strategy.hpp"
 * \brief Class that receives an aa_periodicity_plan object, and produces an aa_periodicity_strategy object.
 * \details A periodicity strategy is required for any pipeline running the periodicity component.
 * \author AstroAccelerate team.
 * \date 23 October 2018.
 */

class aa_periodicity_strategy : public aa_strategy {
public:
	/** \struct aa_psr_chunk
	* \brief A struct for specifying chunks of data which needs to be processed.
	*/
	struct aa_psr_chunk {
		int ddtr_range_id;
		int start_DM_trial;
		int nDM_trials;
		unsigned long int MSD_workarea_size_in_bytes; // memory required for MSD plane profile including memory for decimations, blocks of partial MSD results and memory for mean and stdev of decimated planes
		unsigned long int MSD_profile_size_in_bytes; // memory required for final result of the MSD plane profile with interpolated MSD values for all widths
		unsigned long int MSD_DIT_profile_size_in_bytes; = nDecimations*MSD_PARTIAL_SIZE*sizeof(float);
	};

	// TODO: process ranges individually
	// TODO: separate each range into processable chunks
	// TODO: pre-calculate everything beforehand

	// Steps in periodicity:
	// get chunk from the host
	// fourier transform the chunks
	// calculate power spectrum -- with callbacks
	// apply de-rednning
	// MSD plane profile
	// harmonic summing
	// candidate search

	// Correct this!!!
	/** \brief Constructor for aa_periodicity_strategy that sets all member variables upon construction. */
	aa_periodicity_strategy(const aa_periodicity_plan &plan, size_t available_memory) {
		sigma_cutoff = plan.sigma_cutoff();
		sigma_outlier_rejection_threshold = plan.sigma_outlier_rejection_threshold;
		nHarmonics = plan.nHarmonics;
		candidate_algorithm = plan.candidate_algorithm;
		enable_outlier_rejection = plan.enable_outlier_rejection;
		enable_interpolation = plan.enable_interpolation;
		enable_spectrum_whitening = plan.enable_spectrum_whitening;
		pad_to_nearest_higher_pow2 = plan.pad_to_nearest_higher_pow2;
		
		ready = false;
		
		unsigned long int input_memory_required;
		unsigned long int cuFFT_memory_require;
		unsigned long int spectrum_whitening_memory_required = 0;
		unsigned long int MSD_profile_size_in_bytes;
		unsigned long int MSD_DIT_profile_size_in_bytes;
		unsigned long int MSD_workarea_size_in_bytes;
		unsigned long int HRMS_memory_required = 0;
		unsigned long int candidate_search_memory_required;
		
		unsigned long int max_memory_required = 
		
		if ((nHarmonics > 0) && (sigma_constant < 0)) {
			m_ready = true;
		} 
		else {
			LOG(log_level::warning, "Invalid periodicity strategy parameters. Check the aa_periodicity_plan input parameters.");
			print_info(*this);
		}
		
		
		
		


	bool m_ready; /** Ready state of the instance. */
	}
	

	/** Static member function that prints member variables for a provided aa_periodicity_strategy. */
	static bool print_info(const aa_periodicity_strategy &strategy) {
		LOG(log_level::dev_debug, "--------------> PERIODICITY STRATEGY INFORMATION <--------------");
		LOG(log_level::dev_debug, "PSR - Sigma cutoff for candidate selection:\t\t" + std::to_string(strategy.sigma_cutoff()));
		LOG(log_level::dev_debug, "PSR - Enable outlier rejection:\t\t" + (strategy.candidate_algorithm() ? std::string("true") : std::string("false")));
		LOG(log_level::dev_debug, "PSR - Sigma cutoff for outlier rejection:\t\t" + std::to_string(strategy.sigma_constant()));
		LOG(log_level::dev_debug, "PSR - Number of harmonics:\t\t\t" + std::to_string(strategy.nHarmonics()));
		LOG(log_level::dev_debug, "PSR - Candidate selection algorithm:\t\t" + std::to_string(strategy.export_powers()));
		LOG(log_level::dev_debug, "PSR - Enable interpolation:\t\t" + (strategy.candidate_algorithm() ? std::string("true") : std::string("false")));
		LOG(log_level::dev_debug, "PSR - Enable spectrum whitening:\t" + (strategy.enable_msd_baseline_noise() ? std::string("true") : std::string("false")));
		LOG(log_level::dev_debug, "PSR - Pad to the nearest higher power of 2:\t" + (strategy.enable_msd_baseline_noise() ? std::string("true") : std::string("false")));
		LOG(log_level::dev_debug, "\n");
		LOG(log_level::dev_debug, "PSR - DDTR ranges:\n");
		for (int f = 0; f < ddtr_ranges.size(); f++) {
			LOG(log_level::dev_debug, "DM range: " + std::to_string(ddtr_ranges[f].dm_low()) + " -- " + std::to_string(ddtr_ranges[f].dm_high()) + " step: " + std::to_string(ddtr_ranges[f].dm_step()) + "; binning: " + std::to_string(ddtr_ranges[f].inBin()) + "; time-samples: " + std::to_string(ddtr_ranges[f].nTimesamples()) + "; DM-trials: " + std::to_string(ddtr_ranges[f].nDMs()) + ";\n");
		}
		LOG(log_level::dev_debug, "\n");
		LOG(log_level::dev_debug, "PSR - Processing chunks:\n");
		for (int f = 0; f < psr_chunks.size(); f++) {
			LOG(log_level::dev_debug, "DM range: " + std::to_string(ddtr_ranges[f].dm_low()) + " -- " + std::to_string(ddtr_ranges[f].dm_high()) + " step: " + std::to_string(ddtr_ranges[f].dm_step()) + "; binning: " + std::to_string(ddtr_ranges[f].inBin()) + "; time-samples: " + std::to_string(ddtr_ranges[f].nTimesamples()) + "; DM-trials: " + std::to_string(ddtr_ranges[f].nDMs()) + ";\n");
		}
		return true;
	}
	
	/** \returns The ready state of the instance of the class. */
	bool ready() const {
		return m_ready;
	}

	/** \brief Performs any remaining setup needed.
	 * \returns Whether the operation was successful or not, at the moment this is the ready state of the instance. */
	bool setup() {
		return ready();
	}

	/** \returns The name of the module. */
	std::string name() const {
		return "periodicity_strategy";
	}



	//----------------------------- Getters ---------------------------
	int nRanges() const {
		return((int)ddtr_ranges.size());
	}

	aa_dedispersion_range get_range(int id) const {
		if (id < (int)ddtr_ranges.size()) {
			return(ddtr_ranges[id]);
		} else return(NULL);
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
	std::vector<aa_psr_chunk> psr_chunks; /** Chunks of data to be processed by periodicity search as determined by psr strategy. */
	std::vector<aa_dedispersion_range> ddtr_ranges; /** DDTR - plan as calculated by ddtr strategy. */
	float sigma_cutoff; /** User selected value of sigma cutoff. Any event with SNR (sigma) below this value will not be selected as candidate */
	float sigma_outlier_rejection_threshold; /** User selected sigma for outlier rejection. Any value with sigma greater than this will be rejected from calculation of mean and standard deviation. */
	int   nHarmonics; /** User selected number of harmonics performed by the periodicity search. */
	int   candidate_algorithm; /** User selected flag to select reduction algorithm to select candidates. */
	bool  enable_outlier_rejection; /** User selected flag to enable or disable outlier rejection when calculating mean and standard deviation. */
	bool  enable_interpolation; /** Enable or disable interpolation for periodicity search */
	bool  enable_spectrum_whitening; /** Enable or disable spectrum whitening (removal of the red noise) for periodicity search */
	bool  pad_to_nearest_higher_pow2; /** Whether the periodicity will pad data to the nearest power of two for the Fourier transform. True by default */

	bool ready; /** Ready state of the instance. */
	
	
	void prepare_strategy(){
		// Memory requirements
		int MSD_memory_reguired_per_sample;
		int cuFFT_memory_required_per_sample;
		int HRMS_memory_required_per_sample;
	}
	
	
	
};

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_PERIODICITY_STRATEGY_HPP
