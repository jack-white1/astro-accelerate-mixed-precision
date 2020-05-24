#ifndef ASTRO_ACCELERATE_AA_PERIODICITY_STRATEGY_HPP
#define ASTRO_ACCELERATE_AA_PERIODICITY_STRATEGY_HPP

#include <stdio.h>

#include "aa_strategy.hpp"
#include "aa_periodicity_plan.hpp"
#include "aa_log.hpp"

namespace astroaccelerate {

/**
 * \class aa_periodicity_strategy aa_periodicity_strategy.hpp "include/aa_periodicity_strategy.hpp"
 * \brief Class that receives an aa_periodicity_plan object, and produces an aa_periodicity_strategy object.
 * \details A periodicity strategy is required for any pipeline running the periodicity component.
 * \author Cees Carels.
 * \date 23 October 2018.
 */

class aa_periodicity_strategy : public aa_strategy {
public:

	/** \brief Trivial constructor for aa_periodicity_strategy, which can never have a ready state equal to true. */
	aa_periodicity_strategy() :
		m_sigma_cutoff(0),
		m_sigma_constant(0),
		m_nHarmonics(0),
		m_export_powers(0),
		m_candidate_algorithm(false),
		m_enable_msd_baseline_noise(false),
		m_ready(false) {

	}

	/** \brief Constructor for aa_periodicity_strategy that sets all member variables upon construction. */
	aa_periodicity_strategy(const aa_periodicity_plan &plan, size_t available_memory) :
		m_sigma_cutoff(plan.sigma_cutoff()),
		m_sigma_constant(plan.sigma_constant()),
		m_nHarmonics(plan.nHarmonics()),
		m_export_powers(plan.export_powers()),
		m_candidate_algorithm(plan.candidate_algorithm()),
		m_enable_msd_baseline_noise(plan.enable_msd_baseline_noise()),
		m_ready(false) {
		/** Parse user input, if the user input is not valid, then the ready state will not become true. */
		if ((m_nHarmonics > 0) && (m_sigma_constant > 0) && (m_export_powers >= 0)) {
			m_ready = true;
		} else {
			LOG(log_level::warning, "Invalid periodicity strategy parameters. Check the aa_periodicity_plan input parameters.");
			print_info(*this);
		}
	}
	
		/** \brief Print the member data of the instnace.
	 * \returns A boolean flag to indicate whether the operation completed successfully (true) or unsuccessfully (false). */
	bool print_parameters() const {
		printf("Periodicity - sigma_cutoff %f\n", m_sigma_cutoff);
		printf("Periodicity - sigma_constant %f\n", m_sigma_constant);
		printf("Periodicity - nHarmonics %d\n", m_nHarmonics);
		printf("Periodicity - export_powers %d\n", m_export_powers);
		printf("Periodicity - candidate_algorithm %d\n", m_candidate_algorithm);
		printf("Periodicity - enable_msd_baseline_noise %d\n", m_enable_msd_baseline_noise);
		return true;
	}

	/** Static member function that prints member variables for a provided aa_periodicity_strategy. */
	static bool print_info(const aa_periodicity_strategy &strategy) {
		LOG(log_level::dev_debug, "PERIODICITY STRATEGY INFORMATION:");
		LOG(log_level::dev_debug, "periodicity sigma_cutoff:\t\t" + std::to_string(strategy.sigma_cutoff()));
		LOG(log_level::dev_debug, "periodicity sigma_constant:\t\t" + std::to_string(strategy.sigma_constant()));
		LOG(log_level::dev_debug, "periodicity nHarmonics:\t\t\t" + std::to_string(strategy.nHarmonics()));
		LOG(log_level::dev_debug, "periodicity export_powers:\t\t" + std::to_string(strategy.export_powers()));
		LOG(log_level::dev_debug, "periodicity candidate_algorithm:\t\t" + (strategy.candidate_algorithm() ? std::string("true") : std::string("false")));
		LOG(log_level::dev_debug, "periodicity enable_msd_baseline_noise:\t" + (strategy.enable_msd_baseline_noise() ? std::string("true") : std::string("false")));
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






	//------------------ Getters 
	
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
	float sigma_cutoff; /** User selected value of sigma cutoff. Any event with SNR (sigma) below this value will not be selected as candidate */
	float sigma_outlier_rejection_threshold; /** User selected sigma for outlier rejection. Any value with sigma greater than this will be rejected from calculation of mean and standard deviation. */
	int   nHarmonics; /** User selected number of harmonics performed by the periodicity search. */
	int   candidate_algorithm; /** User selected flag to select reduction algorithm to select candidates. */
	bool  enable_outlier_rejection; /** User selected flag to enable or disable outlier rejection when calculating mean and standard deviation. */
	bool  enable_interpolation; /** Enable or disable interpolation for periodicity search */
	bool  enable_spectrum_whitening; /** Enable or disable spectrum whitening (removal of the red noise) for periodicity search */
	bool  pad_to_nearest_higher_pow2; /** Whether the periodicity will pad data to the nearest power of two for the Fourier transform. True by default */

	bool m_ready; /** Ready state of the instance. */
	
	
	void prepare_strategy(){
		// Memory requirements
		int MSD_memory_reguired_per_sample;
		int cuFFT_memory_required_per_sample;
		int HRMS_memory_required_per_sample;
	}
	
	
	
};

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_PERIODICITY_STRATEGY_HPP
