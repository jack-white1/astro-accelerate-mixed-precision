#ifndef ASTRO_ACCELERATE_AA_HOST_RFI_HPP
#define ASTRO_ACCELERATE_AA_HOST_RFI_HPP

#include <vector>

namespace astroaccelerate {

  /** \brief Function that performs RFI mitigation. */
  void rfi(int nsamp, int nchans, std::vector<unsigned short> &input_buffer);
  void input_data_renormalization(size_t nsamp, size_t nchans, std::vector<unsigned short> &input_buffer, bool enable_outlier_rejection, float sigma);
  void dedisperse_time_renormalization(size_t nDMs, size_t nTimesamples, float *input_buffer, bool enable_outlier_rejection, float sigma);
  void dedisperse_DM_renormalization(size_t nDMs, size_t nTimesamples, float *input_buffer, bool enable_outlier_rejection, float sigma);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_HOST_RFI_HPP

