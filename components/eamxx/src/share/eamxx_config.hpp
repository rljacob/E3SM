#ifndef SCREAM_CONFIG_HPP
#define SCREAM_CONFIG_HPP

#include <string>

// Include this file, not any lower-level configuration file such as that
// generated by CMake from eamxx_config.h.in. The intent is to funnel all
// configuration decisions through this header.

#ifdef SCREAM_CONFIG_IS_CMAKE
# include "eamxx_config.h"
#else
// Purposely error out.
"A non-cmake build of scream is currently not supported."
#endif

namespace scream {

std::string eamxx_config_string();

// Utils to set/get whether leap year is used or not
bool use_leap_year ();
void set_use_leap_year (const bool use_leap);

// Allow downstream code to avoid macros
bool constexpr is_scream_standalone () {
#ifdef SCREAM_CIME_BUILD
  return false;
#else
  return true;
#endif
}

} // namespace scream

#endif // SCREAM_CONFIG_HPP
