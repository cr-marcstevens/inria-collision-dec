#ifndef CUSTOM_BRILLIANTRUSSIAN_H
#define CUSTOM_BRILLIANTRUSSIAN_H

#include "m4ri/brilliantrussian.h"
#include "m4ri/misc.h"

rci_t _mzd_partial_echelonize_m4ri(mzd_t *A, int const full, int k, int heuristic, double const threshold, rci_t nb_col);

#endif
