#ifndef AETHER_WAVE_TRANSMISSION_H
#define AETHER_WAVE_TRANSMISSION_H

#include "symbol_export.h"

SHO_PUBLIC void evaluate_wave_transmission(int n_time, double heartrate, double a0, int no_freq, double *a, double *b, int n_bcparams, double *bc_params);


#endif /* AETHER_WAVE_TRANSMISSION_H */
