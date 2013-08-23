#ifndef MEASURE_H
#define MEASURE_H
#include "isd.h"

void total_stopwatch_start();
void total_stopwatch_stop();
void total_probe_start();
void total_probe_stop();
void pivot_probe_start();
void pivot_probe_stop();
void bday_probe_start();
void bday_probe_stop();
void final_test_probe_start();
void final_test_probe_stop();

void incr_iter_counter();
void incr_final_test_counter();
void incr_collision_counter();

void report(isd_params* params);

#endif
