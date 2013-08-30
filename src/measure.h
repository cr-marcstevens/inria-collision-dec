/**
 * \file measure.h
 * \brief Measure costs of differents steps of the program. Count iterations, collisions and final tests performed. Display summary.
 *
 * time_stopwatch function measure seconds using time()
 * cycles_stopwatch measure cycles using cpu_cycles()
 *
 * \author Grégory Landais <gregory.landais@inria.fr>
 */
#ifndef MEASURE_H
#define MEASURE_H
#include "isd.h"

void total_time_stopwatch_start();
void total_time_stopwatch_stop();
void total_cycle_stopwatch_start();
void total_cycle_stopwatch_stop();
void pivot_cycle_stopwatch_start();
void pivot_cycle_stopwatch_stop();
void bday_cycle_stopwatch_start();
void bday_cycle_stopwatch_stop();
void final_test_cycle_stopwatch_start();
void final_test_cycle_stopwatch_stop();

void incr_iter_counter();
void incr_final_test_counter();
void incr_collision_counter();

void report(isd_params* params);

#endif
