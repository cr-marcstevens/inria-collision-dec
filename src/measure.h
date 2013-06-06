#ifndef MEASURE_H
#define MEASURE_H

void pivot_probe_start();
void pivot_probe_stop();
void bday_probe_start();
void bday_probe_stop();
void final_test_probe_start();
void final_test_probe_stop();
void incr_final_test_counter();
void incr_nb_collision_counter();

void display(unsigned long long nb_iter);

#endif
