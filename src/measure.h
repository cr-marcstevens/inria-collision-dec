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

void get_costs(unsigned long long nb_iter, long long* pivot_cost, long long* bday_cost, long long* final_test_cost);
void display(unsigned long long nb_iter);

#endif
