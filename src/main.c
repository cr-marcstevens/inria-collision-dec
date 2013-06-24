#include <stdlib.h>
#include <stdint.h>
#include <signal.h>
#include "m4ri/m4ri.h"
#include "isd.h"
#include "libisd.h"
#include "io.h"
#include "prng.h"
#include "isd_cmdline.h"

int main(int argc, char **argv)
{
	// Parse command line
	struct gengetopt_args_info args_info;
	if (cmdline_parser (argc, argv, &args_info) != 0)
		exit(EXIT_FAILURE);

	// Initialize prng
	uint64_t seed;
	if (args_info.seed_given) {
		seed = args_info.seed_arg;
	}
	else {
		seed = random_seed();
		seed &= 0x7FFFFFFFFFFFFFFFUL; // gengetopt doesn't accept unsigned longs so we choose a 63 bits seed to simplify the replay
	}

	ranctx* state = (ranctx*) malloc(sizeof(ranctx));
	
	raninit(state, seed);

  // Read CSD instance from file
	// TODO : change file format to png ?
	FILE* fd;
  if (!args_info.input_given) {
    fd=stdin;
  }
  else if((fd=fopen(args_info.input_arg, "r"))==NULL){
    fprintf(stderr, "%s : file not found\n", args_info.input_arg);
    exit(EXIT_FAILURE);
  }
	
	mzd_t* HzeroT;
	unsigned int w, N;
	word** synds;
	CSD_from_file(&HzeroT, &w, &N, &synds, fd);
	fclose(fd);

	// TODO : call workfactor if l not given
	unsigned int l = args_info.l_arg;
	if (l <= 0 || l >= (unsigned int) HzeroT->ncols){
		fprintf(stderr, "l must be >0 and <r\n");
		exit(EXIT_FAILURE);
	}

	unsigned int l2 = args_info.l2_given ? (unsigned int) args_info.l2_arg : l/2;
	/*
	if (args_info.l2_given) {
		l2 = args_info.l2_arg
	}
	else {
		l2 = l/2;
	}
*/
	if (args_info.w_given) {
		w = args_info.w_arg;
	}


	// number of bits to consider if we handle one word of data. Should be word_len except for tiny problem where r is smaller than this.
	unsigned int eff_word_len = min((unsigned int) HzeroT->ncols, word_len);

	// This threshold is an early abort criteria on the weight of the sum.
	// It introduces a probability to miss the solution 
	// p_miss = 1 - (nCr(eff_word_len - l, i)*nCr(r - eff_word_len, w-p-i) / nCr(r-l, w-p));
	// The program compute_threshold tries to compute the optimal value
	// balancing p_miss and the probability of false alarm considering the costs of the different steps.
	// The heuristic threshold = (eff_word_len - l)/4 seems good.
	unsigned int weight_threshold;
	if (args_info.threshold_given) {
		weight_threshold = args_info.threshold_arg;
	}
	else {
		weight_threshold = (eff_word_len - l)/4; // heuristic
	}

	signal(SIGUSR1, status_handler);
	signal(SIGTERM, stop_handler);
	signal(SIGINT, stop_handler);

	printf("seed : 0x%016lx\n", seed);

	isd(HzeroT, l, l2, w, N, synds, weight_threshold, args_info.max_iter_arg, args_info.max_sol_arg, args_info.max_time_arg, state, args_info.skip_arg);

	CSD_free(HzeroT, N, synds);
	free(state);

	cmdline_parser_free (&args_info);
	exit(EXIT_SUCCESS);
}
