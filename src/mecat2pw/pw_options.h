#ifndef MP_OPTIONS_H
#define MP_OPTIONS_H

#include <string>

#include "../common/defs.h"

#define TASK_SEED 0
#define TASK_ALN  1

typedef struct
{
	int task;
    const char* reads;
    const char* output;
    const char* wrk_dir;
    const char* grid_options;
    int         num_threads;
    int         num_candidates;
	int			min_align_size;
	int			min_kmer_match;
    int         output_gapped_start_point;
	int 		tech;
	int 		num_vols;
	int 		job_index;
} options_t;

std::string
make_options(options_t* options);

void
print_options(options_t* options);

void
init_options(options_t* options);

int
parse_arguments(int argc, char* argv[], options_t* options);

void
print_usage(const char* prog);

#endif // MP_OPTIONS_H
