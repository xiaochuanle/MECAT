#ifndef OPTIONS_H
#define OPTIONS_H

#include <string>	// string

#define INPUT_TYPE_CAN 	0
#define INPUT_TYPE_M4	1

struct ReadsCorrectionOptions {
    int         input_type;
    const char* m4;
    const char* reads;
    const char* corrected_reads;
    const char* grid_options;
    const char* grid_options_split;
    int         num_threads;
    size_t	batch_size;
    double      min_mapping_ratio;
    int         min_align_size;
    int         min_cov;
    int		min_size;
    int		print_usage_info;
    int         tech;
    int         num_partition_files;
    int         job_index;
    int         reads_to_correct;
    int		grid_start_delay;
    int		full_reads;
    size_t	read_buffer_size;
    double	error_rate;	// .15 for pacbio, .2 for nanopore
};

void print_usage(const char* prog);

int parse_arguments(int argc, char* argv[], ReadsCorrectionOptions& t);

void make_options(const ReadsCorrectionOptions& t, std::string& new_options_str);

#endif // OPTIONS_H
