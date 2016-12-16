#ifndef READS_CORRECTION_OPTIONS_H
#define READS_CORRECTION_OPTIONS_H

#include "../common/defs.h"


#define INPUT_TYPE_CAN 	0
#define INPUT_TYPE_M4	1

struct ReadsCorrectionOptions
{
	int 		input_type;
    const char* m4;
    const char* reads;
	const char* corrected_reads;
    int         num_threads;
    index_t     batch_size;
    double      min_mapping_ratio;
    int         min_cov;
    index_t     min_size;
    bool        print_usage_info;
};

void init_reads_correction_options(ReadsCorrectionOptions& rco);

void parse_reads_correction_arguments(int argc, char** argv, ReadsCorrectionOptions& rco);

void print_reads_correction_usage(const char* prog_name);

std::ostream& operator<<(std::ostream& out, const ReadsCorrectionOptions& rco);

#endif // READS_CORRECTION_OPTIONS_H
