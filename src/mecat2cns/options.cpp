#include "options.h"

#include <cstring>
#include <unistd.h>

#include <iostream>

using namespace std;

static int input_type_pacbio 		= 1;
static int num_threads_pacbio		= 1;
static index_t batch_size_pacbio	= 100000;
static double mapping_ratio_pacbio  = 0.9;
static int align_size_pacbio		= 2000;
static int cov_pacbio				= 6;
static int min_size_pacbio			= 5000;
static bool print_usage_pacbio		= false;
static int tech_pacbio				= TECH_PACBIO;
static int num_partition_files   	= 10;

static int input_type_nanopore 		    = 1;
static int num_threads_nanopore		    = 1;
static index_t batch_size_nanopore	    = 100000;
static double mapping_ratio_nanopore    = 0.4;
static int align_size_nanopore		    = 400;
static int cov_nanopore				    = 6;
static int min_size_nanopore			= 2000;
static bool print_usage_nanopore		= false;
static int tech_nanopore				= TECH_NANOPORE;

static int default_tech = TECH_PACBIO;

static const char input_type_n    = 'i';
static const char num_threads_n   = 't';
static const char batch_size_n    = 'p';
static const char mapping_ratio_n = 'r';
static const char align_size_n    = 'a';
static const char cov_n           = 'c';
static const char min_size_n      = 'l';
static const char usage_n         = 'h';
static const char tech_n          = 'x';
static const char num_partition_files_n = 'k';

void
print_pacbio_default_options()
{
	cerr << '-' << input_type_n << ' ' << input_type_pacbio 
		 << ' '
		 << '-' << num_threads_n << ' ' << num_threads_pacbio
		 << ' '
		 << '-' << batch_size_n << ' ' << batch_size_pacbio
		 << ' '
		 << '-' << mapping_ratio_n << ' ' << mapping_ratio_pacbio
		 << ' '
		 << '-' << align_size_n << ' ' << align_size_pacbio
		 << ' '
		 << '-' << cov_n << ' ' << cov_pacbio
		 << ' '
		 << '-' << min_size_n << ' ' << min_size_pacbio
		 << ' '
		 << '-' << num_partition_files_n << ' ' << num_partition_files
		 << "\n";
}

void print_nanopore_default_options()
{
	cerr << '-' << input_type_n << ' ' << input_type_nanopore
		 << ' '
		 << '-' << num_threads_n << ' ' << num_threads_nanopore
		 << ' '
		 << '-' << batch_size_n << ' ' << batch_size_nanopore
		 << ' '
		 << '-' << mapping_ratio_n << ' ' << mapping_ratio_nanopore
		 << ' '
		 << '-' << align_size_n << ' ' << align_size_nanopore
		 << ' '
		 << '-' << cov_n << ' ' << cov_nanopore
		 << ' '
		 << '-' << min_size_n << ' ' << min_size_nanopore
		 << ' '
		 << '-' << num_partition_files_n << ' ' << num_partition_files
		 << "\n";
}

void
print_usage(const char* prog)
{
	cerr << "USAGE:\n"
		 << prog << ' ' << "[options]" << ' ' << "input" << ' ' << "reads" << ' ' << "output" << "\n";
	cerr << "\n" << "OPTIONS:" << "\n";
	
	cerr << "-" << tech_n << " <0/1>\t" << "sequencing platform: 0 = PACBIO, 1 = NANOPORE" << "\n"
		 << "\t\t" << "default: 0" << "\n";
	
	cerr << "-" << input_type_n << " <0/1>\t" << "input type: 0 = candidte, 1 = m4" << "\n";
	
	cerr << "-" << num_threads_n << " <Integer>\t" << "number of threads (CPUs)" << "\n";
	
	cerr << "-" << batch_size_n << " <Integer>\t" << "batch size that the reads will be partitioned" << "\n";
	
	cerr << "-" << mapping_ratio_n << " <Real>\t" << "minimum mapping ratio" << "\n";
	
	cerr << "-" << align_size_n << " <Integer>\t" << "minimum overlap size" << "\n";
	
	cerr << "-" << cov_n << " <Integer>\t" << "minimum coverage under consideration" << "\n";
	
	cerr << "-" << min_size_n << " <Integer>\t" << "minimum length of corrected sequence" << "\n";
	
	cerr << "-" << num_partition_files_n << " <Integer>\t" 
		 << "number of partition files when partitioning overlap results" 
		 << " (if < 0, then it will be set to system limit value)"
		 << "\n";
	
	cerr << "-" << usage_n << "\t\t" << "print usage info." << "\n";
	
	cerr << "\n"
		 << "If 'x' is set to be '0' (pacbio), then the other options have the following default values: \n";
	print_pacbio_default_options();
	
	cerr << "\n"
		 << "If 'x' is set to be '1' (nanopore), then the other options have the following default values: \n";
	print_nanopore_default_options();
}

ConsensusOptions
init_consensus_options(int tech)
{
	ConsensusOptions t;
	if (tech == TECH_PACBIO) {
		t.input_type            = input_type_pacbio;
		t.m4                    = NULL;
		t.reads                 = NULL;
		t.corrected_reads       = NULL;
		t.num_threads           = num_threads_pacbio;
		t.batch_size            = batch_size_pacbio;
		t.min_mapping_ratio     = mapping_ratio_pacbio;
		t.min_align_size        = align_size_pacbio;
		t.min_cov               = cov_pacbio;
		t.min_size              = min_size_pacbio;
		t.print_usage_info      = print_usage_pacbio;
		t.num_partition_files 	= num_partition_files;
		t.tech                  = tech_pacbio;
	} else {
		t.input_type            = input_type_nanopore;
		t.m4                    = NULL;
		t.reads                 = NULL;
		t.corrected_reads       = NULL;
		t.num_threads           = num_threads_nanopore;
		t.batch_size            = batch_size_nanopore;
		t.min_mapping_ratio     = mapping_ratio_nanopore;
		t.min_align_size        = align_size_nanopore;
		t.min_cov               = cov_nanopore;
		t.min_size              = min_size_nanopore;
		t.print_usage_info      = print_usage_nanopore;
		t.num_partition_files	= num_partition_files;
		t.tech                  = tech_nanopore;
	}
    return t;
}

int detect_tech(int argc, char* argv[])
{
    int t = default_tech;
    char tech_nstr[64]; tech_nstr[0] = '-'; tech_nstr[1] = tech_n; tech_nstr[2] = '\0';

    for (int i = 0; i < argc; ++i) {
        if (strcmp(tech_nstr, argv[i]) == 0) {
            if (i + 1 == argc) {
                fprintf(stderr, "argument to option '%c' is missing.\n", tech_n);
                t = -1;
            }
			cout << tech_nstr << "\n";
            if (argv[i + 1][0] == '0') {
                t = TECH_PACBIO;
            } else if (argv[i + 1][0] == '1') {
                t = TECH_NANOPORE;
            } else {
                fprintf(stderr, "invalid argument to option '%c': %s\n", tech_n, argv[i + 1]);
                t = -1;
            }
            break;
        }
    }

    return t;
}

int
parse_arguments(int argc, char* argv[], ConsensusOptions& t)
{
	bool parse_success = true;
	int tech = detect_tech(argc, argv);
	if (tech == -1) {
		return 1;
	} 
	t = init_consensus_options(tech);
	
	int opt_char;
    char err_char;
    opterr = 0;
	while((opt_char = getopt(argc, argv, "i:t:p:r:a:c:l:x:k:h")) != -1) {
		switch (opt_char) {
			case input_type_n:
				if (optarg[0] == '0')
					t.input_type = INPUT_TYPE_CAN;
				else if (optarg[0] == '1')
					t.input_type = INPUT_TYPE_M4;
				else {
					fprintf(stderr, "invalid argument to option '%c': %s\n", input_type_n, optarg);
					return 1;
				}
				break;
			case num_threads_n:
				t.num_threads = atoi(optarg);
				break;
			case batch_size_n:
				t.batch_size = atoll(optarg);
				break;
			case mapping_ratio_n:
				t.min_mapping_ratio = atof(optarg);
				break;
			case align_size_n:
				t.min_align_size = atoi(optarg);
				break;
			case cov_n:
				t.min_cov = atoi(optarg);
				break;
			case min_size_n:
				t.min_size = atoll(optarg);
				break;
			case usage_n:
				t.print_usage_info = true;
				break;
			case tech_n:
				break;
			case num_partition_files_n:
				t.num_partition_files = atoi(optarg);
				break;
			case '?':
                err_char = (char)optopt;
				fprintf(stderr, "unrecognised option '%c'\n", err_char);
                return 1;
                break;
            case ':':
                err_char = (char)optopt;
				fprintf(stderr, "argument to option '%c' is missing.\n", err_char);
                return 1;
                break;
		}
	}
	
	if (t.num_threads <= 0)
	{
		std::cerr << "cpu threads must be greater than 0\n";
		parse_success = false;
	}
	if (t.batch_size <= 0)
	{
		std::cerr << "batch size must be greater than 0\n";
		parse_success = false;
	}
	if (t.min_mapping_ratio < 0.0)
	{
		std::cerr << "mapping ratio must be >= 0.0\n";
		parse_success = false;
	}
	if (t.min_cov < 0)
	{
		std::cerr << "coverage must be >= 0\n";
		parse_success = false;
	}
	if (t.min_cov < 0)
	{
		std::cerr << "sequence size must be >= 0\n";
		parse_success = false;
	}
	
	if (argc < 3) return 1;
	
	t.m4 = argv[argc - 3];
	t.reads = argv[argc - 2];
	t.corrected_reads = argv[argc - 1];
	
	if (parse_success) return 0;
	return 1;
}

void
print_options(ConsensusOptions& t)
{
	cout << "input_type:\t"	<< t.input_type << "\n";
	if (t.m4) cout << "reads\t" << t.m4 << "\n";
	if (t.reads) cout << "output\t" << t.reads << "\n";
	if (t.corrected_reads) cout << "m4\t" << t.corrected_reads << "\n";
	cout << "number of threads:\t" << t.num_threads << "\n";
	cout << "batch size:\t" << t.batch_size << "\n";
	cout << "mapping ratio:\t" << t.min_mapping_ratio << "\n";
	cout << "align size:\t" << t.min_align_size << "\n";
	cout << "cov:\t" << t.min_cov << "\n";
	cout << "min size:\t" << t.min_size << "\n";
	cout << "partition files:\t" << t.num_partition_files << "\n";
	cout << "tech:\t" << t.tech << "\n";
}
