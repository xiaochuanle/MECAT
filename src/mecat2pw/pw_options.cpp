#include "pw_options.h"

#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <cstdio>
#include <sstream>
#include <string>

static const int kDefaultNumThreads = 1;
static const int kDefaultNumCandidates = 100;
static const int kDefaultAlignSizePacbio = 1000;
static const int kDefaultAlignSizeNanopore = 500;
static const int kDefaultKmerMatchPacbio = 3;
static const int kDefaultKmerMatchNanopore = 2;

void
print_options(options_t* options)
{
	LOG(stderr, "task\t\t%d", options->task);
	LOG(stderr, "reads\t\t%s", options->reads);
	LOG(stderr, "output\t\t%s", options->output);
	LOG(stderr, "working folder\t%s", options->wrk_dir);
	LOG(stderr, "# of threads\t\t%d", options->num_threads);
	LOG(stderr, "# of candidates\t%d", options->num_candidates);
	LOG(stderr, "min align size\t%d", options->min_align_size);
	LOG(stderr, "min block score\t%d", options->min_kmer_match);
	LOG(stderr, "output gapped start\t%c", options->output_gapped_start_point ? 'Y' : 'N'); 
	LOG(stderr, "tech\t%d", options->tech);
}

// given options, recreate arguments from the command line
std::string
make_options(options_t* options)
{
	std::stringstream cmd;
	if (options->task == TASK_ALN || options->task == TASK_SEED) {
		cmd << " -j " << options->task;
	}
	if (options->reads != NULL) {
		cmd << " -d " << options->reads;
	}
	if (options->output != NULL) {
		cmd << " -o " << options->output;
	}
	if (options->wrk_dir != NULL) {
		cmd << " -w " << options->wrk_dir;
	}
	if (options->grid_options != NULL) {
		cmd << " -G " << options->grid_options;
	}
	if (options->grid_options_split != NULL) {
		cmd << " -S " << options->grid_options_split;
	}
	if (options->num_threads > 0) {
		cmd << " -t " << options->num_threads;
	}
	if (options->num_candidates > 0) {
		cmd << " -n " << options->num_candidates;
	}
	if (options->min_align_size > 0) {
		cmd << " -a " << options->min_align_size;
	}
	if (options->min_kmer_match > 0) {
		cmd << " -k " << options->min_kmer_match;
	}
	if (options->output_gapped_start_point > 0) {
		cmd << " -g " << options->output_gapped_start_point;
	}
	if (options->tech == TECH_PACBIO || options->tech == TECH_NANOPORE) {
		cmd << " -x " << (options->tech == TECH_PACBIO ? 0 : 1);
	}
	if (options->num_vols != -1) {
		cmd << " -N " << options->num_vols;
	}
	if (options->job_index != -1) {
		cmd << " -i " << options->job_index;
	}
	if (options->reads_to_correct) {
		cmd << " -R " << options->reads_to_correct;
	}
	if (options->binary_output) {
		cmd << " -b";
	}
	return cmd.str();
}

void init_options(options_t* const options, const int tech) {
	assert(options);
	options->task = TASK_ALN;
	options->reads = NULL;
	options->output = NULL;
	options->wrk_dir = NULL;
	options->grid_options = NULL;
	options->grid_options_split = NULL;
	options->num_threads = 1;
	options->num_candidates = 100;
	options->output_gapped_start_point = 0;
	options->tech = tech;
	options->job_index = -1;
	options->num_vols = -1;
	options->reads_to_correct = 0;
	options->binary_output = 0;
	if (tech == TECH_PACBIO) {
		options->min_align_size = kDefaultAlignSizePacbio;
		options->min_kmer_match = kDefaultKmerMatchPacbio;
	} else {
		options->min_align_size = kDefaultAlignSizeNanopore;
		options->min_kmer_match = kDefaultKmerMatchNanopore;
	}
}

void print_usage(const char* prog)
{
	fprintf(stderr, "\n\n");
	fprintf(stderr, "usage:\n");
	fprintf(stderr, "%s [-j task] [-d dataset] [-o output] [-w working dir] [-t threads] [-n candidates] [-g 0/1]", prog);
	fprintf(stderr, "\n\n");
	fprintf(stderr, "options:\n");
	fprintf(stderr, "-j <integer>\tjob: %d = seeding, %d = align\n\t\tdefault: %d\n", TASK_SEED, TASK_ALN, TASK_ALN);
	fprintf(stderr, "-d <string>\treads file name\n");
	fprintf(stderr, "-o <string>\toutput file name\n");
	fprintf(stderr, "-w <string>\tworking folder name, will be created if not exist\n");
	fprintf(stderr, "-t <integer>\tnumber of cput threads\n\t\tdefault: 1\n");
	fprintf(stderr, "-n <integer>\tnumber of candidates for gapped extension\n\t\tDefault: 100\n");
	fprintf(stderr, "-a <integer>\tminimum size of overlaps\n\t\t");
	fprintf(stderr, "Default: %d if x = %d, %d if x = %d\n", kDefaultAlignSizePacbio, TECH_PACBIO, kDefaultAlignSizeNanopore, TECH_NANOPORE);
	fprintf(stderr, "-k <integer>\tminimum number of kmer match a matched block has\n\t\t");
	fprintf(stderr, "Default: %d if x = %d, %d if x = %d\n", kDefaultKmerMatchPacbio, TECH_PACBIO, kDefaultKmerMatchNanopore, TECH_NANOPORE);
	fprintf(stderr, "-g <0/1>\twhether print gapped extension start point, 0 = no, 1 = yes\n\t\tDefault: 0\n");
	fprintf(stderr, "-x <0/x>\tsequencing technology: 0 = pacbio, 1 = nanopore\n\t\tDefault: 0\n");
	fprintf(stderr, "-R <integer>\tnumber of reads to error correct [all]\n");
	fprintf(stderr, "-G <string>\tscheduler command/options\n");
	fprintf(stderr, "-S <string>\tscheduler command/options (for split)\n");
	fprintf(stderr, "-b\t\toutput file is binary rather than text\n");
	fprintf(stderr, "(note that grid command/options have to start with \"qsub\" and end with \"-N\")\n");
	fprintf(stderr, "(note that slurm command/options have to start with \"sbatch\" and end with \"-J\")\n");
}

int parse_arguments(int argc, char* argv[], options_t* options) {
	int opt_char;
	char err_char;
	opterr = 0;
	int ret = 0;
	int task = -1;
	const char* reads = NULL;
	const char* output = NULL;
	const char* wrk_dir = NULL;
	const char* grid_options = NULL;
	const char* grid_options_split = NULL;
	int num_threads = -1;
	int num_candidates = -1;
	int min_align_size = -1;
	int min_kmer_match = -1;
	int output_gapped_start_point = -1;
	int tech = TECH_PACBIO;
	int job_index = -1;
	int num_vols = -1;
	int reads_to_correct = 0;
	int binary_output = 0;
	while((opt_char = getopt(argc, argv, "j:d:o:w:t:n:g:x:a:k:G:i:N:R:S:b")) != -1) {
		switch(opt_char) {
			case 'j':
				task = atoi(optarg);
				break;
			case 'd':
				reads = optarg;
				break;
			case 'o':
				output = optarg;
				break;
			case 'w':
				wrk_dir = optarg;
				break;
			case 'G':
				grid_options = optarg;
				break;
			case 'S':
				grid_options_split = optarg;
				break;
			case 't':
				num_threads = atoi(optarg);
				break;
			case 'i':
				job_index = atoi(optarg);
				break;
			case 'N':
				num_vols = atoi(optarg);
				break;
			case 'R':
				reads_to_correct = atoi(optarg);
				break;
			case 'n':
				num_candidates = atoi(optarg);
				break;
			case 'a':
				min_align_size = atoi(optarg);
				break;
			case 'k':
				min_kmer_match = atoi(optarg);
				break;
			case 'g':
				if (optarg[0] == '0')  {
					output_gapped_start_point = 0;
				} else if (optarg[0] == '1') {
					output_gapped_start_point = 1;
				} else {
					LOG(stderr, "argument to option \'-g\' must be either \'0\' or \'1\'");
					return 1;
				}
				break;
			case 'x':
				if (optarg[0] == '0') {
					tech = TECH_PACBIO;
				} else if (optarg[0] == '1') {
					tech = TECH_NANOPORE;
				} else {
					ERROR("invalid argument to option 'x': %s", optarg);
				}
				break;
			case '?':
				err_char = (char)optopt;
				LOG(stderr, "unrecognised option \'%c\'", err_char);
				return 1;
				break;
			case ':':
				err_char = (char)optopt;
				LOG(stderr, "argument to option \'%c\' is not provided!", err_char);
				return 1;
				break;
		}
	}
	init_options(options, tech);
	if (task != -1) {
		options->task = task;
	}
	options->reads = reads;
	options->output = output;
	options->wrk_dir = wrk_dir;
	if (grid_options != NULL) {
		options->grid_options = grid_options;
	}
	if (grid_options_split != NULL) {
		options->grid_options_split = grid_options_split;
	}
	if (num_threads != -1) {
		options->num_threads = num_threads;
	}
	if (num_candidates != -1) {
		options->num_candidates = num_candidates;
	}
	if (min_align_size != -1) {
		options->min_align_size = min_align_size;
	}
	if (min_kmer_match != -1) {
		options->min_kmer_match = min_kmer_match;
	}
	if (output_gapped_start_point != -1) {
		options->output_gapped_start_point = output_gapped_start_point;
	}
	if (job_index != -1) {
		options->job_index = job_index;
	}
	if (num_vols != -1) {
		options->num_vols = num_vols;
	}
	if (reads_to_correct) {
		options->reads_to_correct = reads_to_correct;
	}
	if (binary_output) {
		options->binary_output = 1;
	}
	if (options->task != TASK_SEED && options->task != TASK_ALN) {
		LOG(stderr, "task (-j) must be %d or %d, not %d.", TASK_SEED, TASK_ALN, options->task);
		ret = 1;
	}
	if (!options->reads) {
		LOG(stderr, "dataset must be specified.");
		ret = 1;
	} else if (!options->output) {
		LOG(stderr, "output must be specified.");
		ret = 1;
	} else if (!options->wrk_dir) {
		LOG(stderr, "working directory must be specified.");
		ret = 1;
	} else if (options->num_threads < 1) {
		LOG(stderr, "number of cpu threads must be > 0.");
		ret = 1;
	} else if (options->num_candidates < 1) {
		LOG(stderr, "number of candidates must be > 0.");
		ret = 1;
	} else if (options->reads_to_correct < 0) {
		LOG(stderr, "reads to correct must be > 0.");
		ret = 1;
	}
	if (ret) {
		return ret;
	}
	DIR* dir = opendir(options->wrk_dir);
	if (dir == NULL) {
		if (mkdir(options->wrk_dir, S_IRWXU) == -1) {
			LOG(stderr, "fail to create folder \'%s\'!", options->wrk_dir);
			exit(1);
		}
	} else {
		closedir(dir);
	}
	return ret;
}
