#include "pw_options.h"

#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <cstdio>

static const int kDefaultNumThreads = 1;
static const int kDefaultNumCandidates = 100;
static const int kDefaultAlignSizePacbio = 2000;
static const int kDefaultAlignSizeNanopore = 500;
static const int kDefaultKmerMatchPacbio = 4;
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

void
init_options(options_t* options, int tech)
{
    assert(options);
	options->task = TASK_ALN;
    options->reads = NULL;
    options->output = NULL;
    options->wrk_dir = NULL;
    options->num_threads = 1;
    options->num_candidates = 100;
    options->output_gapped_start_point = 0;
	options->tech = tech;
	
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
}

int
parse_arguments(int argc, char* argv[], options_t* options)
{
    int opt_char;
    char err_char;
    opterr = 0;
    int ret = 0;

	int task = -1;
	const char* reads = NULL;
	const char* output = NULL;
	const char* wrk_dir = NULL;
	int num_threads = -1;
	int num_candidates = -1;
	int min_align_size = -1;
	int min_kmer_match = -1;
	int output_gapped_start_point = -1;
	int tech = TECH_PACBIO;
    
    while((opt_char = getopt(argc, argv, "j:d:o:w:t:n:g:x:a:k:")) != -1)
    {
        switch(opt_char)
        {
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
            case 't':
                num_threads = atoi(optarg);
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
                if (optarg[0] == '0') 
                    output_gapped_start_point = 0;
                else if (optarg[0] == '1')
                    output_gapped_start_point = 1;
                else
                {
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
	if (task != -1) options->task = task;
	options->reads = reads;
	options->output = output;
	options->wrk_dir = wrk_dir;
	if (num_threads != -1) options->num_threads = num_threads;
	if (num_candidates != -1) options->num_candidates = num_candidates;
	if (min_align_size != -1) options->min_align_size = min_align_size;
	if (min_kmer_match != -1) options->min_kmer_match = min_kmer_match;
	if (output_gapped_start_point != -1) options->output_gapped_start_point = output_gapped_start_point;
	
	if (options->task != TASK_SEED && options->task != TASK_ALN)
	{
		LOG(stderr, "task (-j) must be %d or %d, not %d.", TASK_SEED, TASK_ALN, options->task);
		ret = 1;
	}

    if (!options->reads)
    {
        LOG(stderr, "dataset must be specified.");
        ret = 1;
    }
    else if (!options->output)
    {
        LOG(stderr, "output must be specified.");
        ret = 1;
    }
    else if (!options->wrk_dir)
    {
        LOG(stderr, "working directory must be specified.");
        ret = 1;
    }
    else if (options->num_threads < 1)
    {
        LOG(stderr, "number of cpu threads must be > 0.");
        ret = 1;
    }
    else if (options->num_candidates < 1)
    {
        LOG(stderr, "number of candidates must be > 0.");
        ret = 1;
    }

    if (ret) return ret;

    DIR* dir = opendir(options->wrk_dir);
    if (dir == NULL)
    {
        int t = mkdir(options->wrk_dir, S_IRWXU);
        if (t == -1)
        {
            LOG(stderr, "fail to create folder \'%s\'!", options->wrk_dir);
            exit(1);
        }
    }
    else closedir(dir);

    return ret;
}
