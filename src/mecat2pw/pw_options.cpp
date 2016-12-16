#include "pw_options.h"

#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <cstdio>

static const int kDefaultNumThreads = 1;
static const int kDefaultNumCandidates = 100;

void
print_options(options_t* options)
{
	LOG(stderr, "task\t\t%d", options->task);
	LOG(stderr, "reads\t\t%s", options->reads);
	LOG(stderr, "output\t\t%s", options->output);
	LOG(stderr, "working folder\t%s", options->wrk_dir);
	LOG(stderr, "# of threads\t\t%d", options->num_threads);
	LOG(stderr, "# of candidates\t%d", options->num_candidates);
	LOG(stderr, "output gapped start\t%c", options->output_gapped_start_point ? 'Y' : 'N'); 
}

void
init_options(options_t* options)
{
    assert(options);
	options->task = TASK_ALN;
    options->reads = NULL;
    options->output = NULL;
    options->wrk_dir = NULL;
    options->num_threads = 1;
    options->num_candidates = 100;
    options->output_gapped_start_point = 0;
}

void print_usage(const char* prog)
{
	fprintf(stderr, "\n\n");
	fprintf(stderr, "usage:\n");
	fprintf(stderr, "%s [-j task] [-d dataset] [-o output] [-w working dir] [-t threads] [-n candidates] [-p 0/1]", prog);
	fprintf(stderr, "\n\n");
	fprintf(stderr, "options:\n");
	fprintf(stderr, "-j <integer>\tjob: %d = seeding, %d = align\n\t\tdefault: %d\n", TASK_SEED, TASK_ALN, TASK_ALN);
	fprintf(stderr, "-d <string>\treads file name\n");
	fprintf(stderr, "-o <string>\toutput file name\n");
	fprintf(stderr, "-w <string>\tworking folder name, will be created if not exist\n");
	fprintf(stderr, "-t <integer>\tnumber of cput threads\n\t\tdefault: 1\n");
	fprintf(stderr, "-n <integer>\tnumber of candidates for gapped extension\n\t\tDefault: 100\n");
	fprintf(stderr, "-g <0/1>\twhether print gapped extension start point, 0 = no, 1 = yes\n\t\tDefault: 0\n");
}

int
parse_arguments(int argc, char* argv[], options_t* options)
{
    init_options(options);
    int opt_char;
    char err_char;
    opterr = 0;
    int ret = 0;
    
    while((opt_char = getopt(argc, argv, "j:d:o:w:t:n:g:")) != -1)
    {
        switch(opt_char)
        {
			case 'j':
				options->task = atoi(optarg);
				break;
            case 'd':
                options->reads = optarg;
                break;
            case 'o':
                options->output = optarg;
                break;
            case 'w':
                options->wrk_dir = optarg;
                break;
            case 't':
                options->num_threads = atoi(optarg);
                break;
            case 'n':
                options->num_candidates = atoi(optarg);
                break;
            case 'g':
                if (optarg[0] == '0') 
                    options->output_gapped_start_point = 0;
                else if (optarg[0] == '1')
                    options->output_gapped_start_point = 1;
                else
                {
                    LOG(stderr, "argument to option \'-g\' must be either \'0\' or \'1\'");
                    return 1;
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
