#include "reads_correction_option.h"

#include <map>
#include <string>

#include "argument.h"

/// input type
static const std::string input_type_name = "i";
static const std::string input_type_desc = "\tinput type, 0 = can, 1 = m4\n"
										   "\t\tdefault: 1";
static const int 		input_type_default = 1;
static IntegerArgument 	input_type_arg(&input_type_name, &input_type_desc, input_type_default);

/// number of threads
static const std::string num_threads_name = "t";
static const std::string num_threads_desc = "\tnumber of cpu threads used for consensus\n"
                                            "\t\tdefault: 1";
static const int         num_threads_default = 1;
static IntegerArgument num_threads_arg(&num_threads_name, &num_threads_desc, num_threads_default);

/// batch size
static const std::string batch_size_name = "p";
static const std::string batch_size_desc = "\tbatch size that the reads will be partitioned\n"
                                           "\t\tdefault: 100000";
static const index_t     batch_size_default = 100000;
static IntegerArgument batch_size_arg(&batch_size_name, &batch_size_desc, batch_size_default);

/// mapping range ratio
static const std::string min_mapping_ratio_name = "r";
static const std::string min_mapping_ratio_desc = "\tminimum mapping ratio\n"
                                                  "\t\tdefault: 0.9";
static const double      min_mapping_ratio_default = 0.9;
static DoubleArgument min_mapping_ratio_arg(&min_mapping_ratio_name, &min_mapping_ratio_desc, min_mapping_ratio_default);

/// min cov
static const std::string min_cov_name = "c";
static const std::string min_cov_desc = "\tminimum coverage under consideration\n"
                                        "\t\tdefault: 6";
static const int         min_cov_default = 6;
static IntegerArgument min_cov_arg(&min_cov_name, &min_cov_desc, min_cov_default);

/// min size
static const std::string min_size_name = "l";
static const std::string min_size_desc = "\tminimum length of corrected sequence\n"
                                         "\t\tdefault: 5000";
static const index_t     min_size_default = 5000;
static IntegerArgument min_size_arg(&min_size_name, &min_size_desc, min_size_default);

/// help
static const std::string help_name = "h";
static const std::string help_desc = "\t\tprint usage info.";
static const bool help_default = false;
static BooleanArgument help_arg(&help_name, &help_desc, help_default, false);    

void init_reads_correction_options(ReadsCorrectionOptions& rco)
{
	rco.input_type = input_type_default;
    rco.m4 = NULL;
    rco.reads = NULL;
	rco.corrected_reads = NULL;
    rco.num_threads = num_threads_default;
    rco.batch_size = batch_size_default;
    rco.min_mapping_ratio = min_mapping_ratio_default;
    rco.min_cov = min_cov_default;
    rco.min_size = min_size_default;
    rco.print_usage_info = help_default;
}

void parse_reads_correction_arguments(int argc, char** argv, ReadsCorrectionOptions& rco)
{
    std::map<std::string, Argument*> args;
	args.insert(std::make_pair(input_type_name, &input_type_arg));
    args.insert(std::make_pair(num_threads_name, &num_threads_arg));
    args.insert(std::make_pair(batch_size_name, &batch_size_arg));
    args.insert(std::make_pair(min_mapping_ratio_name, &min_mapping_ratio_arg));
    args.insert(std::make_pair(min_cov_name, &min_cov_arg));
    args.insert(std::make_pair(min_size_name, &min_size_arg));
    args.insert(std::make_pair(help_name, &help_arg));

	bool parse_success = true;
	int ind = 1;
	while (ind < argc)
	{
		if (argv[ind][0] != '-') break;
		std::string arg(argv[ind] + 1);
		if (args.find(arg) == args.end())
		{
			std::cerr << "Unrecognised option \'" << argv[ind] << "\'\n";
			parse_success = false;
			break;
		}
		++ind;
		int t = args[arg]->ProcessArgument(argc - ind, argv + ind);
		ind += t;
	}
	
	if (parse_success && help_arg.value())
	{
		print_reads_correction_usage(argv[0]);
		exit(0);
	}
	
	if (parse_success && (argc - ind != 3))
	{
		std::cerr << "Error: You must provide the 'm4-file', the 'reads-file', and the 'corrected-reads-file'\n";
		parse_success = false;
	}

	rco.input_type = input_type_arg.value();
    rco.num_threads = num_threads_arg.value();
    rco.batch_size = batch_size_arg.value();
    rco.min_mapping_ratio = min_mapping_ratio_arg.value();
    rco.min_cov = min_cov_arg.value();
    rco.min_size = min_size_arg.value();
    rco.print_usage_info = help_arg.value();
	
	if (rco.input_type != INPUT_TYPE_CAN && rco.input_type != INPUT_TYPE_M4)
	{
		std::cerr << "input type must be either " << INPUT_TYPE_CAN << " or " << INPUT_TYPE_M4 << "\n";
		parse_success = false;
	}
	if (rco.num_threads <= 0)
	{
		std::cerr << "cpu threads must be greater than 0\n";
		parse_success = false;
	}
	if (rco.batch_size <= 0)
	{
		std::cerr << "batch size must be greater than 0\n";
		parse_success = false;
	}
	if (rco.min_mapping_ratio < 0.0)
	{
		std::cerr << "mapping ratio must be >= 0.0\n";
		parse_success = false;
	}
	if (rco.min_cov < 0)
	{
		std::cerr << "coverage must be >= 0\n";
		parse_success = false;
	}
	if (rco.min_cov < 0)
	{
		std::cerr << "sequence size must be >= 0\n";
		parse_success = false;
	}
	
	if (!parse_success)
	{
		print_reads_correction_usage(argv[0]);
		exit(1);
	}
	
	rco.m4 = argv[ind++];
	rco.reads = argv[ind++];
	rco.corrected_reads = argv[ind++];
}

std::ostream& operator<<(std::ostream& out, const ReadsCorrectionOptions& rco)
{
	out << "m4\t" << rco.m4 << "\n"
		<< "reads\t" << rco.reads << "\n"
		<< "num_threads\t" << rco.num_threads << "\n"
		<< "batch_size\t" << rco.batch_size << "\n"
		<< "min_mapping_ratio\t" << rco.min_mapping_ratio << "\n"
		<< "min_cov\t" << rco.min_cov << "\n"
		<< "min_size\t" << rco.min_size << "\n"
		<< "print_usage_info\t" << rco.print_usage_info << "\n";
	return out;
}

void print_reads_correction_usage(const char* prog_name)
{
	std::cerr << "usage:\n"
			  << prog_name << " [options] m4-file reads-file corrected-reads-file\n"
			  << "\noptions:\n\n"
			  << "-" << input_type_name << " <Integer>" << input_type_desc << "\n"
			  << "-" << num_threads_name << " <Integer>" << num_threads_desc << "\n"
			  << "-" << batch_size_name << " <Integer>" << batch_size_desc << "\n"
			  << "-" << min_mapping_ratio_name << " <Real>" << min_mapping_ratio_desc << "\n"
			  << "-" << min_cov_name << " <Integer>" << min_cov_desc << "\n"
			  << "-" << min_size_name << " <Integer>" << min_size_desc << "\n"
			  << "-" << help_name << help_desc << "\n";
}
