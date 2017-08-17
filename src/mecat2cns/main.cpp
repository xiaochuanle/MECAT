#include "reads_correction_can.h"
#include "reads_correction_m4.h"
#include "overlaps_partition.h"
#include "options.h"

#include <list>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>

static void grid_start(char *prog, const ReadsCorrectionOptions &options, const int i) {
	// create grid command
	ReadsCorrectionOptions new_options = options;
	new_options.job_index = i;
	new_options.grid_options = NULL;
	std::ostringstream cmd;
	cmd << "qsub -N m2cns." << i << " " << options.grid_options;
	cmd << " \"" << prog << make_options(new_options) << "\"";
	assert(system(cmd.str().c_str()) == 0);
}

// pass by value so we can modify list to easily avoid checking previously
// found results files
static void wait_for_files(std::list<std::string> results) {
	const std::list<std::string>::const_iterator end_a = results.end();
	while (!results.empty()) {
		sleep(60);
		std::list<std::string>::iterator a = results.begin();
		while (a != end_a) {
			if (access(a->c_str(), F_OK) == 0) {
				a = results.erase(a);
			} else {
				++a;
			}
		}
	}
}

static void merge_results(const char *output, const std::list<std::string> &results) {
	std::list<std::string>::const_iterator a = results.begin();
	const std::list<std::string>::const_iterator end_a = results.end();
	for (; a != end_a; ++a) {
		std::ostringstream cmd;
		if (a == results.begin()) {
			cmd << "cat " << *a << " >" << output;
		} else {
			cmd << "cat " << *a << " >> " << output;
		}
		assert(system(cmd.str().c_str()) == 0);
		unlink(a->c_str());
	}
}

int main(int argc, char** argv) {
	ReadsCorrectionOptions rco;
	if (parse_arguments(argc, argv, rco)) {
		print_usage(argv[0]);
		exit(1);
	}
	if (rco.print_usage_info) {
		print_usage(argv[0]);
		exit(0);
	}
	
	// partition once up front, to allow for grid jobs
	if (rco.job_index == -1) {
		if (rco.input_type == INPUT_TYPE_CAN) {
			partition_candidates(rco.m4, rco.batch_size, rco.min_size);
		} else {
			partition_m4records(rco.m4, rco.min_mapping_ratio - 0.02, rco.batch_size, rco.min_size);
		}
	}
	if (rco.grid_options == NULL) {		// single process run
		if (rco.input_type == INPUT_TYPE_CAN) {
			return reads_correction_can(rco);
		} else {
			return reads_correction_m4(rco);
		}
	} else {				// grid run
		// get number of partitions
		std::string idx_file_name;
		generate_partition_index_file_name(rco.m4, idx_file_name);
		std::vector<PartitionFileInfo> partition_file_vec;
		load_partition_files_info(idx_file_name.c_str(), partition_file_vec);
		std::list<std::string> results;
		for (size_t i = 0; i != partition_file_vec.size(); ++i) {
			std::ostringstream os;
			os << rco.corrected_reads << "." << i;
			if (access(os.str().c_str(), F_OK) != 0) {
				grid_start(argv[0], rco, i);
				results.push_back(os.str());
			}
		}
		wait_for_files(results);
		merge_results(rco.corrected_reads, results);
		return 0;
	}
}
