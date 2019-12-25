#include "pw_options.h"
#include "pw_impl.h"
#include "../common/split_database.h"

#include <cstdio>
#include <fstream>
#include <unistd.h>	// ... unlink()
#include <limits.h>	// PATH_MAX
#include <fcntl.h>	// S_IRUSR, S_IXUSR
#include <sys/stat.h>	// chmod()

#include <sstream>
#include <string>
#include <list>

static void create_volume_results_name(const int vid, const char* const wrk_dir, std::string& name) {
	name = wrk_dir;
	if (name[name.size() - 1] != '/') {
		name += '/';
	}
	name += "r_";
	std::ostringstream x;
	x << vid;
	name += x.str();
}

static void merge_results(const char* const output, const std::list<std::string>& files) {
	std::string out_tmp(output);
	out_tmp += ".tmp";
	std::ofstream out(out_tmp.c_str());
	std::list<std::string>::const_iterator a(files.begin());
	const std::list<std::string>::const_iterator end_a(files.end());
	for (; a != end_a; ++a) {
		std::ifstream in(a->c_str());
		out << in.rdbuf();
		if (!in) {
			std::cerr << "Error reading from " << *a << "\n";
			exit(1);
		} else if (!out) {
			std::cerr << "Error writing to " << out_tmp << "\n";
			exit(1);
		}
	}
	out.close();
	if (rename(out_tmp.c_str(), output) == -1) {
		std::cerr << "Could not rename concatenated output file: " << out_tmp << "\n";
		exit(1);
	}
}

static void grid_start(std::string prog, const options_t &options, const int i, const int num_vols) {
	// create grid script, have grid run it
	options_t new_options = options;
	new_options.job_index = i;
	new_options.grid_options = NULL;
	new_options.grid_options_split = NULL;
	new_options.num_vols = num_vols;
	std::string name("m2pw.");
	if (i == -1) {
		name += "split";
	} else {
		std::ostringstream x;
		x << i;
		name += x.str();
	}
	std::string script_file(options.wrk_dir);
	script_file += "/" + name + ".sh";
	std::ofstream out;
	// as we might not have write permission, delete it
	unlink(script_file.c_str());
	open_fstream(out, script_file.c_str(), std::ios::out);
	out << "#!/bin/bash\nulimit -c 0\n" << prog << make_options(&new_options) << "\n";
	if (!out) {
		std::cerr << "Error writing to " << script_file << "\n";
		exit(1);
	}
	close_fstream(out);
	chmod(script_file.c_str(), S_IRUSR | S_IXUSR);
	std::string cmd(i == -1 && options.grid_options_split ? options.grid_options_split : options.grid_options);
	cmd += " " + name + " " + script_file;
	assert(system(cmd.c_str()) == 0);
}

// pass by value so we can modify list to easily avoid checking previously
// found results files
static void wait_for_files(std::list<std::string> results) {
	const std::list<std::string>::const_iterator end_a(results.end());
	while (!results.empty()) {
		sleep(60);
		std::list<std::string>::iterator a(results.begin());
		while (a != end_a) {
			if (access(a->c_str(), F_OK) == 0) {
				a = results.erase(a);
			} else {
				++a;
			}
		}
	}
}

int main(int argc, char* argv[]) {
	options_t options;
	if (parse_arguments(argc, argv, &options)) {
		print_usage(argv[0]);
		return 1;
	}
	char vol_idx_file_name[PATH_MAX];
	generate_idx_file_name(options.wrk_dir, vol_idx_file_name);
	volume_names_t* vn;
	if (options.num_vols == 0) {		// just do the split
		split_raw_dataset(options.reads, options.wrk_dir);
		return 0;
	} else if (options.grid_options == NULL && options.grid_options_split == NULL) {
	} else if (options.num_vols == -1) {	// spin off split and wait
		// check to see if split is already done
		vn = load_volume_names(vol_idx_file_name, 0);
		if (vn->num_vols == 0) {
			grid_start(argv[0], options, -1, 0);
			std::list<std::string> split_results;
			split_results.push_back(vol_idx_file_name);
			wait_for_files(split_results);
		}
		delete_volume_names_t(vn);
	}
	int num_vols = options.num_vols != -1 ? options.num_vols : split_raw_dataset(options.reads, options.wrk_dir);
	vn = load_volume_names(vol_idx_file_name, num_vols);
	r_assert(num_vols == vn->num_vols);
	if (options.job_index != -1) {
		std::string volume_results_name;
		create_volume_results_name(options.job_index, options.wrk_dir, volume_results_name);
		process_one_volume(&options, options.job_index, vn->num_vols, volume_results_name, vn);
		return 0;
	}
	std::cout << vol_idx_file_name << "\n";
	std::list<std::string> results;
	for (int i(0); i < vn->num_vols; ++i) {
		if (options.reads_to_correct) {
			const char *volume_name(get_vol_name(vn, i));
			volume_t *v(load_volume_header(volume_name));
			if (v->start_read_id > options.reads_to_correct) {
				break;
			}
			delete_volume_t(v);
		}
		std::string volume_results_name;
		create_volume_results_name(i, options.wrk_dir, volume_results_name);
		results.push_back(volume_results_name);
		if (access(volume_results_name.c_str(), F_OK) == 0) {
			LOG(stderr, "volume %d has been finished", i);
		} else if (options.grid_options == NULL) {
			process_one_volume(&options, i, vn->num_vols, volume_results_name, vn);
		} else {
			grid_start(argv[0], options, i, vn->num_vols);
		}
	}
	delete_volume_names_t(vn);
	if (options.grid_options != NULL) {
		wait_for_files(results);
	}
	merge_results(options.output, results);
	return 0;
}
