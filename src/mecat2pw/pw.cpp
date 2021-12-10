#include "pw_options.h"
#include "pw_impl.h"
#include "../common/split_database.h"

#include <cstdio>
#include <fstream>
#include <unistd.h>

#include <sstream>
#include <string>
#include <list>

using namespace std;

void
create_volume_results_name_working(int vid, const char* wrk_dir, string& name)
{
	name = wrk_dir;
	if (name[name.size() - 1] != '/') name += '/';
	ostringstream os;
	os << "r_" << vid << ".working";
	name += os.str();
}

void
create_volume_results_name_finished(int vid, const char* wrk_dir, string& name)
{
	name = wrk_dir;
	if (name[name.size() - 1] != '/') name += '/';
	ostringstream os;
	os << "r_" << vid;
	name += os.str();
}

void
merge_results(const char* output, const char* wrk_dir, const std::list<std::string>& files)
{
	string vrn;
	std::list<std::string>::const_iterator a = files.begin();
	const std::list<std::string>::const_iterator end_a = files.end();
	for (; a != end_a; ++a)
	{
		ostringstream cmd;
		if (a == files.begin()) cmd << "cat " << *a << " >" << output;
		else cmd << "cat " << *a << " >> " << output;
		assert(system(cmd.str().c_str()) == 0);
	}
}

void
grid_start(std::string prog, const options_t &options, const int i, const int num_vols, const std::string &volume_results_name_working, const std::string &volume_results_name_finished)
{
	// create grid command
	options_t new_options = options;
	new_options.job_index = i;
	new_options.grid_options = NULL;
	new_options.num_vols = num_vols;
	new_options.reads = volume_results_name_working.c_str();
	new_options.output = volume_results_name_finished.c_str();
	std::ostringstream cmd;
	cmd << "qsub -N m2pw." << i << " " << options.grid_options;
	cmd << " \"" << prog << make_options(&new_options) << "\"";
	assert(system(cmd.str().c_str()) == 0);
}

// pass by value so we can modify list to easily avoid checking previously
// found results files
void
wait_for_files(std::list<std::string> results)
{
	const std::list<std::string>::const_iterator end_a = results.end();
	while (!results.empty())
	{
		sleep(60);
		std::list<std::string>::iterator a = results.begin();
		while (a != end_a)
		{
			if (access(a->c_str(), F_OK) == 0)
			{
				a = results.erase(a);
			}
			else
			{
				++a;
			}
		}
	}
}

int main(int argc, char* argv[])
{
    options_t options;
    int r = parse_arguments(argc, argv, &options);
	if (r)
	{
		print_usage(argv[0]);
		return 1;
	}
	
	int num_vols = options.num_vols != -1 ? options.num_vols : split_raw_dataset(options.reads, options.wrk_dir);
	
	char vol_idx_file_name[1024];
	generate_idx_file_name(options.wrk_dir, vol_idx_file_name);
	volume_names_t* vn = load_volume_names(vol_idx_file_name, 0);
	r_assert(num_vols == vn->num_vols);
	if (options.job_index != -1)
	{
		int i = options.job_index;
		// repurposing these two for the wrapper
		std::string volume_results_name_working = options.reads;
		std::string volume_results_name_finished = options.output;
		std::ofstream out;
		open_fstream(out, volume_results_name_working.c_str(), std::ios::out);
		process_one_volume(&options, i, vn->num_vols, vn, &out);
		close_fstream(out);
		assert(rename(volume_results_name_working.c_str(), volume_results_name_finished.c_str()) == 0);
		return 1;
	}
	cout << vol_idx_file_name << "\n";
	std::list<std::string> results;
	for (int i = 0; i < vn->num_vols; ++i)
	{
		string volume_results_name_finished;
		create_volume_results_name_finished(i, options.wrk_dir, volume_results_name_finished);
		results.push_back(volume_results_name_finished);
		if (access(volume_results_name_finished.c_str(), F_OK) == 0) 
		{
			LOG(stderr, "volume %d has been finished\n", i);
			continue;
		}
		string volume_results_name_working;
		create_volume_results_name_working(i, options.wrk_dir, volume_results_name_working);
		if (options.grid_options == NULL)
		{
			ofstream out;
			open_fstream(out, volume_results_name_working.c_str(), ios::out);
			process_one_volume(&options, i, vn->num_vols, vn, &out);
			close_fstream(out);
			assert(rename(volume_results_name_working.c_str(), volume_results_name_finished.c_str()) == 0);
		}
		else
		{
			grid_start(argv[0], options, i, vn->num_vols, volume_results_name_working, volume_results_name_finished);
		}
	}
	vn = delete_volume_names_t(vn);

	if (options.grid_options != NULL)
	{
		wait_for_files(results);
	}
	
	merge_results(options.output, options.wrk_dir, results);
}
