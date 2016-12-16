#include "pw_options.h"
#include "pw_impl.h"
#include "../common/split_database.h"

#include <cstdio>
#include <fstream>
#include <unistd.h>

#include <sstream>
#include <string>

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
merge_results(const char* output, const char* wrk_dir, const int num_volumes)
{
	string vrn;
	for (int i = 0; i < num_volumes; ++i)
	{
		create_volume_results_name_finished(i, wrk_dir, vrn);
		ostringstream cmd;
		if (i == 0) cmd << "cat " << vrn << " >" << output;
		else cmd << "cat " << vrn << " >> " << output;
		assert(system(cmd.str().c_str()) == 0);
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
	
	int num_vols = split_raw_dataset(options.reads, options.wrk_dir);
	
	char vol_idx_file_name[1024];
	generate_idx_file_name(options.wrk_dir, vol_idx_file_name);
	cout << vol_idx_file_name << "\n";
	volume_names_t* vn = load_volume_names(vol_idx_file_name, 0);
	r_assert(num_vols == vn->num_vols);
	for (int i = 0; i < vn->num_vols; ++i)
	{
		string volume_results_name_finished;
		create_volume_results_name_finished(i, options.wrk_dir, volume_results_name_finished);
		if (access(volume_results_name_finished.c_str(), F_OK) == 0) 
		{
			LOG(stderr, "volume %d has been finished\n", i);
			continue;
		}
		string volume_results_name_working;
		create_volume_results_name_working(i, options.wrk_dir, volume_results_name_working);
		ofstream out;
		open_fstream(out, volume_results_name_working.c_str(), ios::out);
		process_one_volume(&options, i, vn->num_vols, vn, &out);
		close_fstream(out);
		assert(rename(volume_results_name_working.c_str(), volume_results_name_finished.c_str()) == 0);
	}
	vn = delete_volume_names_t(vn);
	
	merge_results(options.output, options.wrk_dir, num_vols);
}
