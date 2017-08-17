#ifndef _OVERLAPS_STORE_H
#define _OVERLAPS_STORE_H

#include <fstream>
#include <limits>
#include <string>
#include <vector>
#include <unistd.h>
#include <iostream>

#include "../common/defs.h"
#include "../common/pod_darr.h"

template <class T>
class PartitionResultsWriter
{
public:
	typedef void (*file_name_generator)(const char* prefix, const idx_t id, std::string& name);
	
public:
	PartitionResultsWriter()
	{
		file_is_open = false;
		num_open_files = false;

                kStoreSize = 0;
		// leave room for stdin, stdout, stderr, a few others
		kNumFiles = sysconf(_SC_OPEN_MAX) - 10;
		files = NULL;
		results = NULL;
	}
	~PartitionResultsWriter()
	{
		CloseFiles();
	}
	void OpenFiles(const idx_t sfid, const idx_t efid, const std::string& prefix, file_name_generator fng)
	{
		CloseFiles();
		const int nf = efid - sfid;
		if (nf == 0) return;
		// allocate about a gb of memory as buffer
		kStoreSize = 1073741824 / sizeof(ExtensionCandidate) / nf;
		file_names.assign(nf, "");
		min_seq_ids.assign(nf, std::numeric_limits<idx_t>::max());
		max_seq_ids.assign(nf, std::numeric_limits<idx_t>::min());
		files = new std::ofstream[nf];
		results = new PODArray<T>[nf];
		for (int i = 0; i < nf; ++i)
		{
			fng(prefix.data(), i + sfid, file_names[i]);
			open_fstream(files[i], file_names[i].c_str(), std::ios::binary);
			results[i].reserve(kStoreSize);
		}
		num_open_files = nf;
		file_is_open = true;
	}
	void CloseFiles()
	{
		if (!file_is_open) return;
		for (int i = 0; i < num_open_files; ++i)
		{
			if (results[i].size())
			{
				char* buf = (char*)results[i].data();
				std::streamsize s = sizeof(T) * results[i].size();
				files[i].write(buf, s);
			}
		}
		for (int i = 0; i < num_open_files; ++i) close_fstream(files[i]);
		file_names.clear();
		min_seq_ids.clear();
		max_seq_ids.clear();
		delete[] files;
		delete[] results;
		files = NULL;
		results = NULL;
		kStoreSize = 0;
		file_is_open = false;
		num_open_files = 0;
	}
	void WriteOneResult(const int fid, const idx_t seq_id, const T& r)
	{
		r_assert(fid < num_open_files);
		min_seq_ids[fid] = std::min(min_seq_ids[fid], seq_id);
        max_seq_ids[fid] = std::max(max_seq_ids[fid], seq_id);
		results[fid].push_back(r);
		if (results[fid].size() == kStoreSize)
		{
			char* buf = (char*)results[fid].data();
			std::streamsize s = sizeof(T) * results[fid].size();
			files[fid].write(buf, s);
			results[fid].clear();
		}
	}
	
public:
	int kNumFiles;
	int kStoreSize;
	
	PODArray<T> *results;	// can't use vector<>, causes memory corruption
	bool file_is_open;
	int num_open_files;
	std::ofstream *files;	// can't use vector<>, non-copyable
	std::vector<std::string> file_names;
	std::vector<idx_t> min_seq_ids;
	std::vector<idx_t> max_seq_ids;
};

template <class T>
T* load_partition_data(const char* path, idx_t& num_results)
{
	std::ifstream in;
	open_fstream(in, path, std::ios::binary);
	in.seekg(0, std::ios::end);
	std::streampos fs = in.tellg();
	in.seekg(0, std::ios::beg);
	num_results = fs / sizeof(T);
	T* arr = new T[num_results];
	in.read((char*)arr, fs);
	close_fstream(in);
	return arr;
}

#endif // _OVERLAPS_STORE_H
