#ifndef _OVERLAPS_STORE_H
#define _OVERLAPS_STORE_H

#include <fstream>
#include <limits>
#include <string>

#include "../common/defs.h"
#include "../common/pod_darr.h"

template <class T>
class PartitionResultsWriter
{
public:
	typedef void (*file_name_generator)(const char* prefix, const idx_t id, std::string& name);
	
public:
	PartitionResultsWriter(const int max_num_files)
	{
		file_is_open = 0;
		num_open_files = 0;
		MaxNumFiles = max_num_files;
		results = new PODArray<T>[MaxNumFiles];
		files = new std::ofstream[MaxNumFiles];
		file_names = new std::string[MaxNumFiles];
		min_seq_ids = new idx_t[MaxNumFiles];
		max_seq_ids = new idx_t[MaxNumFiles];
		for (int i = 0; i < MaxNumFiles; ++i) results[i].reserve(kStoreSize);
	}
	
	~PartitionResultsWriter()
	{
		delete[] results;
		delete[] files;
		delete[] file_names;
		delete[] min_seq_ids;
		delete[] max_seq_ids;
	}
	
	void OpenFiles(const idx_t sfid, const idx_t efid, const std::string& prefix, file_name_generator fng)
	{
		if (file_is_open) CloseFiles();
		
		const int nf = efid - sfid;
        if (nf == 0) return;
        for (int i = 0; i < nf; ++i)
        {
            fng(prefix.data(), i + sfid, file_names[i]);
            open_fstream(files[i], file_names[i].c_str(), std::ios::binary);
            min_seq_ids[i] = std::numeric_limits<index_t>::max();
            max_seq_ids[i] = std::numeric_limits<index_t>::min();
			results[i].clear();
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
        file_is_open = false;
        num_open_files = 0;
	}
	
	void WriteOneResult(const int fid, const idx_t seq_id, const T& r)
	{
		if (fid >= num_open_files) {
			std::cout << "fid = " << fid
					  << ", num_open_files = "
					  << num_open_files
					  << "\n";
		}
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
	int MaxNumFiles;
	static const int kStoreSize = 500000;
	PODArray<T>* results;
	bool file_is_open;
	int num_open_files;
	std::ofstream* files;
	std::string* file_names;
	idx_t* min_seq_ids;
	idx_t* max_seq_ids;
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
