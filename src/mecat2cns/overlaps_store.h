#ifndef _OVERLAPS_STORE_H
#define _OVERLAPS_STORE_H

#include <fstream>	// ifstream, ofstream, streampos, streamsize
#include <string>	// string
#include <stdint.h>	// int64_t
#include <time.h>	// time(), time_t
#include <unistd.h>	// _SC_OPEN_MAX, F_OK, access(), off_t, rename(), sysconf(), unlink()
#include <vector>	// vector<>

#include "../common/defs.h"	// ERROR(), LOG()
#include "../common/pod_darr.h"	// PODArray<>

template <class T> class PartitionResultsWriter {
    public:
	typedef void (*file_name_generator)(const std::string& prefix, int id, std::string& name);
    public:
	const int kNumFiles;	// effective open file limit
	int kStoreSize;
	int num_open_files;
	PODArray<T>* results;	// can't use vector<>, causes memory corruption
	std::ofstream* files;	// can't use vector<>, non-copyable
	std::vector<std::string> file_names;
	std::vector<int64_t> counts;
    public:
	// can't make kNumFiles static, as sysconf() is run-time only;
	// leave room for stdin, stdout, stderr, a few others
	explicit PartitionResultsWriter(const int max_files_per_batch) : kNumFiles(max_files_per_batch > 0 ? max_files_per_batch : sysconf(_SC_OPEN_MAX) - 10), kStoreSize(0), num_open_files(0), results(0), files(0) { }
	~PartitionResultsWriter() {
		CloseFiles();
	}
	void OpenFiles(const int sfid, const int efid, const std::string& prefix, file_name_generator fng, const std::string& done_file) {
		CloseFiles();
		if (efid <= sfid) {
			return;
		}
		done_file_ = done_file;
		ckpt_file_ = done_file_ + ".ckpt";
		ckpt_file_tmp_ = ckpt_file_ + ".tmp";
		num_open_files = efid - sfid;
		batch_start_ = sfid;
		allocate_data(prefix, fng, 0);
	}
	void CloseFiles() {
		if (num_open_files == 0) {
			return;
		}
		for (int i(0); i != num_open_files; ++i) {	// flush buffers
			if (results[i].size()) {
				write_buffer_to_disk(i);
			}
			files[i].close();
			// don't attempt to re-finish finished files
			if (!file_names[i].empty()) {
				const std::string tmp_file(file_names[i] + ".tmp");
				if (rename(tmp_file.c_str(), file_names[i].c_str()) == -1) {
					ERROR("Could not rename %s", tmp_file.c_str());
				}
			}
		}
		delete[] results;
		delete[] files;
		file_names.clear();
		counts.clear();
		kStoreSize = num_open_files = 0;
		results = 0;
		files = 0;
	}
	void finalize() {
		// touch done file
		std::ofstream done(done_file_.c_str());
		if (!done) {
			ERROR("Could not create done file %s", done_file_.c_str());
		}
		done.close();
		// remove checkpoint file
		unlink(std::string(done_file_ + ".ckpt").c_str());
	}
	int WriteOneResult(const int i, const int64_t seq_id, const T& r) {
		++counts[i];
		results[i].push_back(r);
		if (results[i].size() == kStoreSize) {
			write_buffer_to_disk(i);
			results[i].clear();
			if (time(0) >= next_checkpoint_time_) {
				return 1;
			}
		}
		return 0;
	}
	int restart(const std::string& prefix, file_name_generator fng, const std::string& done_file, int& sfid, off_t& input_pos) {
		CloseFiles();
		done_file_ = done_file;
		ckpt_file_ = done_file_ + ".ckpt";
		ckpt_file_tmp_ = ckpt_file_ + ".tmp";
		std::ifstream in(ckpt_file_.c_str());
		if (!in) {
			return 0;
		}
		in >> batch_start_ >> num_open_files >> input_pos;
		if (!in) {
			ERROR("Read error while restoring checkpoint from %s", ckpt_file_.c_str());
		}
		allocate_data(prefix, fng, 1);
		for (int i(0); i != num_open_files; ++i) {
			off_t file_pos;
			in >> file_pos >> counts[i];
			if (!in) {
				ERROR("Read error while restore checkpoint from %s (%d)", ckpt_file_.c_str(), i);
			}
			// don't need to restart finished files
			if (!file_names[i].empty() && !files[i].seekp(file_pos)) {
				ERROR("Seek failed while restoring checkpoint from %s (%d)", ckpt_file_.c_str(), i);
			}
		}
		sfid = batch_start_;
		return 1;
	}
	void checkpoint(const off_t input_pos) {
		std::ofstream out(ckpt_file_tmp_.c_str());
		if (!out) {
			LOG(stderr, "Checkpoint failed: couldn't open %s", ckpt_file_tmp_.c_str());
			return;
		}
		out << batch_start_ << " " << num_open_files << " " << input_pos << "\n";
		if (!out) {
			LOG(stderr, "Checkpoint failed: write failed: %s", ckpt_file_tmp_.c_str());
			return;
		}
		// flush buffers
		for (int i(0); i != num_open_files; ++i) {
			if (results[i].size()) {
				write_buffer_to_disk(i);
				results[i].clear();
			}
			files[i].flush();
			out << off_t(files[i].tellp()) << " " << counts[i] << "\n";
			if (!out) {
				LOG(stderr, "Checkpoint failed: write failed: %s", ckpt_file_tmp_.c_str());
				return;
			}
		}
		out.close();
		if (rename(ckpt_file_tmp_.c_str(), ckpt_file_.c_str()) == -1) {
			LOG(stderr, "Checkpoint failed: rename failed: %s", ckpt_file_.c_str());
		}
		next_checkpoint_time_ = time(0) + 300;
	}
	int64_t total_count() const {
		int64_t total(0);
		std::vector<int64_t>::const_iterator a(counts.begin());
		const std::vector<int64_t>::const_iterator end_a(counts.end());
		for (; a != end_a; ++a) {
			total += *a;
		}
		return total;
	}
    private:
	int batch_start_ ;
	time_t next_checkpoint_time_;
	std::string done_file_;
	std::string ckpt_file_;
	std::string ckpt_file_tmp_;
    private:
	void allocate_data(const std::string& prefix, file_name_generator fng, const int is_restart) {
		// allocate about a gb of memory as buffer, split among num_open_files
		kStoreSize = (1 << 30) / sizeof(T) / num_open_files;
		results = new PODArray<T>[num_open_files];
		file_names.assign(num_open_files, "");
		counts.assign(num_open_files, 0);
		files = new std::ofstream[num_open_files];
		for (int i(0); i != num_open_files; ++i) {
			fng(prefix, i + batch_start_, file_names[i]);
			const std::string tmp_file(file_names[i] + ".tmp");
			if (is_restart) {
				if (access(file_names[i].c_str(), F_OK) == 0) {
					// mark as already finished
					file_names[i].clear();
					// use /dev/null to prevent write errors
					files[i].open("/dev/null", std::ios::binary);
				} else {
					// don't truncate on restart
					files[i].open(tmp_file.c_str(), std::ios::binary | std::ios::in);
				}
			} else {
				files[i].open(tmp_file.c_str(), std::ios::binary);
			}
			if (!files[i]) {
				ERROR("Open failed on %s", tmp_file.c_str());
			}
			results[i].reserve(kStoreSize);
		}
		next_checkpoint_time_ = time(0) + 300;
	}
	void write_buffer_to_disk(const int i) {
		// can't use static_cast<>
		const char* const buf((char*)results[i].data());
		const std::streamsize s(sizeof(T) * results[i].size());
		if (!files[i].write(buf, s)) {
			ERROR("Error writing to %s", file_names[i].c_str());
		}
	}
};

template <class T> T* load_partition_data(const char* const path, int64_t& num_results) {
	std::ifstream in;
	open_fstream(in, path, std::ios::binary);
	in.seekg(0, std::ios::end);
	const std::streampos fs(in.tellg());
	in.seekg(0, std::ios::beg);
	num_results = fs / sizeof(T);
	T* const arr(new T[num_results]);
	in.read((char*)arr, fs);	// can't use static_cast<>
	if (!in) {
		ERROR("Error reading partition data: %s", path);
	}
	close_fstream(in);
	return arr;
}

#endif // _OVERLAPS_STORE_H
