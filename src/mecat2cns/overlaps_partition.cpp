#include "overlaps_partition.h"

#include <cassert>	// assert()
#include <fstream>
#include <limits>
#include <set>
#include <sys/stat.h>	// stat(), struct stat
#include <vector>
#include <utility>	// make_pair(), pair<>
#include <math.h>	// ceil()

#include "overlaps_store.h"
#include "reads_correction_aux.h"
#include "packed_db.h"	// PackedDB

inline static bool query_is_contained(const M4Record& m4, const double min_cov_ratio) {
	return m4qend(m4) - m4qoff(m4) >= m4qsize(m4) * min_cov_ratio;
}

inline static bool subject_is_contained(const M4Record& m4, const double min_cov_ratio) {
	return m4send(m4) - m4soff(m4) >= m4ssize(m4) * min_cov_ratio;
}

inline static bool check_m4record_mapping_range(const M4Record& m4, const double min_cov_ratio) {
	return query_is_contained(m4, min_cov_ratio) || subject_is_contained(m4, min_cov_ratio);
}

static int64_t get_qualified_m4record_counts(const char* const m4_file_name, const double min_cov_ratio) {
	std::ifstream in;
	open_fstream(in, m4_file_name, std::ios::in);
	int64_t num_reads(-1), num_records(0), num_qualified_records(0);
	M4Record m4;
	m4qext(m4) = m4sext(m4) = INVALID_IDX;
	while (in >> m4) {
		if (m4qext(m4) == INVALID_IDX || m4sext(m4) == INVALID_IDX) {
			ERROR("no gapped start position is provided, please make sure that you have run \'meap_pairwise\' with option \'-g 1\'");
		}
		++num_records;
		if (check_m4record_mapping_range(m4, min_cov_ratio)) {
			++num_qualified_records;
		}
		num_reads = std::max(num_reads, m4qid(m4));
		num_reads = std::max(num_reads, m4sid(m4));
	}
	close_fstream(in);
	LOG(stderr, "there are %ld overlaps, %ld are qualified.", num_records, num_qualified_records);
	return num_reads + 1;
}

// not in use at the moment
#if 0
static void get_repeat_reads(const char* const m4_file_name, const double min_cov_ratio, const int64_t num_reads, std::set<int64_t>& repeat_reads) {
	const int max_contained(100);
	// used to increment to a max of max_contained
	char count_table[max_contained + 1];
	for (int i(0); i < max_contained; ++i) {
		count_table[i] = i + 1;
	}
	count_table[max_contained] = max_contained;
	std::vector<char> counts(num_reads, 0);
	std::ifstream in;
	open_fstream(in, m4_file_name, std::ios::in);
	M4Record m4;
	m4qext(m4) = m4sext(m4) = INVALID_IDX;
	while (in >> m4) {
		if (query_is_contained(m4, min_cov_ratio)) {
			const int64_t qid = m4qid(m4);
			// increments count up to max_contained
			counts[qid] = count_table[static_cast<int>(counts[qid])];
		}
		if (subject_is_contained(m4, min_cov_ratio)) {
			const int64_t sid = m4sid(m4);
			// increments count up to max_contained
			counts[sid] = count_table[static_cast<int>(counts[sid])];
		}
	}
	close_fstream(in);
	for (int64_t i(0); i < num_reads; ++i) {
		if (counts[i] == max_contained) {
			std::cerr << "repeat read " << i << "\n";
			repeat_reads.insert(i);
		}
	}
	LOG(stderr, "number of repeat reads: %lu", repeat_reads.size());
}
#endif

void generate_partition_index_file_name(const std::string& input_file_name, std::string& ret) {
	ret = input_file_name + ".partition_files";
}

static void generate_partition_file_name(const std::string& input_file_name, const int part, std::string& ret) {
	std::ostringstream os;
	os << part;
	ret = input_file_name + ".part" + os.str();
}

void partition_candidates(const std::string& input, const std::string& pac_prefix, const size_t file_size, const int max_files_per_batch, const int64_t num_reads) {
	DynamicTimer dtimer(__func__);
	struct stat buf;
	if (stat(input.c_str(), &buf) == -1) {
		ERROR("Could not get file size: %s", input.c_str());
	}
	// each candidate line takes on average 44 characters, but go with 32;
	// each one produces two candidates (forward and reverse)
	const int num_files(ceil(double(buf.st_size / 32 * 2) * sizeof(ExtensionCandidateCompressed) / file_size));
	// separate them by read id as we don't know how many candidates each
	// read is part of; we're assuming an even distribution on average
	const int64_t reads_per_file((num_reads + num_files - 1) / num_files);
	std::vector<int64_t> read_sizes;
	PackedDB::read_sizes(pac_prefix, read_sizes);
	PartitionResultsWriter<ExtensionCandidateCompressed> prw(max_files_per_batch);
	int i(0);
	off_t input_pos;
	// if restarting, both i and input_pos can be changed
	int is_restart(prw.restart(input, generate_partition_file_name, "partition.done", i, input_pos));
	std::string idx_file_name;
	generate_partition_index_file_name(input, idx_file_name);
	std::ofstream idx_file;
	open_fstream(idx_file, idx_file_name.c_str(), std::ios::out);
	ExtensionCandidate ec;
	ExtensionCandidateCompressed nec;
	// and here we go through the input file to write num_files,
	// being limited by the number of open output files we can have
	for (; i < num_files; i += prw.kNumFiles) {
		const int sfid(i);
		const int efid(std::min(sfid + prw.kNumFiles, num_files));
		const int nf(efid - sfid);
		const int64_t L(sfid * reads_per_file);
		const int64_t R(std::min(efid * reads_per_file, num_reads));
		std::ifstream in;
		open_fstream(in, input.c_str(), std::ios::in);
		if (is_restart) {
			if (!in.seekg(input_pos)) {
				ERROR("Input seek failed while restoring checkpoint: %s", input.c_str());
			}
			is_restart = 0;
		} else {
			prw.OpenFiles(sfid, efid, input, generate_partition_file_name, "partition.done");
		}
		while (in >> ec) {
			// size has been set to zero for too-small reads
			if (!read_sizes[ec.sid] || !read_sizes[ec.qid]) {
				continue;
			}
			// make sure values match before tossing the values
			assert(ec.ssize == read_sizes[ec.sid] && ec.qsize == read_sizes[ec.qid]);
			// as variable sizes may be different, make sure the ones
			// we read in can be safely stored
			assert(ec.sid <= ExtensionCandidateCompressed::max_value);
			assert(ec.qid <= ExtensionCandidateCompressed::max_value);
			assert(ec.sext <= ExtensionCandidateCompressed::max_value);
			assert(ec.score <= ExtensionCandidateCompressed::max_value);
			// qext is one bit smaller than the others
			assert(ec.qext <= ExtensionCandidateCompressed::max_qext);
			if (L <= ec.sid && ec.sid < R) {
				nec.set(ec);
				if (prw.WriteOneResult((nec.sid - L) / reads_per_file, nec.sid, nec)) {
					prw.checkpoint(in.tellg());
				}
			}
			if (L <= ec.qid && ec.qid < R) {
				nec.set_swap(ec);
				if (prw.WriteOneResult((nec.sid - L) / reads_per_file, nec.sid, nec)) {
					prw.checkpoint(in.tellg());
				}
			}
		}
		for (int k(0); k < nf; ++k) {
			std::cerr << prw.file_names[k] << " contains " << prw.counts[k] << " overlaps\n";
			if (prw.counts[k] != 0) {
				idx_file << prw.file_names[k] << "\n";
			}
		}
		prw.CloseFiles();
	}
	close_fstream(idx_file);
	prw.finalize();
}

void partition_m4records(const char* const m4_file_name, const double min_cov_ratio, const size_t file_size, const int min_read_size, const int max_files_per_batch) {
	DynamicTimer dtimer(__func__);
	int64_t num_reads(get_qualified_m4record_counts(m4_file_name, min_cov_ratio));
	std::set<int64_t> repeat_reads;
	//get_repeat_reads(m4_file_name, min_cov_ratio, num_reads, repeat_reads);
	struct stat buf;
	if (stat(m4_file_name, &buf) == -1) {
		ERROR("Could not get file size: %s", m4_file_name);
	}
	const int num_files(ceil(double(buf.st_size / 32 * 2) * sizeof(ExtensionCandidate) / file_size));
	const int64_t reads_per_file((num_reads + num_files - 1) / num_files);
	std::string idx_file_name;
	generate_partition_index_file_name(m4_file_name, idx_file_name);
	std::ofstream idx_file;
	open_fstream(idx_file, idx_file_name.c_str(), std::ios::out);
	M4Record m4, nm4;
	ExtensionCandidate ec;
	PartitionResultsWriter<ExtensionCandidate> prw(max_files_per_batch);
	for (int i(0); i < num_files; i += prw.kNumFiles) {
		const int sfid(i);
		const int efid(std::min(sfid + prw.kNumFiles, num_files));
		const int nf(efid - sfid);
		const int64_t L(reads_per_file * sfid);
		const int64_t R(efid < num_files ? reads_per_file * efid : num_reads);
		std::ifstream in;
		open_fstream(in, m4_file_name, std::ios::in);
		prw.OpenFiles(sfid, efid, m4_file_name, generate_partition_file_name, "partition.done");
		while (in >> m4) {
			if (m4qsize(m4) < min_read_size || m4ssize(m4) < min_read_size) {
				continue;
			} else if (!check_m4record_mapping_range(m4, min_cov_ratio)) {
				continue;
			} else if (repeat_reads.find(m4qid(m4)) != repeat_reads.end() || repeat_reads.find(m4sid(m4)) != repeat_reads.end()) {
				continue;
			}
			if (m4qid(m4) >= L && m4qid(m4) < R) {
				normalize_m4record(m4, false, nm4);
				m4_to_candidate(nm4, ec);
				prw.WriteOneResult((m4qid(m4) - L) / reads_per_file, m4qid(m4), ec);
			}
			if (m4sid(m4) >= L && m4sid(m4) < R) {
				normalize_m4record(m4, true, nm4);
				m4_to_candidate(nm4, ec);
				prw.WriteOneResult((m4sid(m4) - L) / reads_per_file, m4sid(m4), ec);
			}
		}
		for (int k(0); k < nf; ++k) {
			std::cerr << prw.file_names[k] << " contains " << prw.counts[k] << " overlaps\n";
			if (prw.counts[k] != 0) {
				idx_file << prw.file_names[k] << "\n";
			}
		}
		prw.CloseFiles();
	}
	close_fstream(idx_file);
}

void load_partition_files_info(const char* const idx_file_name, std::vector<std::string>& file_info_vec) {
	file_info_vec.clear();
	std::ifstream in;
	open_fstream(in, idx_file_name, std::ios::in);
	std::string line;
	while (in >> line) {
		file_info_vec.push_back(line);
	}
	close_fstream(in);
}
