#ifndef PACKED_DB_H
#define PACKED_DB_H

#include <fstream>	// ifstream
#include <string>	// string
#include <stdint.h>	// uint8_t
#include <unistd.h>	// off_t
#include <vector>	// vector<>

#include "../common/alignment.h"	// ExtensionCandidateCompressed
#include "../common/sequence.h"		// Sequence

class PackedDB {
    private:
	struct SeqIndex {
		off_t file_offset;
		int64_t memory_offset, size;
		explicit SeqIndex() { }
		explicit SeqIndex(const off_t i, const int64_t j, const int64_t k) : file_offset(i), memory_offset(j), size(k) { }
	};
    public:
	explicit PackedDB() : pac_(0), db_size(0), max_db_size(0), max_read_size_(0) { }
	~PackedDB() {
		delete[] pac_;
	}
	// returns number of reads
	static size_t convert_fasta_to_db(const std::string& fasta, const std::string& output_prefix, int64_t min_size);
	static void read_sizes(const std::string& output_prefix, std::vector<int64_t>& sizes);
	// only call one of load_fasta_db and open_db exactly once
	void load_fasta_db(const char* fasta);
	// opens data file, reads in index file
	void open_db(const std::string& filename, int64_t memory_footprint);
	// returns number of candidates that can be processed
	int64_t load_reads(const ExtensionCandidateCompressed* ec_list, int64_t nec);
	void GetSequence(const int64_t id, const int forward, std::string& seq) const {
		const SeqIndex &si(seq_idx[id]);
		seq.resize(si.size);
		if (forward) {
			const int64_t offset(si.memory_offset);
			for (int64_t i(0); i < si.size; ++i) {
				seq[i] = get_char(offset + i);
			}
		} else {
			const int64_t offset(si.memory_offset + si.size - 1);
			for (int64_t i(0); i < si.size; ++i) {
				seq[i] = 3 - get_char(offset - i);
			}
		}
	}
	int64_t num_reads() const {
		return seq_idx.size();
	}
	int64_t read_size(const int64_t read_id) const {
		return seq_idx[read_id].size;
	}
	int64_t max_read_size() const {
		return max_read_size_;
	}
    private:
	static void set_char(std::vector<uint8_t>& p, const int64_t idx, const uint8_t c) {
		p[idx >> 2] |= c << ((~idx & 3) << 1);
	}
	void set_char(const int64_t idx, const uint8_t c) {
		// use ~x instead of 3 - x for speed, since we have to & 3 anyway
		pac_[idx >> 2] |= c << ((~idx & 3) << 1);
	}
	uint8_t get_char(const int64_t idx) const {
		return pac_[idx >> 2] >> ((~idx & 3) << 1) & 3;
	}
	void add_one_seq(const Sequence& seq);
	void destroy() {
		delete[] pac_;
		seq_idx.clear();
		max_db_size = db_size = 0;
	}
    private:
	uint8_t* pac_;
	int64_t db_size, max_db_size, max_read_size_;
	std::vector<SeqIndex> seq_idx;
	std::ifstream pstream;
};

#endif // PACKED_DB_H
