#ifndef _READS_CORRECTION_AUX_H
#define _READS_CORRECTION_AUX_H

#include <cassert>	// assert()
#include <fstream>	// ifstream, ofstream
#include <iostream>	// ostream
#include <mutex>	// lock_guard<>, mutex
#include <string>	// string
#include <stdint.h>	// int64_t, uint8_t
#include <unistd.h>	// F_OK, access(), rename()
#include <vector>	// vector<>

#include "../common/alignment.h"	// ExtensionCandidate, ExtensionCandidateCompressed
#include "dw.h"		// DiffRunningData, ERROR(), LOG(), M5Record
#include "options.h"	// ReadsCorrectionOptions
#include "packed_db.h"	// PackedDB

#define MAX_CNS_OVLPS 100
#define GAP_CHAR '-'

// 1k seems to work a bit better than 10k - perhaps less time waiting for
// another thread to finish writing?
#define MAX_CNS_RESULTS 1000

struct CnsResult {
	int64_t id, range[2];
	std::string seq;
};

struct CnsTableItem {
	char base;
	uint8_t mat_cnt, ins_cnt, del_cnt;
	explicit CnsTableItem() : base('N'), mat_cnt(0), ins_cnt(0), del_cnt(0) { }
};

struct MappingRange {
	int start, end;
	explicit MappingRange(const int s, const int e) : start(s), end(e) { }
};

class CnsAln : public MappingRange {
    public:
	explicit CnsAln(const int i, const int j, const std::string& q, const std::string& s) : MappingRange(i, j), aln_idx_(0), qaln_(s), saln_(q) { }
	~CnsAln() { }
	int retrieve_aln_subseqs(const int sb, const int se, std::string& qstr, std::string& tstr, int& sb_out) {
		const int aln_size(saln_.size() - 1);
		if (se <= start || sb >= end || aln_idx_ == aln_size) {
			return 0;
		}
		sb_out = std::max(start, sb);
		while (start < sb && aln_idx_ != aln_size) {
			if (saln_[++aln_idx_] != GAP_CHAR) {
				++start;
			}
		}
		// should we test for start < sb, and return 0 if so,
		// rather than returning the last basepair of the alignment?
		const int aln_start(aln_idx_);
		while (start < se && aln_idx_ != aln_size) {
			if (saln_[++aln_idx_] != GAP_CHAR) {
				++start;
			}
		}
		// this looks like it could return the same basepair twice -
		// once at the end of a call, once at the start of the next
		const int aln_length(aln_idx_ - aln_start + 1);
		qstr.assign(qaln_, aln_start, aln_length);
		tstr.assign(saln_, aln_start, aln_length);
		return 1;
	}
    private:
	int aln_idx_;
	// despite never changing, we can't make these const because this class
	// is used in a vector which uses default construction and copy mechanics
	std::string qaln_, saln_;
};

class CnsAlns {
    public:
	explicit CnsAlns() { }
	~CnsAlns() { }
	void clear() {
		cns_alns_.clear();
	}
	size_t num_alns() const {
		return cns_alns_.size();
	}
	std::vector<CnsAln>::iterator begin() {
		return cns_alns_.begin();
	}
	std::vector<CnsAln>::const_iterator end() const {
		return cns_alns_.end();
	}
	void add_aln(const int soff, const int send, const std::string& qstr, const std::string& tstr) {
		assert(qstr.size() == tstr.size());
		cns_alns_.push_back(CnsAln(soff, send, qstr, tstr));
	}
	void get_mapping_ranges(std::vector<MappingRange>& ranges) const {
		ranges.assign(cns_alns_.begin(), cns_alns_.end());
	}
    private:
	std::vector<CnsAln> cns_alns_;
};

struct CmpExtensionCandidateBySidAndScore {
	bool operator()(const ExtensionCandidate& a, const ExtensionCandidate& b) {
		if (a.sid != b.sid) {			// primary sort
			return a.sid < b.sid;
		} else if (a.score != b.score) {	// secondary sort
			return b.score < a.score;	// process best ones first
		} else if (a.qid != b.qid) {
			return a.qid < b.qid;
		} else if (a.qext != b.qext) {
			return a.qext < b.qext;
		} else if (a.sext != b.sext) {		// tertiary sort
			return a.sext < b.sext;		// make sorting consistent
		} else if (a.soff != b.soff) {
			return a.soff < b.soff;
		} else if (a.qoff != b.qoff) {
			return a.qoff < b.qoff;
		} else if (a.send != b.send) {
			return a.send < b.send;
		} else if (a.qend != b.qend) {
			return a.qend < b.qend;
		} else {
			return a.qdir < b.qdir;
		}
		// sdir, qsize, ssize are all the same by this point
	}
};

struct CmpExtensionCandidateCompressedBySidAndScore {
	bool operator()(const ExtensionCandidateCompressed& a, const ExtensionCandidateCompressed& b) {
		if (a.sid != b.sid) {			// primary sort
			return a.sid < b.sid;		// for splitting up in allocate_ecs()
		} else if (a.score != b.score) {	// secondary sort
			return b.score < a.score;	// process best ones first
		} else if (a.qid != b.qid) {
			return a.qid < b.qid;
		} else if (a.qext() != b.qext()) {
			return a.qext() < b.qext();
		} else if (a.sext != b.sext) {		// tertiary sort
			return a.sext < b.sext;		// make sorting consistent
		} else {
			return a.qdir() < b.qdir();
		}
	}
};

// like CmpExtensionCandidateCompressedBySidAndScore, but use a new ordering
// (don't store actual vector in this class, as it gets copied a lot by sort())

class CmpExtensionCandidateCompressedNewOrder {
    public:
	explicit CmpExtensionCandidateCompressedNewOrder(const std::vector<int64_t>& order) : order_(order) { }
	bool operator()(const ExtensionCandidateCompressed& a, const ExtensionCandidateCompressed& b) {
		if (a.sid != b.sid) {			// primary sort
							// for splitting up in allocate_ecs()
			return order_[a.sid] < order_[b.sid];
		} else if (a.score != b.score) {	// secondary sort
			return b.score < a.score;	// process best ones first
		} else if (a.qid != b.qid) {		// group ties by qid
			return a.qid < b.qid;		// (doesn't matter if new or old)
		} else if (a.qext() != b.qext()) {
			return a.qext() < b.qext();
		} else if (a.sext != b.sext) {		// tertiary sort
			return a.sext < b.sext;		// make sorting consistent
		} else {
			return a.qdir() < b.qdir();
		}
	}
    private:
	const std::vector<int64_t>& order_;	// [read_id] = new_order
};

class ConsensusPerThreadData {
    public:
	// this is ExtensionCandidate for m4 runs, ExtensionCandidateCompressed
	// for candidate runs (to reduce memory usage)
	void* candidates;
	// num_candidates, candidates initialized by allocate_ecs()
	// next_candidate initialized by ConsensusThreadData::restart()
	int64_t num_candidates, next_candidate;
	DiffRunningData drd;
	M5Record m5;
	CnsAlns cns_alns;
	std::vector<CnsTableItem> cns_table;
	std::vector<uint8_t> id_list;
	std::vector<CnsResult> cns_results;
	std::string query, target, qaln, saln;
    public:
	explicit ConsensusPerThreadData(const double error_rate, const int64_t max_read_size) : drd(error_rate, max_read_size) {
		// we'll definitely be seeing at least this much use,
		// so might as well preallocate
		cns_results.reserve(MAX_CNS_RESULTS);
	}
	~ConsensusPerThreadData() { }
};

class ConsensusThreadData {
    public:
	ReadsCorrectionOptions& rco;
	PackedDB& reads;
	std::ostream& out;
	std::mutex out_lock;
	int64_t ec_offset;
	std::vector<ConsensusPerThreadData> data;
    public:
	ConsensusThreadData(ReadsCorrectionOptions& prco, PackedDB& r, std::ostream& output, const std::string& input_file_name) :
		rco(prco),
		reads(r),
		out(output),
		ec_offset(0),
		data(rco.num_threads, ConsensusPerThreadData(rco.error_rate, reads.max_read_size())),
		last_thread_id_(-1),
		num_threads_written_(0),
		done_file_(input_file_name + ".done"),
		ckpt_file_(input_file_name + ".ckpt"),
		ckpt_file_tmp_(ckpt_file_ + ".tmp") { }
	~ConsensusThreadData() { }
	int get_thread_id() {
		std::lock_guard<std::mutex> lock(id_lock_);
		return ++last_thread_id_;
	}
	void reset_threads() {
		last_thread_id_ = -1;
		std::vector<ConsensusPerThreadData>::iterator a(data.begin());
		const std::vector<ConsensusPerThreadData>::const_iterator end_a(data.end());
		for (; a != end_a; ++a) {
			a->next_candidate = 0;
		}
	}
	void write_buffer(const int tid, const int64_t i) {
		ConsensusPerThreadData& pdata(data[tid]);
		{						// scoping for lock_guard
			std::vector<CnsResult>::const_iterator a(pdata.cns_results.begin());
			const std::vector<CnsResult>::const_iterator end_a(pdata.cns_results.end());
			std::lock_guard<std::mutex> lock(out_lock);
			for (; a != end_a; ++a) {
				out << ">" << a->id << "_" << a->range[0] << "_" << a->range[1] << "_" << a->seq.size() << "\n" << a->seq << "\n";
				if (!out) {
					ERROR("Error writing output");
				}
			}
			pdata.next_candidate = i;
			if (++num_threads_written_ >= rco.num_threads) {
				checkpoint();
				num_threads_written_ = 0;
			}
		}
		pdata.cns_results.clear();
	}
	int restart(off_t& output_pos) {
		std::vector<ConsensusPerThreadData>::iterator a(data.begin());
		const std::vector<ConsensusPerThreadData>::const_iterator end_a(data.end());
		if (access(ckpt_file_.c_str(), F_OK) != 0) {		// fresh start
			for (; a != end_a; ++a) {
				a->next_candidate = 0;
			}
			return 0;
		}
		std::ifstream ckpt_in(ckpt_file_.c_str());
		if (!ckpt_in) {
			ERROR("Restart failed: could not open checkpoint file: %s", ckpt_file_.c_str());
		}
		ckpt_in >> rco.job_index >> output_pos >> ec_offset;
		if (!ckpt_in) {
			ERROR("Restart failed: could not read checkpoint file: %s", ckpt_file_.c_str());
		}
		for (; a != end_a; ++a) {
			ckpt_in >> a->next_candidate;
			if (!ckpt_in) {
				ERROR("Restart failed: could not read checkpoint file: %s", ckpt_file_.c_str());
			}
		}
		return 1;
	}
    private:
	void checkpoint() {
		std::ofstream ckpt_out(ckpt_file_tmp_.c_str());
		if (!ckpt_out) {
			LOG(stderr, "Checkpoint failed: couldn't open %s", ckpt_file_tmp_.c_str());
			return;
		}
		out.flush();
		ckpt_out << rco.job_index << " " << off_t(out.tellp()) << " " << ec_offset << "\n";
		if (!ckpt_out) {
			LOG(stderr, "Checkpoint failed: write failed: %s", ckpt_file_tmp_.c_str());
			return;
		}
		std::vector<ConsensusPerThreadData>::iterator a(data.begin());
		const std::vector<ConsensusPerThreadData>::const_iterator end_a(data.end());
		for (; a != end_a; ++a) {
			ckpt_out << a->next_candidate << "\n";
			if (!ckpt_out) {
				LOG(stderr, "Checkpoint failed: write failed: %s", ckpt_file_tmp_.c_str());
				return;
			}
		}
		ckpt_out.close();
		if (rename(ckpt_file_tmp_.c_str(), ckpt_file_.c_str()) == -1) {
			LOG(stderr, "Checkpoint failed: rename failed: %s", ckpt_file_.c_str());
		}
	}
    private:
	std::mutex id_lock_;
	int last_thread_id_, num_threads_written_;
	std::string done_file_, ckpt_file_, ckpt_file_tmp_;
};

void normalize_gaps(const std::string& qstr, const std::string& tstr, std::string& qnorm, std::string& tnorm, int push);

void allocate_ecs(ConsensusThreadData &data, ExtensionCandidate* ec_list, int64_t nec);
void allocate_ecs(ConsensusThreadData &data, ExtensionCandidateCompressed* ec_list, int64_t nec);

#endif // _READS_CORRECTION_AUX_H
