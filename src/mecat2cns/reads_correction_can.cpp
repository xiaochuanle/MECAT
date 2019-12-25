#include "reads_correction_can.h"

#include <cassert>	// assert()
#include <cstring>
#include <limits>	// numeric_limits<>
#include <sstream>
#include <string>	// string

#include "../common/defs.h"	// TECH_PACBIO
#include "MECAT_AlnGraphBoost.H"
#include "mecat_correction.h"
#include "overlaps_partition.h"
#include "overlaps_store.h"
#include "options.h"		// ReadsCorrectionOptions

static void* reads_correction_func_can(void* const arg) {
	ConsensusThreadData& data(*(static_cast<ConsensusThreadData*>(arg)));
	const int tid(data.get_thread_id());
	ConsensusPerThreadData& pdata(data.data[tid]);
	const ExtensionCandidateCompressed* const candidates((const ExtensionCandidateCompressed*)pdata.candidates);
	int64_t i(pdata.next_candidate);
	if (data.rco.tech == TECH_PACBIO) {
		while (i != pdata.num_candidates) {
			const int64_t start(i);
			const int64_t sid(candidates[start].sid);
			for (++i; i != pdata.num_candidates && candidates[i].sid == sid; ++i) { }
			// still have to check this here, as it's not always checked earlier
			if (i - start < data.rco.min_cov) {
				continue;
			}
			consensus_one_read_can_pacbio(data, pdata, sid, start, i);
			if (pdata.cns_results.size() >= MAX_CNS_RESULTS) {
				data.write_buffer(tid, i);
			}
		}
	} else {
		while (i != pdata.num_candidates) {
			const int64_t start(i);
			const int64_t sid(candidates[start].sid);
			for (++i; i != pdata.num_candidates && candidates[i].sid == sid; ++i) { }
			if (i - start < data.rco.min_cov) {
				continue;
			}
			consensus_one_read_can_nanopore(data, pdata, sid, start, i);
			if (pdata.cns_results.size() >= MAX_CNS_RESULTS) {
				data.write_buffer(tid, i);
			}
		}
	}
	data.write_buffer(tid, i);
	return NULL;
}

struct CmpExtensionCandidateCompressedBySid {
	bool operator()(const ExtensionCandidateCompressed& a, const ExtensionCandidateCompressed& b) {
		return a.sid < b.sid;
	}
};

class EC_Index {	// offset into ec_list (and number of ecs) for each read id
    public:
	int64_t offset, count;
	explicit EC_Index() : offset(0), count(0) { }
	explicit EC_Index(int64_t i, int64_t j) : offset(i), count(j) { }
	~EC_Index() { }
};

// reorder list so reads are more concentrated and we can process more
// ecs per pass, with limited read space; also filters list to exclude
// candidates of reads with low coverage; we filter out alignments of
// low-coverage reads, so nec may be reduced

static void reorder_candidates(ExtensionCandidateCompressed* const ec_list, int64_t& nec, const int64_t reads_to_correct, const int64_t num_reads, const int min_cov) {
	// allow us to easily access a given sid's aligns
	std::sort(ec_list, ec_list + nec, CmpExtensionCandidateCompressedBySid());
	int64_t total_ec(0);
	// count will be zero for any read id beyond reads_to_correct
	std::vector<EC_Index> index(reads_to_correct);
	// index existing list by sid
	for (int64_t i(0); i != nec;) {
		// we are guaranteed sid is < reads_to_correct by partitioning
		const int64_t sid(ec_list[i].sid);
		const int64_t start(i);
		for (++i; i != nec && ec_list[i].sid == sid; ++i) { }
		// make sure we have enough coverage
		const int64_t count(i - start);
		if (count >= min_cov) {
			index[sid] = EC_Index(start, count);
			total_ec += count;
		}
	}
	// generate the new read order
	std::vector<char> used(num_reads, 0);
	std::vector<int64_t> new_order;
	new_order.reserve(num_reads);	// possible overestimate, but whatever
	int64_t next_unused(0);
	size_t next_search(0);
	for (;;) {
		// skip over used reads, reads with no alignments
		for (; next_unused != reads_to_correct && (used[next_unused] || index[next_unused].count == 0); ++next_unused) { }
		if (next_unused == reads_to_correct) {
			break;
		}
		used[next_unused] = 1;
		new_order.push_back(next_unused);
		// add all reads aligned to, and aligned to those, and so on
		for (; next_search != new_order.size(); ++next_search) {
			const int64_t sid(new_order[next_search]);
			// make sure read has index entry
			if (sid < reads_to_correct) {
				const EC_Index& a(index[sid]);
				int64_t i(a.offset);
				const int64_t end_i(i + a.count);
				for (; i != end_i; ++i) {
					const int64_t qid(ec_list[i].qid);
					if (!used[qid]) {
						used[qid] = 1;
						new_order.push_back(qid);
					}
				}
			}
		}
	}
	// re-sort ec_list with new order (note that coverage-excluded
	// ec's will sort last, so we also change nec to exclude them);
	// sort reads not used (no new_order entry) at end
	std::vector<int64_t> rid_to_order(num_reads, std::numeric_limits<int64_t>::max());
	for (size_t i(0); i != new_order.size(); ++i) {
		rid_to_order[new_order[i]] = i;
	}
	CmpExtensionCandidateCompressedNewOrder new_rid_order_cmp(rid_to_order);
	std::sort(ec_list, ec_list + nec, new_rid_order_cmp);
	nec = total_ec;
}

// load and sort partition data, assign to threads, start threads
static void consensus_one_partition_can(const char* const m4_file_name, ConsensusThreadData& data) {
	int64_t nec;
	ExtensionCandidateCompressed* ec_list(load_partition_data<ExtensionCandidateCompressed>(m4_file_name, nec));
	// if we're memory limited spend some cpu time to speed up passes
	// (~16s for ~300s speedup in test case)
	if (data.rco.read_buffer_size) {
		reorder_candidates(ec_list, nec, data.rco.reads_to_correct ? data.rco.reads_to_correct : data.reads.num_reads(), data.reads.num_reads(), data.rco.min_cov);
	} else {
		std::sort(ec_list, ec_list + nec, CmpExtensionCandidateCompressedBySidAndScore());
	}
	pthread_t thread_ids[data.rco.num_threads];
	while (data.ec_offset != nec) {
		// see how many candidates we can run, given
		// how much read sequence we can load into memory
		// (unless we're not limited, in which case load 'em all)
		const int64_t ecs(data.rco.read_buffer_size ? data.reads.load_reads(ec_list + data.ec_offset, nec - data.ec_offset) : nec);
		allocate_ecs(data, ec_list + data.ec_offset, ecs);
		for (int i(0); i != data.rco.num_threads; ++i) {
			pthread_create(&thread_ids[i], NULL, reads_correction_func_can, static_cast<void*>(&data));
		}
		for (int i(0); i != data.rco.num_threads; ++i) {
			pthread_join(thread_ids[i], NULL);
		}
		// update offset after threading to avoid checkpointing with new value
		data.ec_offset += ecs;
		data.reset_threads();
	}
	delete[] ec_list;
}

// set up output file and thread data, handle restart if needed
static int reads_correction_can_p(ReadsCorrectionOptions& rco, std::vector<std::string> &partition_file_vec, PackedDB& reads) {
	const std::string& p(partition_file_vec[rco.job_index]);
	std::ostringstream os;
	os << rco.corrected_reads << "." << rco.job_index;
	const std::string results_file(os.str());
	const std::string results_file_tmp(results_file + ".tmp");
	std::ofstream out;
	ConsensusThreadData data(rco, reads, out, results_file.c_str());
	off_t output_pos;
	if (data.restart(output_pos)) {					// is a restart
		open_fstream(out, results_file_tmp.c_str(), std::ios::out | std::ios::in);
		if (!out.seekp(output_pos)) {
			ERROR("Restart failed: output file seek failed: %s", results_file_tmp.c_str());
		}
	} else {
		open_fstream(out, results_file_tmp.c_str(), std::ios::out);
	}
	char process_info[1024];
	sprintf(process_info, "processing %s", p.c_str());
	DynamicTimer dtimer(process_info);
	consensus_one_partition_can(p.c_str(), data);
	close_fstream(out);
	assert(rename(results_file_tmp.c_str(), results_file.c_str()) == 0);
	return 0;
}

// load reads, partition file index
int reads_correction_can(ReadsCorrectionOptions& rco) {
	std::string idx_file_name;
	generate_partition_index_file_name(rco.m4, idx_file_name);
	std::vector<std::string> partition_file_vec;
	load_partition_files_info(idx_file_name.c_str(), partition_file_vec);
	PackedDB reads;
	reads.open_db("fasta.db", rco.read_buffer_size);
	if (rco.job_index != -1) {
		return reads_correction_can_p(rco, partition_file_vec, reads);
	} else {
		rco.job_index = 0;	// may get changed by data.restart() below
		std::string tmp_file(rco.corrected_reads);
		tmp_file += ".tmp";
		std::ofstream out;
		ConsensusThreadData data(rco, reads, out, rco.corrected_reads);
		off_t output_pos;
		if (data.restart(output_pos)) {				// is a restart
			open_fstream(out, rco.corrected_reads, std::ios::out | std::ios::in);
			out.seekp(output_pos);
			if (!out.seekp(output_pos)) {
				ERROR("Restart failed: output file seek failed: %s", rco.corrected_reads);
			}
		} else {
			open_fstream(out, rco.corrected_reads, std::ios::out);
		}
		char process_info[1024];
		const int job_end(partition_file_vec.size());
		for (; rco.job_index != job_end; ++rco.job_index) {
			const std::string& p(partition_file_vec[rco.job_index]);
			sprintf(process_info, "processing %s", p.c_str());
			DynamicTimer dtimer(process_info);
			consensus_one_partition_can(p.c_str(), data);
		}
		close_fstream(out);
		assert(rename(tmp_file.c_str(), rco.corrected_reads) == 0);
	}
	return 0;
}
