#include "mecat_correction.h"

#include <cassert>	// assert()
#include <limits>	// numeric_limits::max()
#include <stdint.h>	// uint8_t

#include "MECAT_AlnGraphBoost.H"
#include "reads_correction_aux.h"	// GAP_CHAR

#define FMAT 1
#define FDEL 2
#define FINS 4
#define UNDS 8

// returns type of coverage present

static inline uint8_t identify_one_consensus_item(const CnsTableItem& cns_item) {
	const int cov((cns_item.mat_cnt + cns_item.ins_cnt) * 0.8);
	uint8_t ident;
	if (cns_item.mat_cnt >= cov) {		// coverage is 80% or more matches
		ident = FMAT;
	} else if (cns_item.ins_cnt >= cov) {	// coverage is 80% or more inserts
		ident = FINS;
	} else {				// neither of the above
		ident = UNDS;
	}
	if (cns_item.del_cnt * 2 >= cov) {	// deletes are 40% or more than the coverage
		ident |= FDEL;
	}
	return ident;
}

struct CompareOverlapByOverlapSize {
	bool operator()(const Overlap& a, const Overlap& b) {
		const int64_t ovlp_a(std::max(a.qend - a.qoff, a.send - a.soff));
		const int64_t ovlp_b(std::max(b.qend - b.qoff, b.send - b.soff));
		return ovlp_b < ovlp_a;
	}
};

static void meap_add_one_aln(const std::string& qaln, const std::string& saln, int64_t start_soff, std::vector<CnsTableItem>& cns_table) {
	assert(qaln.size() == saln.size());
	const int64_t aln_size(qaln.size());
	for (int64_t i(0); i < aln_size;) {
		const char q(qaln[i]);
		const char s(saln[i]);
		if (q == GAP_CHAR && s == GAP_CHAR) {	// skip
			++i;
		} else if (q == s) {		// match
			++cns_table[start_soff].mat_cnt;
			cns_table[start_soff].base = s;
			++start_soff;
			++i;
		} else if (q == GAP_CHAR) {		// insert
			++cns_table[start_soff].ins_cnt;
			++start_soff;
			++i;
		} else {			// delete
			assert(s == GAP_CHAR);
			for (++i; i < aln_size && saln[i] == GAP_CHAR; ++i) { }
			++cns_table[start_soff - 1].del_cnt;
		}
	}
}

static void meap_cns_one_indel(const int sb, const int se, CnsAlns& cns_vec, const int min_cov, std::string& aux_qstr, std::string& aux_tstr, std::string& cns) {
	ns_meap_cns::AlnGraphBoost ag(se - sb + 1);
	int sb_out;
	for (std::vector<CnsAln>::iterator a(cns_vec.begin()); a != cns_vec.end(); ++a) {
		if (a->retrieve_aln_subseqs(sb, se, aux_qstr, aux_tstr, sb_out)) {
			ag.addAln(aux_qstr, aux_tstr, sb_out - sb + 1);
		}
	}		
	ag.mergeNodes();
	ag.consensus(min_cov * 0.4, cns);
}

static void meap_consensus_one_segment(const std::vector<CnsTableItem>& cns_list, const int cns_list_size, std::vector<uint8_t>& cns_id_vec, const int start_soff, CnsAlns& cns_vec, std::string& aux_qstr, std::string& aux_tstr, std::string& target) {
	if (cns_id_vec.size() < static_cast<size_t>(cns_list_size)) {
		cns_id_vec.resize(cns_list_size);
	}
	// get types of coverage
	for (int i(0); i < cns_list_size; ++i) {
		cns_id_vec[i] = identify_one_consensus_item(cns_list[start_soff + i]);
	}
	std::string cns;
	target.clear();
	const uint8_t unds_or_fdel(UNDS | FDEL);	// questionable coverage types
	int i(0); 
	// advance to matching coverage
	for (; i < cns_list_size && !(cns_id_vec[i] & FMAT); ++i) { }
	while (i < cns_list_size) {
		target.push_back(cns_list[start_soff + i].base);
		const int start(i);
		// advance to next matching coverage
		for (++i; i < cns_list_size && !(cns_id_vec[i] & FMAT); ++i) { }
		int need_refinement(0);
		// check to see if anything in-between has questionable coverage
		for (int k(start); k != i; ++k) {
			if (cns_id_vec[k] & unds_or_fdel) {
				need_refinement = 1;
				break;
			}
		}
		if (need_refinement) {
			meap_cns_one_indel(start_soff + start, start_soff + i, cns_vec, cns_list[start_soff + start].mat_cnt + cns_list[start_soff + start].ins_cnt, aux_qstr, aux_tstr, cns);
			// trim first and last as they have good coverage
			if (cns.size() > 2) {
				target.append(cns.data() + 1, cns.size() - 2);
			}
		}
	}
}

struct CmpMappingRangeBySoff {
	bool operator()(const MappingRange& a, const MappingRange& b) {
		if (a.start != b.start) {
			return a.start < b.start;
		} else {
			return b.end < a.end;
		}
	}
};

static void get_effective_ranges(std::vector<MappingRange>& mranges, std::vector<MappingRange>& eranges, const int read_size, const int min_size) {
	eranges.clear();
	if (mranges.size() == 0) {
		return;
	}
	// see if any ranges are effectively the whole read
	std::vector<MappingRange>::const_iterator a(mranges.begin());
	const std::vector<MappingRange>::const_iterator end_a(mranges.end());
	for (; a != end_a; ++a) {
		if (a->start <= 500 && read_size - a->end <= 500) {
			eranges.push_back(MappingRange(0, read_size));
			return;
		}
	}
	std::sort(mranges.begin(), mranges.end(), CmpMappingRangeBySoff());
	const int nr(mranges.size());
	// -1 so we can use > instead of >=
	const int min_size_95(ceil(min_size * 0.95) - 1);
	int i(0), left(mranges[0].start);
	for (;;) {
		const int right(mranges[i].end);
		// include ranges we completely overlap
		for (++i; i < nr && mranges[i].end <= right; ++i) { }
		if (i == nr) {
			if (right - left > min_size_95) {
				eranges.push_back(MappingRange(left, right));
			}
			break;
		}
		// if the next non-contained range overlaps by 1k or more,
		// treat as part of this range and keep extending
		if (right - mranges[i].start < 1000) {
			// truncate current range at overlap start
			const int eff_right(std::min(right, mranges[i].start));
			if (eff_right - left > min_size_95) {
				eranges.push_back(MappingRange(left, eff_right));
			}
			// put some space between the effective ranges
			left = std::max(right, mranges[i].start);
		}
	}
}

// breaks output sequence into chunks of no more than MaxSeqSize (if needed);
// split chunks will have an OvlpSize overlap

static void output_cns_result(std::vector<CnsResult>& cns_results, CnsResult& cr, const int64_t beg, const int64_t end, const std::string& cns_seq) {
	const size_t MaxSeqSize = 60000;
	const size_t OvlpSize = 10000;
	// BlkSize must be >= OvlpSize
	const size_t BlkSize = MaxSeqSize - OvlpSize - 1000;
	
	const size_t size = cns_seq.size();
	if (size <= MaxSeqSize) {
		cr.range[0] = beg;
		cr.range[1] = end;
		cr.seq = cns_seq;
		cns_results.push_back(cr);
	} else {
		const size_t cutoff = size - OvlpSize - 1000;
		size_t L = 0, R;
		do {
			R = L + BlkSize;
			if (R >= cutoff) {
				R = size;
			}
			cr.range[0] = L + beg;
			cr.range[1] = R < size && R + beg < static_cast<size_t>(end) ? R + beg : end;
			cr.seq = cns_seq.substr(L, R - L);
			cns_results.push_back(cr);
			L = R - OvlpSize;
		} while (R < size);
	}
}

static inline bool check_ovlp_mapping_range(const int qb, const int qe, const int qs, const int sb, const int se, const int ss, double ratio) {
	return qe - qb >= qs * ratio || se - sb >= ss * ratio;
}

// look for areas of high coverage of about min_size or more,
// improve them and stick on the results pile

static void consensus_worker(const std::vector<CnsTableItem>& cns_table, std::vector<uint8_t>& id_list, CnsAlns& cns_vec, std::string& aux_qstr, std::string& aux_tstr, const std::vector<MappingRange>& eranges, const int min_cov, const int min_size, const int read_id, std::vector<CnsResult>& cns_results) {
	CnsResult cns_result;
	cns_result.id = read_id;
	const int64_t min_size_95(ceil(0.95 * min_size));
	std::string cns_seq;
	std::vector<MappingRange>::const_iterator a(eranges.begin());
	const std::vector<MappingRange>::const_iterator end_a(eranges.end());
	for (; a != end_a; ++a) {
		const int end_i(a->end);
		for (int64_t i(a->start); i < end_i;) {
			// find start of next high coverage area
			for (; i < end_i && cns_table[i].mat_cnt + cns_table[i].ins_cnt < min_cov; ++i) { }
			const int64_t start(i);
			// find end of high coverage area
			for (++i; i < end_i && cns_table[i].mat_cnt + cns_table[i].ins_cnt >= min_cov; ++i) { }
			if (i - start >= min_size_95) {
				meap_consensus_one_segment(cns_table, i - start, id_list, start, cns_vec, aux_qstr, aux_tstr, cns_seq);
				if (cns_seq.size() >= static_cast<size_t>(min_size)) {
					output_cns_result(cns_results, cns_result, start, i, cns_seq);
				}
			}
		}
	}
}

static void decode_and_append_sequence(std::string& s, const std::string& seq, int64_t i, const int64_t end_i) {
	s.reserve(s.size() + end_i - i);
	for (; i < end_i; ++i) {
		s += "ACGT"[static_cast<int>(seq[i])];
	}
}

// same as consensus_worker, but produces entire read as one entry;
// uncorrected sections are just copied as is;

static void consensus_worker_one_read(const std::vector<CnsTableItem>& cns_table, std::vector<uint8_t>& id_list, CnsAlns& cns_vec, std::string& aux_qstr, std::string& aux_tstr, const std::vector<MappingRange>& eranges, const int min_cov, const int min_size, const int read_id, const std::string& tstr, std::vector<CnsResult>& cns_results) {
	CnsResult cns_result;
	cns_result.id = read_id;
	const int64_t min_size_95(std::max<int>(ceil(0.95 * min_size), 1) - 1);
	std::string cns_seq;
	std::vector<MappingRange>::const_iterator a(eranges.begin());
	const std::vector<MappingRange>::const_iterator end_a(eranges.end());
	std::vector<MappingRange>::const_iterator last_a(eranges.end());
	for (; a != end_a; last_a = a++) {
		if (last_a != end_a) {		// add in-between range to cns_result
			decode_and_append_sequence(cns_result.seq, tstr, last_a->end, a->start);
		} else if (a->start > 0) {	// add beginning of read
			decode_and_append_sequence(cns_result.seq, tstr, 0, a->start);
		}
		const int begin_i(a->start - 1);
		const int end_i(a->end);
		for (int64_t i(begin_i); i < end_i;) {
			// find start of next high coverage area
			const int64_t last_end(i != begin_i ? i : a->start);
			for (++i; i < end_i && cns_table[i].mat_cnt + cns_table[i].ins_cnt < min_cov; ++i) { }
			// add low coverage area as-is
			decode_and_append_sequence(cns_result.seq, tstr, last_end, i);
			if (i == end_i) {
				break;
			}
			const int64_t start(i);
			// find end of high coverage area
			for (++i; i < end_i && cns_table[i].mat_cnt + cns_table[i].ins_cnt >= min_cov; ++i) { }
			if (i - start > min_size_95) {
				meap_consensus_one_segment(cns_table, i - start, id_list, start, cns_vec, aux_qstr, aux_tstr, cns_seq);
				if (cns_seq.size() >= static_cast<size_t>(min_size)) {
					// add corrected sequence
					cns_result.seq += cns_seq;
					continue;
				}
			}
			// add uncorrected sequence
			decode_and_append_sequence(cns_result.seq, tstr, start, i);
		}
	}
	// add end of read
	if (last_a != end_a) {
		decode_and_append_sequence(cns_result.seq, tstr, last_a->end, tstr.size());
	}
	cns_result.range[0] = 0;
	cns_result.range[1] = cns_result.seq.size();
	cns_results.push_back(cns_result);
}

void consensus_one_read_m4_pacbio(ConsensusThreadData& ctd, ConsensusPerThreadData &pctd, const int64_t read_id, const int64_t sid, const int64_t eid) {
	PackedDB& reads(ctd.reads);
	ExtensionCandidate* overlaps((ExtensionCandidate*)pctd.candidates);
	DiffRunningData& drd(pctd.drd);
	M5Record& m5(pctd.m5);
	CnsAlns& cns_vec(pctd.cns_alns);
	std::vector<CnsResult>& cns_results(pctd.cns_results);
	const int64_t read_size(overlaps[read_id].ssize);
	std::string& qstr(pctd.query);
	std::string& tstr(pctd.target);
	reads.GetSequence(read_id, 1, tstr);
	std::string& nqstr(pctd.qaln);
	std::string& ntstr(pctd.saln);
	const int min_align_size(ctd.rco.min_align_size);
	const double error_rate(ctd.rco.error_rate);
	const int max_added(60);
	std::vector<CnsTableItem>& cns_table(pctd.cns_table);
	cns_table.assign(read_size, CnsTableItem());	// reset table
	cns_vec.clear();
	const int64_t L(sid);
	const int64_t R(eid - sid <= max_added ? eid : L + max_added);
	if (eid - sid > max_added) {	// only use largest max_added overlaps
		std::sort(overlaps + sid, overlaps + eid, CompareOverlapByOverlapSize());
	}
	for (int64_t i(L); i < R; ++i) {
		Overlap& ovlp(overlaps[i]);
		reads.GetSequence(ovlp.qid, ovlp.qdir == FWD, qstr);
		const int64_t qext(ovlp.qdir == FWD ? ovlp.qext : ovlp.qsize - 1 - ovlp.qext);
		const int r(GetAlignment(qstr, qext, tstr, ovlp.sext, drd, m5, error_rate, min_align_size));
		if (r) {
			normalize_gaps(m5.qaln, m5.saln, nqstr, ntstr, 1);
			meap_add_one_aln(nqstr, ntstr, m5.soff, cns_table);
			cns_vec.add_aln(m5.soff, m5.send, nqstr, ntstr);
		}
	}
	std::vector<MappingRange> mranges, eranges;
	cns_vec.get_mapping_ranges(mranges);
	get_effective_ranges(mranges, eranges, read_size, ctd.rco.min_size);
	consensus_worker(cns_table, pctd.id_list, cns_vec, nqstr, ntstr, eranges, ctd.rco.min_cov, ctd.rco.min_size, read_id, cns_results);
}

void consensus_one_read_m4_nanopore(ConsensusThreadData& ctd, ConsensusPerThreadData &pctd, const int64_t read_id, const int64_t sid, const int64_t eid) {
	PackedDB& reads = ctd.reads;
	ExtensionCandidate* overlaps = (ExtensionCandidate*)pctd.candidates;
	DiffRunningData& drd = pctd.drd;
	M5Record& m5 = pctd.m5;
	CnsAlns& cns_vec = pctd.cns_alns;
	std::vector<CnsResult>& cns_results = pctd.cns_results;
	const int64_t read_size = overlaps[read_id].ssize;
	std::string& qstr = pctd.query;
	std::string& tstr = pctd.target;
	reads.GetSequence(read_id, 1, tstr);
	std::string& nqstr = pctd.qaln;
	std::string& ntstr = pctd.saln;
	const int min_align_size = ctd.rco.min_align_size;
	const double error_rate(ctd.rco.error_rate);
	const double min_mapping_ratio = ctd.rco.min_mapping_ratio - 0.02;

	int64_t L, R;
	if (eid - sid <= MAX_CNS_OVLPS)
	{
		L = sid;
		R = eid;
	}
	else
	{
		L = sid;
		R = L + MAX_CNS_OVLPS;    
		std::sort(overlaps + sid, overlaps + eid, CompareOverlapByOverlapSize());
	}

	std::vector<CnsTableItem>& cns_table = pctd.cns_table;
	cns_table.assign(read_size, CnsTableItem());	// reset table
	cns_vec.clear();
	for (int64_t i(L); i < R; ++i) {
		Overlap& ovlp(overlaps[i]);
		reads.GetSequence(ovlp.qid, ovlp.qdir == FWD, qstr);
		const int64_t qext(ovlp.qdir == FWD ? ovlp.qext : ovlp.qsize - 1 - qext);
		const int r(GetAlignment(qstr, qext, tstr, ovlp.sext, drd, m5, error_rate, min_align_size));
		if (r && check_ovlp_mapping_range(m5.qoff, m5.qend, ovlp.qsize, m5.soff, m5.send, ovlp.ssize, min_mapping_ratio)) {
			normalize_gaps(m5.qaln, m5.saln, nqstr, ntstr, 1);
			meap_add_one_aln(nqstr, ntstr, m5.soff, cns_table);
			cns_vec.add_aln(m5.soff, m5.send, nqstr, ntstr);
		}
	}
	
	std::vector<MappingRange> mranges, eranges;
	eranges.push_back(MappingRange(0, read_size));

	consensus_worker(cns_table, pctd.id_list, cns_vec, nqstr, ntstr, eranges, ctd.rco.min_cov, ctd.rco.min_size, read_id,  cns_results);
}

// limit coverage - if there's not yet sufficient coverage, add coverage and return true;
// once coverage gets high enough, stop adding coverage and return false

static inline int check_cov_stats(std::vector<uint8_t>& cov_stats, const int soff, const int send) {
	if (cov_stats.size() < static_cast<size_t>(send)) {
		cov_stats.resize(send, 0);
	}
	int n(0);
	for (int i(soff); i < send; ++i) {
		if (cov_stats[i] >= 20) {	// max coverage
			++n;
		}
	}
	if (send - soff >= n + 200) {
		for (int i(soff); i < send; ++i) {
			// don't let small redundant region cause overflow
			if (cov_stats[i] != std::numeric_limits<uint8_t>::max()) {
				++cov_stats[i];
			}
		}
		return 1;
	}
	return 0;
}

void consensus_one_read_can_pacbio(ConsensusThreadData& ctd, ConsensusPerThreadData& pctd, const int64_t read_id, const int64_t sid, int64_t eid) {
	const PackedDB& reads(ctd.reads);
	ExtensionCandidateCompressed* candidates((ExtensionCandidateCompressed*)pctd.candidates);
	DiffRunningData& drd(pctd.drd);
	M5Record& m5(pctd.m5);
	CnsAlns& cns_vec(pctd.cns_alns);
	std::vector<CnsResult>& cns_results(pctd.cns_results);
	const int64_t ssize(reads.read_size(read_id));
	std::string& qstr(pctd.query);
	std::string& tstr(pctd.target);
	reads.GetSequence(read_id, 1, tstr);
	std::string& nqstr(pctd.qaln);
	std::string& ntstr(pctd.saln);
	const int min_align_size(ctd.rco.min_align_size);
	const double error_rate(ctd.rco.error_rate);
	const double min_mapping_ratio(ctd.rco.min_mapping_ratio - 0.02);
	int num_added(0);
	const int max_added(60);
	eid = std::min(eid, sid + 200);		// max of 200 extents
	std::vector<CnsTableItem>& cns_table(pctd.cns_table);
	cns_table.assign(ssize, CnsTableItem());		// reset table
	cns_vec.clear();
	std::set<int64_t> used_ids;
	std::vector<uint8_t>& id_list(pctd.id_list);		// used to be called cov_stats
	id_list.clear();
	for (int64_t i(sid); i < eid && num_added < max_added; ++i) {
		const ExtensionCandidateCompressed& ec(candidates[i]);
		if (used_ids.find(ec.qid) != used_ids.end()) {
			continue;
		}
		const int64_t qsize(reads.read_size(ec.qid));
		reads.GetSequence(ec.qid, ec.qdir() == FWD, qstr);
		const int64_t qext(ec.qdir() == FWD ? ec.qext() : qsize - 1 - ec.qext());
		const int r(GetAlignment(qstr, qext, tstr, ec.sext, drd, m5, error_rate, min_align_size));
		if (r && check_ovlp_mapping_range(m5.qoff, m5.qend, qsize, m5.soff, m5.send, ssize, min_mapping_ratio)) {
			if (check_cov_stats(id_list, m5.soff, m5.send)) {
				++num_added;
				used_ids.insert(ec.qid);
				normalize_gaps(m5.qaln, m5.saln, nqstr, ntstr, 1);
				meap_add_one_aln(nqstr, ntstr, m5.soff, cns_table);
				cns_vec.add_aln(m5.soff, m5.send, nqstr, ntstr);
			}
		}
	}
	std::vector<MappingRange> mranges, eranges;
	cns_vec.get_mapping_ranges(mranges);
	get_effective_ranges(mranges, eranges, ssize, ctd.rco.min_size);
	if (ctd.rco.full_reads) {
		consensus_worker_one_read(cns_table, pctd.id_list, cns_vec, nqstr, ntstr, eranges, ctd.rco.min_cov, ctd.rco.min_size, read_id, tstr, cns_results);
	} else {
		consensus_worker(cns_table, pctd.id_list, cns_vec, nqstr, ntstr, eranges, ctd.rco.min_cov, ctd.rco.min_size, read_id, cns_results);
	}
}

void consensus_one_read_can_nanopore(ConsensusThreadData& ctd, ConsensusPerThreadData &pctd, const int64_t read_id, const int64_t sid, const int64_t eid) {
	PackedDB& reads(ctd.reads);
	ExtensionCandidateCompressed* candidates((ExtensionCandidateCompressed*)pctd.candidates);
	DiffRunningData& drd(pctd.drd);
	M5Record& m5(pctd.m5);
	CnsAlns& cns_vec(pctd.cns_alns);
	std::vector<CnsResult>& cns_results(pctd.cns_results);
	const int64_t ssize(reads.read_size(read_id));
	std::string& qstr(pctd.query);
	std::string& tstr(pctd.target);
	reads.GetSequence(read_id, 1, tstr);
	std::string& nqstr(pctd.qaln);
	std::string& ntstr(pctd.saln);
	const int min_align_size(ctd.rco.min_align_size);
	const double error_rate(ctd.rco.error_rate);
	const double min_mapping_ratio(ctd.rco.min_mapping_ratio - 0.02);
	int num_added(0);
	int num_ext(0);
	const int max_ext(200);
	std::vector<CnsTableItem>& cns_table(pctd.cns_table);
	cns_table.assign(ssize, CnsTableItem());		// reset table
	cns_vec.clear();
	std::set<int> used_ids;
	std::vector<uint8_t>& id_list(pctd.id_list);	// used to be called cov_stats
	id_list.clear();
	for (int64_t i(sid); i < eid && num_added < MAX_CNS_OVLPS && num_ext < max_ext; ++i) {
		++num_ext;
		const ExtensionCandidateCompressed& ec(candidates[i]);
		if (used_ids.find(ec.qid) != used_ids.end()) {
			continue;
		}
		const int64_t qsize(reads.read_size(ec.qid));
		reads.GetSequence(ec.qid, ec.qdir() == FWD, qstr);
		const int64_t qext(ec.qdir() == FWD ? ec.qext() : qsize - 1 - ec.qext());
		const int r(GetAlignment(qstr, qext, tstr, ec.sext, drd, m5, error_rate, min_align_size));
		if (r && check_ovlp_mapping_range(m5.qoff, m5.qend, qsize, m5.soff, m5.send, ssize, min_mapping_ratio)) {
			if (check_cov_stats(id_list, m5.soff, m5.send)) {
				++num_added;
				used_ids.insert(ec.qid);
				normalize_gaps(m5.qaln, m5.saln, nqstr, ntstr, 1);
				meap_add_one_aln(nqstr, ntstr, m5.soff, cns_table);
				cns_vec.add_aln(m5.soff, m5.send, nqstr, ntstr);
			}
		}
	}
	std::vector<MappingRange> mranges, eranges;
	eranges.push_back(MappingRange(0, ssize));
	consensus_worker(cns_table, pctd.id_list, cns_vec, nqstr, ntstr, eranges, ctd.rco.min_cov, ctd.rco.min_size, read_id,  cns_results);
}
