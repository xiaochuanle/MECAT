#include "overlaps_partition.h"

#include <fstream>
#include <limits>
#include <set>
#include <vector>

#include <unistd.h>

#include "overlaps_store.h"
#include "reads_correction_aux.h"

using namespace std;

#define error_and_exit(msg) { std::cerr << msg << "\n"; abort(); }

inline bool
check_m4record_mapping_range(const M4Record& m4, const double min_cov_ratio)
{
    const index_t qm = m4qend(m4) - m4qoff(m4);
    const index_t qs = m4qsize(m4) * min_cov_ratio;
    const index_t sm = m4send(m4) - m4soff(m4);
    const index_t ss = m4ssize(m4) * min_cov_ratio;
    return qm >= qs || sm >= ss;
}

inline bool
query_is_contained(const M4Record& m4, const double min_cov_ratio)
{
	const index_t qm = m4qend(m4) - m4qoff(m4);
    const index_t qs = m4qsize(m4) * min_cov_ratio;
	return qm >= qs;
}

inline bool
subject_is_contained(const M4Record& m4, const double min_cov_ratio)
{
	const index_t sm = m4send(m4) - m4soff(m4);
    const index_t ss = m4ssize(m4) * min_cov_ratio;
	return sm >= ss;
}

void
get_qualified_m4record_counts(const char* m4_file_name, const double min_cov_ratio, index_t& num_qualified_records, index_t& num_reads)
{
    std::ifstream in;
    open_fstream(in, m4_file_name, std::ios::in);
    num_qualified_records = 0;
    num_reads = -1;
    index_t num_records = 0;
    M4Record m4;
	m4qext(m4) = m4sext(m4) = INVALID_IDX;
    while (in >> m4)
    {
		if (m4qext(m4) == INVALID_IDX || m4sext(m4) == INVALID_IDX)
		{
			ERROR("no gapped start position is provided, please make sure that you have run \'meap_pairwise\' with option \'-g 1\'");
		}
        ++num_records;
        if (check_m4record_mapping_range(m4, min_cov_ratio)) ++num_qualified_records;
        num_reads = std::max(num_reads, m4qid(m4));
        num_reads = std::max(num_reads, m4sid(m4));
    }
    close_fstream(in);

	LOG(stderr, "there are %d overlaps, %d are qualified.", (int)num_records, (int)num_qualified_records);
    ++num_reads;
}

void
get_repeat_reads(const char* m4_file_name, const double min_cov_ratio, const index_t num_reads, set<index_t>& repeat_reads)
{
	const char MaxContained = 100;
	char cnt_table[MaxContained + 1];
	for (int i = 0; i < MaxContained; ++i) cnt_table[i] = i + 1;
	cnt_table[MaxContained] = MaxContained;
	vector<char> cnts(num_reads, 0);
	std::ifstream in;
    open_fstream(in, m4_file_name, std::ios::in);
	M4Record m4;
	m4qext(m4) = m4sext(m4) = INVALID_IDX;
    while (in >> m4)
    {
		if (query_is_contained(m4, min_cov_ratio))
		{
			index_t qid = m4qid(m4);
			cnts[qid] = cnt_table[cnts[qid]];
		}
		if (subject_is_contained(m4, min_cov_ratio))
		{
			index_t sid = m4sid(m4);
			cnts[sid] = cnt_table[cnts[sid]];
		}
    }
    close_fstream(in);
	
	for(index_t i = 0; i < num_reads; ++i)
		if (cnts[i] >= MaxContained) 
		{
			cerr << "repeat read " << i << "\n";
			repeat_reads.insert(i);
		}

	LOG(stderr, "number of repeat reads: %d", (int)repeat_reads.size());
}

void 
generate_partition_index_file_name(const char* m4_file_name, std::string& ret)
{
    ret = m4_file_name;
    ret += ".partition_files";
}

void
generate_partition_file_name(const char* m4_file_name, const index_t part, std::string& ret)
{
    ret = m4_file_name;
    ret += ".part";
    std::ostringstream os;
    os << part;
    ret += os.str();
}

idx_t
get_num_reads(const char* candidates_file)
{
	ifstream in;
	open_fstream(in, candidates_file, ios::in);
	ExtensionCandidate ec;
	int max_id = -1;
	while (in >> ec)
	{
		max_id = std::max(ec.qid, max_id);
		max_id = std::max(ec.sid, max_id);
	}
	close_fstream(in);
	return max_id + 1;
}

void
normalise_candidate(ExtensionCandidate& src, ExtensionCandidate& dst, const bool subject_is_target)
{
	if (subject_is_target)
	{
		dst = src;
	}
	else
	{
		dst.qdir = src.sdir;
		dst.qid = src.sid;
		dst.qext = src.sext;
		dst.qsize = src.ssize;
		dst.sdir = src.qdir;
		dst.sid = src.qid;
		dst.sext = src.qext;
		dst.ssize = src.qsize;
		dst.score = src.score;
	}
	
	if (dst.sdir == REV)
	{
		dst.qdir = REVERSE_STRAND(dst.qdir);
		dst.sdir = REVERSE_STRAND(dst.sdir);
	}
}

int
fix_file_counts(int num_files) {
	if (num_files < 0) {
		num_files = sysconf(_SC_OPEN_MAX) - 10;
	}
	return num_files;
}

void
partition_candidates(const char* input, const idx_t batch_size, const int min_read_size, int num_files)
{
	DynamicTimer dt(__func__);
	
	num_files = fix_file_counts(num_files);
	const idx_t num_reads = get_num_reads(input);
	const idx_t num_batches = (num_reads + batch_size - 1) / batch_size;
	string idx_file_name;
	generate_partition_index_file_name(input, idx_file_name);
	ofstream idx_file;
	open_fstream(idx_file, idx_file_name.c_str(), ios::out);
	
	ExtensionCandidate ec, nec;
	PartitionResultsWriter<ExtensionCandidate> prw(num_files);
	for (idx_t i = 0; i < num_batches; i += num_files) {
		const idx_t sfid = i;
		const idx_t efid = min(sfid + num_files, num_batches);
		const int nf = efid - sfid;
		const idx_t Lid = batch_size * sfid;
		const idx_t Rid = batch_size * efid;
		cout << "Lid = " << Lid
			 << ", Rid = " << Rid
			 << "\n";
		ifstream in;
		open_fstream(in, input, ios::in);
		prw.OpenFiles(sfid, efid, input, generate_partition_file_name);
		
		while (in >> ec) {
			if (ec.qsize < min_read_size || ec.ssize < min_read_size) continue;
			if (ec.qid >= Lid && ec.qid < Rid) {
				normalise_candidate(ec, nec, false);
				prw.WriteOneResult((ec.qid - Lid) / batch_size, ec.qid, nec);
			}
			if (ec.sid >= Lid && ec.sid < Rid)
			{
				normalise_candidate(ec, nec, true);
				prw.WriteOneResult((ec.sid - Lid) / batch_size, ec.sid, nec);
			}
		}
		for (int k = 0; k < nf; ++k)
		{
			if (prw.max_seq_ids[k] == std::numeric_limits<index_t>::min()) continue;
			idx_file << prw.file_names[k] << "\t" << prw.min_seq_ids[k] << "\t" << prw.max_seq_ids[k] << "\n";
			fprintf(stderr, "%s contains reads %d --- %d\n", prw.file_names[k].c_str(), (int)prw.min_seq_ids[k], (int)prw.max_seq_ids[k]);
		}
		prw.CloseFiles();
	}
	close_fstream(idx_file);
}

/*
void
partition_candidates(const char* input, const idx_t batch_size, const int min_read_size, const int num_files)
{
	DynamicTimer dtimer(__func__);
	
	idx_t num_reads = get_num_reads(input);
	const index_t num_batches = (num_reads + batch_size - 1) / batch_size;
	std::string idx_file_name;
    generate_partition_index_file_name(input, idx_file_name);
    std::ofstream idx_file;
    open_fstream(idx_file, idx_file_name.c_str(), std::ios::out);
	
	ExtensionCandidate ec, nec;
	PartitionResultsWriter<ExtensionCandidate> prw;
	for (index_t i = 0; i < num_batches; i += PartitionResultsWriter<ExtensionCandidate>::kNumFiles)
	{
		const index_t sfid = i;
        const index_t efid = std::min(sfid + PartitionResultsWriter<ExtensionCandidate>::kNumFiles, num_batches);
        const int nf = efid - sfid;
        const index_t L = batch_size * sfid;
        const index_t R = batch_size * efid;
        std::ifstream in;
        open_fstream(in, input, std::ios::in);
		prw.OpenFiles(sfid, efid, input, generate_partition_file_name);
		
		while (in >> ec)
		{
			if (ec.qsize < min_read_size || ec.ssize < min_read_size) continue;
			if (ec.qid >= L && ec.qid < R)
			{
				normalise_candidate(ec, nec, false);
				prw.WriteOneResult((ec.qid - L) / batch_size, ec.qid, nec);
			}
			if (ec.sid >= L && ec.sid < R)
			{
				normalise_candidate(ec, nec, true);
				prw.WriteOneResult((ec.sid - L) / batch_size, ec.sid, nec);
			}
		}
		
		for (int k = 0; k < nf; ++k)
        {
            if (prw.max_seq_ids[k] == std::numeric_limits<index_t>::min()) continue;
            idx_file << prw.file_names[k] << "\t" << prw.min_seq_ids[k] << "\t" << prw.max_seq_ids[k] << "\n";
			fprintf(stderr, "%s contains reads %d --- %d\n", prw.file_names[k].c_str(), (int)prw.min_seq_ids[k], (int)prw.max_seq_ids[k]);
        }
	}
    prw.CloseFiles();
	close_fstream(idx_file);
}
*/

/*
void
partition_m4records(const char* m4_file_name, const double min_cov_ratio, const index_t batch_size, const int min_read_size)
{
	DynamicTimer dtimer(__func__);
	
    index_t num_reads, num_qualified_records;
    get_qualified_m4record_counts(m4_file_name, min_cov_ratio, num_qualified_records, num_reads);
	set<index_t> repeat_reads;
	//get_repeat_reads(m4_file_name, min_cov_ratio, num_reads, repeat_reads);
    const index_t num_batches = (num_reads + batch_size - 1) / batch_size;
    std::string idx_file_name;
    generate_partition_index_file_name(m4_file_name, idx_file_name);
    std::ofstream idx_file;
    open_fstream(idx_file, idx_file_name.c_str(), std::ios::out);

    M4Record m4, nm4;
	ExtensionCandidate ec;
    PartitionResultsWriter<ExtensionCandidate> prw;
    for (index_t i = 0; i < num_batches; i += PartitionResultsWriter<ExtensionCandidate>::kNumFiles)
    {
        const index_t sfid = i;
        const index_t efid = std::min(sfid + PartitionResultsWriter<ExtensionCandidate>::kNumFiles, num_batches);
        const int nf = efid - sfid;
        const index_t L = batch_size * sfid;
        const index_t R = batch_size * efid;
        std::ifstream in;
        open_fstream(in, m4_file_name, std::ios::in);
		prw.OpenFiles(sfid, efid, m4_file_name, generate_partition_file_name);

        while (in >> m4)
        {
			if (m4qsize(m4) < min_read_size || m4ssize(m4) < min_read_size) continue;
            if (!check_m4record_mapping_range(m4, min_cov_ratio)) continue;
			if (repeat_reads.find(m4qid(m4)) != repeat_reads.end()
				||
				repeat_reads.find(m4sid(m4)) != repeat_reads.end()) continue;
			
            if (m4qid(m4) >= L && m4qid(m4) < R)
            {
                normalize_m4record(m4, false, nm4);
				m4_to_candidate(nm4, ec);
                prw.WriteOneResult((m4qid(m4) - L) / batch_size, m4qid(m4), ec);
            }
            if (m4sid(m4) >= L && m4sid(m4) < R)
            {
                normalize_m4record(m4, true, nm4);
				m4_to_candidate(nm4, ec);
                prw.WriteOneResult((m4sid(m4) - L) / batch_size, m4sid(m4), ec);
            }
        }

        for (int k = 0; k < nf; ++k)
        {
            if (prw.max_seq_ids[k] == std::numeric_limits<index_t>::min()) continue;
            idx_file << prw.file_names[k] << "\t" << prw.min_seq_ids[k] << "\t" << prw.max_seq_ids[k] << "\n";
			fprintf(stderr, "%s contains reads %d --- %d\n", prw.file_names[k].c_str(), (int)prw.min_seq_ids[k], (int)prw.max_seq_ids[k]);
        }

        prw.CloseFiles();
    }
    close_fstream(idx_file);
}
*/

void
partition_m4records(const char* m4_file_name, 
					const double min_cov_ratio, 
					const index_t batch_size, 
					const int min_read_size,
				    int num_files)
{
	DynamicTimer dtimer(__func__);
	
	num_files = fix_file_counts(num_files);
    index_t num_reads, num_qualified_records;
    get_qualified_m4record_counts(m4_file_name, min_cov_ratio, num_qualified_records, num_reads);
	set<index_t> repeat_reads;
    const index_t num_batches = (num_reads + batch_size - 1) / batch_size;
    std::string idx_file_name;
    generate_partition_index_file_name(m4_file_name, idx_file_name);
    std::ofstream idx_file;
    open_fstream(idx_file, idx_file_name.c_str(), std::ios::out);

    M4Record m4, nm4;
	ExtensionCandidate ec;
    PartitionResultsWriter<ExtensionCandidate> prw(num_files);
    for (index_t i = 0; i < num_batches; i += num_files)
    {
        const index_t sfid = i;
        const index_t efid = std::min(sfid + num_files, num_batches);
        const int nf = efid - sfid;
        const index_t L = batch_size * sfid;
        const index_t R = batch_size * efid;
        std::ifstream in;
        open_fstream(in, m4_file_name, std::ios::in);
		prw.OpenFiles(sfid, efid, m4_file_name, generate_partition_file_name);

        while (in >> m4)
        {
			if (m4qsize(m4) < min_read_size || m4ssize(m4) < min_read_size) continue;
            if (!check_m4record_mapping_range(m4, min_cov_ratio)) continue;
			if (repeat_reads.find(m4qid(m4)) != repeat_reads.end()
					||
					repeat_reads.find(m4sid(m4)) != repeat_reads.end()) continue;
			
            if (m4qid(m4) >= L && m4qid(m4) < R)
            {
                normalize_m4record(m4, false, nm4);
				m4_to_candidate(nm4, ec);
                prw.WriteOneResult((m4qid(m4) - L) / batch_size, m4qid(m4), ec);
            }
            if (m4sid(m4) >= L && m4sid(m4) < R)
            {
                normalize_m4record(m4, true, nm4);
				m4_to_candidate(nm4, ec);
                prw.WriteOneResult((m4sid(m4) - L) / batch_size, m4sid(m4), ec);
            }
        }

        for (int k = 0; k < nf; ++k)
        {
            if (prw.max_seq_ids[k] == std::numeric_limits<index_t>::min()) continue;
            idx_file << prw.file_names[k] << "\t" << prw.min_seq_ids[k] << "\t" << prw.max_seq_ids[k] << "\n";
			fprintf(stderr, "%s contains reads %d --- %d\n", prw.file_names[k].c_str(), (int)prw.min_seq_ids[k], (int)prw.max_seq_ids[k]);
        }

        prw.CloseFiles();
    }
    close_fstream(idx_file);
}

void
load_partition_files_info(const char* idx_file_name, std::vector<PartitionFileInfo>& file_info_vec)
{
    file_info_vec.clear();
    std::ifstream in;
    open_fstream(in, idx_file_name, std::ios::in);
    PartitionFileInfo pfi;
    while (in >> pfi.file_name)
    {
        in >> pfi.min_seq_id >> pfi.max_seq_id;
        file_info_vec.push_back(pfi);
    }
    close_fstream(in);
}
