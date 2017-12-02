#include "reads_correction_m4.h"

#include <cstring>

#include "mecat_correction.h"
#include "overlaps_partition.h"
#include "overlaps_store.h"

void*
reads_correction_func_m4(void* arg)
{
    ConsensusThreadData& cns_data = *static_cast<ConsensusThreadData*>(arg);
	ExtensionCandidate* overlaps = cns_data.candidates;
	const index_t num_ovlps = cns_data.num_candidates;
    index_t i = 0, j;
    while (i < num_ovlps)
    {
        const index_t sid = overlaps[i].sid;
        j = i + 1;
        while (j < num_ovlps && overlaps[j].sid == sid) ++j;
        if (j - i < cns_data.rco.min_cov) { i = j; continue; }
        if (overlaps[i].ssize < cns_data.rco.min_size * 0.95) { i = j; continue; }
		if (cns_data.rco.tech == TECH_PACBIO) {
			ns_meap_cns::consensus_one_read_m4_pacbio(&cns_data, sid, i, j);
		} else {
			ns_meap_cns::consensus_one_read_m4_nanopore(&cns_data, sid, i, j);
		}
		if (cns_data.cns_results.size() >= MAX_CNS_RESULTS)
		{
			pthread_mutex_lock(&cns_data.out_lock);
			for (std::vector<CnsResult>::iterator iter = cns_data.cns_results.begin(); iter != cns_data.cns_results.end(); ++iter)
			{
				(*cns_data.out) << ">" << iter->id << "_" << iter->range[0] << "_" << iter->range[1] << "_" << iter->seq.size() << "\n";
				std::string& seq = iter->seq;
				(*cns_data.out) << seq << "\n";
			}
			cns_data.cns_results.clear();
			pthread_mutex_unlock(&cns_data.out_lock);
		}
        i = j;
    }
    return NULL;
}

void
consensus_one_partition_m4(const char* m4_file_name,
        const index_t min_read_id,
        const index_t max_read_id,
        ReadsCorrectionOptions& rco,
        PackedDB& reads,
        std::ostream& out)
{
	idx_t num_ec;
	ExtensionCandidate* ec_list = load_partition_data<ExtensionCandidate>(m4_file_name, num_ec);
    ConsensusThreadData* pctds[rco.num_threads];
	build_cns_thrd_data_can(ec_list, num_ec, min_read_id, max_read_id, &rco, &reads, &out, pctds);
    pthread_t thread_ids[rco.num_threads];
    for (int i = 0; i < rco.num_threads; ++i)
        pthread_create(&thread_ids[i], NULL, reads_correction_func_m4, static_cast<void*>(pctds[i]));
    for (int i = 0; i < rco.num_threads; ++i)
        pthread_join(thread_ids[i], NULL);
	for (int i = 0; i < rco.num_threads; ++i)
	{
		std::vector<CnsResult>& cns_results = pctds[i]->cns_results;
		for (std::vector<CnsResult>::iterator iter = cns_results.begin(); iter != cns_results.end(); ++iter)
		{
			out << ">" << iter->id << "_" << iter->range[0] << "_" << iter->range[1] << "_" << iter->seq.size() << "\n";
			std::string& seq = iter->seq;
			out << seq << "\n";
		}
	}

    delete[] ec_list;
    for (int i = 0; i < rco.num_threads; ++i) delete pctds[i];
}

int reads_correction_m4(ReadsCorrectionOptions& rco)
{
    double mapping_ratio = rco.min_mapping_ratio - 0.02;
	partition_m4records(rco.m4, mapping_ratio, rco.batch_size, rco.min_size, rco.num_partition_files);
	std::string idx_file_name;
	generate_partition_index_file_name(rco.m4, idx_file_name);
	std::vector<PartitionFileInfo> partition_file_vec;
	load_partition_files_info(idx_file_name.c_str(), partition_file_vec);
	PackedDB reads;
	reads.load_fasta_db(rco.reads);
	std::ofstream out;
	open_fstream(out, rco.corrected_reads, std::ios::out);
	char process_info[1024];
	for (std::vector<PartitionFileInfo>::iterator iter = partition_file_vec.begin(); iter != partition_file_vec.end(); ++iter)
	{
		sprintf(process_info, "processing %s", iter->file_name.c_str());
		DynamicTimer dtimer(process_info);
		consensus_one_partition_m4(iter->file_name.c_str(), iter->min_seq_id, iter->max_seq_id, rco, reads, out);
	}
	
	return 0;
}
