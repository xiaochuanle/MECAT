#include "reads_correction_m4.h"

#include <algorithm>	// sort()
#include <cassert>	// assert()
#include <sstream>
#include <string>	// string

#include "../common/defs.h"	// TECH_PACBIO
#include "mecat_correction.h"
#include "overlaps_partition.h"
#include "overlaps_store.h"
#include "reads_correction_aux.h"

void* reads_correction_func_m4(void* arg) {
	ConsensusThreadData& data(*(static_cast<ConsensusThreadData*>(arg)));
	const int tid(data.get_thread_id());
	ConsensusPerThreadData& pdata(data.data[tid]);
	ExtensionCandidate* const overlaps((ExtensionCandidate*)pdata.candidates);
	const int64_t num_ovlps(pdata.num_candidates);
	int64_t i = 0, j;
	while (i < num_ovlps) {
		const int64_t sid = overlaps[i].sid;
		j = i + 1;
		while (j < num_ovlps && overlaps[j].sid == sid) {
			++j;
		}
		if (j - i < data.rco.min_cov) {
			i = j;
			continue;
		}
		if (overlaps[i].ssize < data.rco.min_size * 0.95) {
			i = j;
			continue;
		}
		if (data.rco.tech == TECH_PACBIO) {
			consensus_one_read_m4_pacbio(data, pdata, sid, i, j);
		} else {
			consensus_one_read_m4_nanopore(data, pdata, sid, i, j);
		}
		if (pdata.cns_results.size() >= MAX_CNS_RESULTS) {
			{
				std::lock_guard<std::mutex> lock(data.out_lock);
				for (std::vector<CnsResult>::iterator iter = pdata.cns_results.begin(); iter != pdata.cns_results.end(); ++iter) {
					data.out << ">" << iter->id << "_" << iter->range[0] << "_" << iter->range[1] << "_" << iter->seq.size() << "\n" << iter->seq << "\n";
					if (!data.out) {
						ERROR("Error writing output");
					}
				}
			}
			pdata.cns_results.clear();
		}
		i = j;
	}
	return NULL;
}

void consensus_one_partition_m4(const char* m4_file_name, ReadsCorrectionOptions& rco, PackedDB& reads, std::ostream& out) {
	int64_t nec;
	ExtensionCandidate* ec_list = load_partition_data<ExtensionCandidate>(m4_file_name, nec);
	std::sort(ec_list, ec_list + nec, CmpExtensionCandidateBySidAndScore());
	ConsensusThreadData data(rco, reads, out, m4_file_name);
	allocate_ecs(data, ec_list, nec);
	pthread_t thread_ids[rco.num_threads];
	for (int i = 0; i < rco.num_threads; ++i) {
		pthread_create(&thread_ids[i], NULL, reads_correction_func_m4, static_cast<void*>(&data));
	}
	for (int i = 0; i < rco.num_threads; ++i) {
		pthread_join(thread_ids[i], NULL);
	}
	for (int i = 0; i < rco.num_threads; ++i) {
		std::vector<CnsResult>& cns_results = data.data[i].cns_results;
		for (std::vector<CnsResult>::iterator iter = cns_results.begin(); iter != cns_results.end(); ++iter) {
			out << ">" << iter->id << "_" << iter->range[0] << "_" << iter->range[1] << "_" << iter->seq.size() << "\n";
			std::string& seq = iter->seq;
			out << seq << "\n";
		}
	}
	delete[] ec_list;
}

static int reads_correction_m4_p(ReadsCorrectionOptions& rco, std::vector<std::string> &partition_file_vec, PackedDB &reads) {
	std::ostringstream os;
	os << rco.corrected_reads << "." << rco.job_index;
	std::string results_file = os.str();
	std::string working_file = results_file + ".working";
	std::ofstream out;
	open_fstream(out, working_file.c_str(), std::ios::out);
	const std::string &p = partition_file_vec[rco.job_index];
	char process_info[1024];
	sprintf(process_info, "processing %s", p.c_str());
	DynamicTimer dtimer(process_info);
	consensus_one_partition_m4(p.c_str(), rco, reads, out);
	assert(rename(working_file.c_str(), results_file.c_str()) == 0);
	return 0;
}

int reads_correction_m4(ReadsCorrectionOptions& rco) {
	std::string idx_file_name;
	generate_partition_index_file_name(rco.m4, idx_file_name);
	std::vector<std::string> partition_file_vec;
	load_partition_files_info(idx_file_name.c_str(), partition_file_vec);
	PackedDB reads;
	reads.load_fasta_db(rco.reads);
	if (rco.job_index != -1) {
		return reads_correction_m4_p(rco, partition_file_vec, reads);
	} else {
		std::ofstream out;
		open_fstream(out, rco.corrected_reads, std::ios::out);
		char process_info[1024];
		for (std::vector<std::string>::iterator iter = partition_file_vec.begin(); iter != partition_file_vec.end(); ++iter)
		{
			sprintf(process_info, "processing %s", iter->c_str());
			DynamicTimer dtimer(process_info);
			consensus_one_partition_m4(iter->c_str(), rco, reads, out);
		}
		
		return 0;
	}
}
