#ifndef OVERLAPS_PARTITION_H
#define OVERLAPS_PARTITION_H

#include <string>	// string
#include <vector>	// vector<>
#include <stdint.h>	// int64_t

void generate_partition_index_file_name(const std::string& input_file_name, std::string& ret);

void partition_m4records(const char* m4_file_name, double min_cov_ratio, size_t batch_size, int min_read_size, int num_files);

void partition_candidates(const std::string& input, const std::string& pac_prefix, size_t batch_size, int num_files, int64_t num_reads);

void load_partition_files_info(const char* idx_file_name, std::vector<std::string>& file_info_vec);

#endif // OVERLAPS_PARTITION_H
