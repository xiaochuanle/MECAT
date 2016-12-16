#ifndef SPLIT_DATABASE_H
#define SPLIT_DATABASE_H

#include "../common/defs.h"

#define MCS (2140000000L) // max chunk size
//#define MCS 50000000L

typedef struct {
    int offset, size;
} offset_t;

typedef struct {
    int curr, max_size;
    offset_t* offset_list;
} offset_list_t;

typedef struct {
    int num_reads;
    int curr, max_size;
	int start_read_id;
    uint8_t* data;
    offset_list_t* offset_list;
} volume_t;

volume_t*
new_volume_t(int num_reads, int num_bases);

void
clear_volume_t(volume_t* v);

volume_t*
delete_volume_t(volume_t* v);

typedef struct {
	int num_vols, max_num_vols;
	int* name_offsets;
	int buf_size, max_buf_size;
	char* vn_buffer;
} volume_names_t;

volume_names_t*
new_volume_names_t(int num_vols);

volume_names_t*
delete_volume_names_t(volume_names_t* vn);

const char*
get_vol_name(volume_names_t* vn, const int vid);

void
add_one_volume_name(volume_names_t* vn, const char* name, const int ns);

volume_names_t*
load_volume_names(const char* idx_file_name, int num_vols);

void
print_volume_names(volume_names_t* vn);

void 
dump_volume(const char* vol_name, volume_t* v);

int
get_read_id_from_offset_list(offset_list_t* list, const int offset);

volume_t*
load_volume(const char* vol_name);

void
generate_vol_file_name(const char* wrk_dir, int vol, char* vol_file_name);

void
generate_idx_file_name(const char* wrk_dir, char* idx_file_name);

void
split_dataset(const char* reads, const char* wrk_dir, int* num_vols);

void
extract_one_seq(volume_t* v, const int id, char* s);

int
split_raw_dataset(const char* reads, const char* wrk_dir);

#endif // SPLIT_DATABASE_H
