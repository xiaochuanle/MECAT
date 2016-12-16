#include "split_database.h"

#include <string.h>

#include <fstream>
#include <string>

#include "packed_db.h"
#include "fasta_reader.h"

#define MSS MAX_SEQ_SIZE

using namespace std;

int
get_read_id_from_offset_list(offset_list_t* list, const int offset)
{
	int n = list->curr;
	offset_t* a = list->offset_list;
	int left = 0, right = n - 1, mid = (left + right) / 2;
	if (a[right].offset < offset) return right;
	while (left <= right)
	{
		if (a[mid].offset <= offset && a[mid].offset + a[mid].size > offset) return mid;
		if (a[mid].offset + a[mid].size <= offset) left = mid + 1;
		else if (a[mid].offset > offset) right = mid - 1;
		else
		{
			LOG(stderr, "Error!");
			exit(1);
		}
		mid = (left + right) / 2;
	}
	return mid;
}

offset_list_t*
new_offset_list_t(int size)
{
    offset_list_t* list = (offset_list_t*)malloc(sizeof(offset_list_t));
    list->curr = 0;
	if (size == 0) size = 100000;
    list->max_size = size;
    safe_malloc(list->offset_list, offset_t, size);
    return list;
}

offset_list_t*
delete_offset_list_t(offset_list_t* list)
{
    if (!list) return list;
    if (list->offset_list) free(list->offset_list);
    free(list);
    return NULL;
}

void
insert_one_offset(offset_list_t* list, const int offset, const int size)
{
    if (list->curr >= list->max_size)
    {
        list->max_size *= 2;
        safe_realloc(list->offset_list, offset_t, list->max_size);
    }
    list->offset_list[list->curr].offset = offset;
	list->offset_list[list->curr].size = size;
	++list->curr;
}

volume_t*
new_volume_t(int num_reads, int num_bases)
{
    volume_t* volume = (volume_t*)malloc(sizeof(volume_t));
    volume->num_reads = 0;
    volume->curr = 0;
	if (num_bases == 0) num_bases = (MCS + MSS);
    volume->max_size = num_bases;
    idx_t vol_bytes = (num_bases + 3) / 4;
    safe_calloc(volume->data, uint8_t, vol_bytes);
    volume->offset_list = new_offset_list_t(num_reads);
    return volume;
}

void
clear_volume_t(volume_t* v)
{
    assert(v);
    v->num_reads = 0;
    v->curr = 0;
    v->offset_list->curr = 0;
	memset(v->data, 0, v->max_size / 4);
}

volume_t*
delete_volume_t(volume_t* v)
{
    v->offset_list = delete_offset_list_t(v->offset_list);
    free(v->data);
    free(v);
    return NULL;
}

void
add_one_seq(volume_t* volume, const char* s, const int size)
{
	++volume->num_reads;
	insert_one_offset(volume->offset_list, volume->curr, size);
	const uint8_t* encode_table = get_dna_encode_table();
	int i;
	for (i = 0; i < size; ++i) 
	{
		uint8_t* d = volume->data;
		int idx = volume->curr;
		uint8_t c = s[i];
		c = encode_table[c];
		PackedDB::set_char(d, idx, c);
		++volume->curr;
	}
}

void
extract_one_seq(volume_t* v, const int id, char* s)
{
	assert(id < v->num_reads);
	int offset = v->offset_list->offset_list[id].offset;
	int size = v->offset_list->offset_list[id].size;
	int i = 0;
	for (i = 0; i < size; ++i)
	{
		int k = offset + i;
		s[i] = PackedDB::get_char(v->data, k);
	}
}

void 
dump_volume(const char* vol_name, volume_t* v)
{
	FILE* out = fopen(vol_name, "wb");
	assert(out);
	// 1) number of reads
	SAFE_WRITE(&v->num_reads, int, 1, out);
	// 2) number of bases
	SAFE_WRITE(&v->curr, int, 1, out);
	// 3) start read id
	SAFE_WRITE(&v->start_read_id, int, 1, out);
	// 4) offset list
	assert(v->offset_list->curr == v->num_reads);
	SAFE_WRITE(v->offset_list->offset_list, offset_t, v->num_reads, out);
	// 5) pac
	int vol_bytes = (v->curr + 3) / 4;
	SAFE_WRITE(v->data, uint8_t, vol_bytes, out);
	fclose(out);
}

volume_t*
load_volume(const char* vol_name)
{
	int num_reads, num_bases;
	FILE* in = fopen(vol_name, "rb");
	if (!in) { LOG(stderr, "failed to open file \'%s\'.", vol_name); exit(1); }
	
	// 1) number of reads
	SAFE_READ(&num_reads, int, 1, in);
	// 2) number of bases
	SAFE_READ(&num_bases, int, 1, in);
	
	volume_t* v = new_volume_t(num_reads, num_bases);
	v->num_reads = num_reads;
	v->curr = num_bases;
	// 3) start read id
	SAFE_READ(&v->start_read_id, int, 1, in);
	// 4) offset list
	SAFE_READ(v->offset_list->offset_list, offset_t, num_reads, in);
	v->offset_list->curr = num_reads;
	// 5) pac
	int vol_bytes = (num_bases + 3) / 4;
	SAFE_READ(v->data, uint8_t, vol_bytes, in);
	
	fclose(in);
	return v;
}

void
generate_vol_file_name(const char* wrk_dir, int vol, char* vol_file_name)
{
	char buffer[64];
	strcpy(vol_file_name, wrk_dir);
	if (vol_file_name[strlen(vol_file_name) - 1] != '/') strcat(vol_file_name, "/");
	strcat(vol_file_name, "vol");
	sprintf(buffer, "%d", vol);
	strcat(vol_file_name, buffer);
}

void
generate_idx_file_name(const char* wrk_dir, char* idx_file_name)
{
	strcpy(idx_file_name, wrk_dir);
	if (idx_file_name[strlen(idx_file_name) - 1] != '/') strcat(idx_file_name, "/");
	strcat(idx_file_name, "fileindex.txt");
}

void
extract_one_seq(ifstream& pac_file, PackedDB::SeqIndex& si, u1_t* buffer, char* seq)
{
	idx_t offset = si.offset / 4;
	idx_t bytes = (si.size + 3) / 4;
	pac_file.seekg(offset, ios::beg);
	memset(buffer, 0, MAX_SEQ_SIZE);
	pac_file.read((char*)buffer, bytes);
	const char* dt = get_dna_decode_table();
	idx_t i = 0;
	for(i = 0; i < si.size; ++i)
	{
		u1_t c = PackedDB::get_char(buffer, i);
		c = dt[c];
		seq[i] = c;
	}
	seq[i] = '\0';
}

int
split_raw_dataset(const char* reads, const char* wrk_dir)
{
	DynamicTimer dtimer(__func__);
	volume_t* v = new_volume_t(0, 0);
	int vol = 0;
	int rid = 0;
	char idx_file_name[1024], vol_file_name[1024];
	generate_idx_file_name(wrk_dir, idx_file_name);
	FILE* idx_file = fopen(idx_file_name, "w");
	FastaReader fr(reads);
	Sequence read;
	idx_t num_reads = 0, num_nucls = 0;
	while (1)
	{
		idx_t rsize = fr.read_one_seq(read);
		if (rsize == -1) break;
		++num_reads;
		num_nucls += rsize;
		if (v->curr + rsize + 1 > MCS)
		{
			v->start_read_id = rid;
			rid += v->num_reads;
			generate_vol_file_name(wrk_dir, vol++, vol_file_name);
			fprintf(idx_file, "%s\n", vol_file_name);
			dump_volume(vol_file_name, v);
			clear_volume_t(v);
		}
		add_one_seq(v, read.sequence().data(), rsize);
		++v->curr;
	}
	
	if (v->curr > 0)
	{
		v->start_read_id = rid;
		rid += v->num_reads;
		generate_vol_file_name(wrk_dir, vol++, vol_file_name);
		fprintf(idx_file, "%s\n", vol_file_name);
		dump_volume(vol_file_name, v);
		clear_volume_t(v);
	}
	fclose(idx_file);
	delete_volume_t(v);
	LOG(stderr, "split \'%s\' (%lld reads, %lld nucls) into %d volumes.", reads, (long long)num_reads, (long long)num_nucls, vol);
	return vol;
}

void
split_dataset(const char* reads, const char* wrk_dir, int* num_vols)
{
	string name;
	PackedDB::generate_idx_name(reads, name);
	ifstream in_idx_file;
	open_fstream(in_idx_file, name.c_str(), ios::in);
	ifstream pac_file;
	PackedDB::generate_pac_name(reads, name);
	open_fstream(pac_file, name.c_str(), ios::in | ios::binary);
	u1_t* buffer;
	char* seq;
	safe_malloc(buffer, u1_t, MAX_SEQ_SIZE);
	safe_malloc(seq, char, MAX_SEQ_SIZE);
	volume_t* v = new_volume_t(0, 0);
	int vol = 0;
	int rid = 0;
	char idx_file_name[1024], vol_file_name[1024];
	generate_idx_file_name(wrk_dir, idx_file_name);
	FILE* idx_file = fopen(idx_file_name, "w");
	PackedDB::SeqIndex si;
	while (in_idx_file >> si.id >> si.offset >> si.size)
	{
		if (v->curr + si.size + 1 > MCS)
		{
			v->start_read_id = rid;
			rid += v->num_reads;
			generate_vol_file_name(wrk_dir, vol++, vol_file_name);
			fprintf(idx_file, "%s\n", vol_file_name);
			dump_volume(vol_file_name, v);
			clear_volume_t(v);
		}
		extract_one_seq(pac_file, si, buffer, seq);
		add_one_seq(v, seq, si.size);
		++v->curr;
	}
	
	if (v->curr > 0)
	{
		v->start_read_id = rid;
		rid += v->num_reads;
		generate_vol_file_name(wrk_dir, vol++, vol_file_name);
		fprintf(idx_file, "%s\n", vol_file_name);
		dump_volume(vol_file_name, v);
		clear_volume_t(v);
	}
	
	fclose(idx_file);
	safe_free(buffer);
	safe_free(seq);
	delete_volume_t(v);
	*num_vols = vol;
	close_fstream(pac_file);
	close_fstream(in_idx_file);
}

volume_names_t*
new_volume_names_t(int num_vols)
{
	if (num_vols == 0) num_vols = 1000;
	volume_names_t* vn = (volume_names_t*)malloc(sizeof(volume_names_t));
	vn->num_vols = 0;
	vn->max_num_vols = num_vols;
	safe_malloc(vn->name_offsets, int, num_vols);
	vn->buf_size = 0;
	vn->max_buf_size = 100000;
	safe_malloc(vn->vn_buffer, char, vn->max_buf_size);
	return vn;
}

volume_names_t*
delete_volume_names_t(volume_names_t* vn)
{
	free(vn->name_offsets);
	free(vn->vn_buffer);
	free(vn);
	return NULL;
}

const char*
get_vol_name(volume_names_t* vn, const int vid)
{
	return vn->vn_buffer + vn->name_offsets[vid];
}

void
add_one_volume_name(volume_names_t* vn, const char* name, const int ns)
{
	if (vn->num_vols + 1 >= vn->max_num_vols)
	{
		vn->max_num_vols *= 2;
		safe_realloc(vn->name_offsets, int, vn->max_num_vols);
	}
	vn->name_offsets[vn->num_vols++] = vn->buf_size;
	
	if (vn->buf_size + ns + 1 >= vn->max_buf_size)
	{
		vn->max_buf_size *= 2;
		safe_realloc(vn->vn_buffer, char, vn->max_buf_size);
	}
	memcpy(vn->vn_buffer + vn->buf_size, name, ns);
	vn->buf_size += ns;
	vn->vn_buffer[vn->buf_size++] = '\0';
}

volume_names_t*
load_volume_names(const char* idx_file_name, int num_vols)
{
	volume_names_t* vn = new_volume_names_t(num_vols);
	FILE* fvn = fopen(idx_file_name, "r");
	assert(fvn);
	char* name;
	size_t ns = 1024;
	safe_malloc(name, char, ns);
	ssize_t ls;
	while(-1 != (ls = getline(&name, &ns, fvn))) 
	{
		if (ls > 0 && name[ls - 1] == '\n') --ls;
		if (ls > 0 && name[ls - 1] == '\r') --ls;
		if (ls > 0) add_one_volume_name(vn, name, ls);
	}
	fclose(fvn);
	free(name);
	return vn;
}

void
print_volume_names(volume_names_t* vn)
{
	int num_vols = vn->num_vols;
	int i;
	for (i = 0; i < num_vols; ++i)
	{
		const char* name = get_vol_name(vn, i);
		fprintf(stderr, "%s\n", name);
	}
}
