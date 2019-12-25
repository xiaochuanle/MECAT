#include "split_database.h"

// make off_t 64 bit (from ftello man page)
#define _FILE_OFFSET_BITS 64

#include <string.h>

#include <fstream>
#include <string>
#include <limits.h>	// PATH_MAX
#include <unistd.h>	// unlink()

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

volume_t* new_volume_t(const int num_reads, int num_bases, const int no_allocate) {
	volume_t* volume = (volume_t*)malloc(sizeof(volume_t));
	volume->num_reads = 0;
	volume->curr = 0;
	if (no_allocate) {
		volume->max_size = 0;
		volume->data = 0;
		volume->offset_list = 0;
	} else {
		if (num_bases == 0) {
			num_bases = MCS + MSS;
		}
		idx_t vol_bytes = (num_bases + 3) / 4;
		volume->max_size = vol_bytes * 4;
		safe_calloc(volume->data, uint8_t, vol_bytes);
		volume->offset_list = new_offset_list_t(num_reads);
	}
	return volume;
}

void
clear_volume_t(volume_t* v)
{
    assert(v);
    v->num_reads = 0;
    v->curr = 0;
    if (v->offset_list) {
        v->offset_list->curr = 0;
    }
	if (v->data) {
		memset(v->data, 0, v->max_size / 4);
	}
}

volume_t*
delete_volume_t(volume_t* v)
{
    delete_offset_list_t(v->offset_list);
    if (v->data) {
        free(v->data);
    }
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
		SET_CHAR(d, idx, c);
		++volume->curr;
	}
}

void
extract_one_seq(const volume_t* v, const int id, char* s)
{
	assert(id < v->num_reads);
	int offset = v->offset_list->offset_list[id].offset;
	int size = v->offset_list->offset_list[id].size;
	int i = 0;
	for (i = 0; i < size; ++i)
	{
		int k = offset + i;
		s[i] = GET_CHAR(v->data, k);
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
load_volume_header(const char* vol_name)
{
	int num_reads, num_bases;
	FILE* in = fopen(vol_name, "rb");
	if (!in) { LOG(stderr, "failed to open file \'%s\'.", vol_name); exit(1); }
	
	// 1) number of reads
	SAFE_READ(&num_reads, int, 1, in);
	// 2) number of bases
	SAFE_READ(&num_bases, int, 1, in);
	
	volume_t* v = new_volume_t(num_reads, num_bases, 1);
	v->num_reads = num_reads;
	v->curr = num_bases;
	// 3) start read id
	SAFE_READ(&v->start_read_id, int, 1, in);
	
	fclose(in);
	return v;
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

void extract_one_seq(ifstream& pac_file, const idx_t offset, const idx_t size, u1_t* const buffer, char* const seq) {
	const idx_t bytes((size + 3) / 4);
	pac_file.seekg(offset / 4, ios::beg);
	pac_file.read((char*)buffer, bytes);
	const char* const dt(get_dna_decode_table());
	idx_t i(0);
	for (; i < size; ++i) {
		seq[i] = dt[GET_CHAR(buffer, i)];
	}
	seq[i] = '\0';
}

class SplitState {
    public:
	int vol, rid;
	idx_t num_reads, num_nucls;
	FILE *idx_file;
	SplitState(const char * const reads, const char * const wrk_dir) : done_(0), rd_pos_(0), fr_(reads) {
		generate_idx_file_name(wrk_dir, idx_file_name_);
		idx_file = fopen(idx_file_name_, "r");
		if (idx_file) {		// number of vols = lines in index file
			vol = 0;
			char *line;
			size_t length(1024);
			safe_malloc(line, char, length);
			ssize_t i;
			while ((i = getline(&line, &length, idx_file)) != -1) {
				if (i > 0 && line[i - 1] == '\n') {
					--i;
				}
				if (i > 0 && line[i - 1] == '\r') {
					--i;
				}
				if (i > 0) {
					++vol;
				}
			}
			fclose(idx_file);
			done_ = 1;
			return;
		}
		idx_file_name_tmp_ = ckpt_file_name_ = idx_file_name_;
		ckpt_file_name_ += ".ckpt";
		idx_file_name_tmp_ += ".tmp";
		std::ifstream ckpt_in(ckpt_file_name_.c_str());
		if (ckpt_in) {
			off_t idx_pos;
			// using off_t instead of std::streampos for rd_pos to allow
			// >> for reading from file
			ckpt_in >> vol >> num_reads >> num_nucls >> idx_pos >> rd_pos_;
			if (!ckpt_in) {
				ERROR("Error reading checkpoint file %s", ckpt_file_name_.c_str());
			}
			rid = num_reads;
			fr_.seekg(rd_pos_);
			idx_file = fopen(idx_file_name_tmp_.c_str(), "r+");
			if (fseeko(idx_file, idx_pos, SEEK_SET) == -1) {
				ERROR("Could not restore checkpoint, fseeko failed");
			}
		} else {
			vol = rid = 0;
			num_reads = num_nucls = 0;
			idx_file = fopen(idx_file_name_tmp_.c_str(), "w");
		}
		if (!idx_file) {
			ERROR("Could not open index file: %s", idx_file_name_tmp_.c_str());
		}
	}
	~SplitState() { }
	idx_t read_one_seq(Sequence& seq) {
		// have to get position before reading, as read might not be
		// covered by next checkpoint
		rd_pos_ = fr_.tellg();
		return fr_.read_one_seq(seq);
	}
	void checkpoint() {
		std::string ckpt_file_name_tmp(ckpt_file_name_ + ".tmp");
		std::ofstream ckpt_out(ckpt_file_name_tmp.c_str());
		if (!ckpt_out) {
			LOG(stderr, "Checkpoint failed at %d, could not open: %s", vol, ckpt_file_name_tmp.c_str());
			return;
		}
		ckpt_out << vol << "\n" << num_reads << "\n" << num_nucls << "\n" << ftello(idx_file) << "\n" << rd_pos_ << "\n";
		if (!ckpt_out) {
			LOG(stderr, "Warning: error writing checkpoint at %d", vol);
			return;
		}
		ckpt_out.close();
		if (rename(ckpt_file_name_tmp.c_str(), ckpt_file_name_.c_str()) == -1) {
			LOG(stderr, "Checkpoint failed at %d, could not rename: %s", vol, ckpt_file_name_tmp.c_str());
		}
	}
	void finish() {
		fclose(idx_file);
		if (rename(idx_file_name_tmp_.c_str(), idx_file_name_) == -1) {
			ERROR("Could not rename index file: %s", idx_file_name_tmp_.c_str());
		}
		unlink(ckpt_file_name_.c_str());
	}
	bool already_done() const {
		return done_;
	}
    private:
	bool done_;
	off_t rd_pos_;
	std::string ckpt_file_name_;	// idx_file_name + ".ckpt"
	std::string idx_file_name_tmp_;	// idx_file_name + ".tmp"
	char idx_file_name_[PATH_MAX];
	FastaReader fr_;
};

int split_raw_dataset(const char* reads, const char* wrk_dir) {
	DynamicTimer dtimer(__func__);
	SplitState state(reads, wrk_dir);
	if (state.already_done()) {
		return state.vol;
	}
	volume_t* v(new_volume_t(0, 0));
	char vol_file_name[PATH_MAX];
	Sequence read;
	for (;;) {
		const idx_t rsize(state.read_one_seq(read));
		if (v->curr + rsize >= MCS || rsize == -1) {
			v->start_read_id = state.rid;
			state.rid += v->num_reads;
			generate_vol_file_name(wrk_dir, state.vol++, vol_file_name);
			if (fprintf(state.idx_file, "%s\n", vol_file_name) < 0) {
				ERROR("Failed to write to index file");
			}
			fflush(state.idx_file);
			dump_volume(vol_file_name, v);
			if (rsize == -1) {
				break;
			}
			state.checkpoint();
			clear_volume_t(v);
		}
		add_one_seq(v, read.sequence().data(), rsize);
		++state.num_reads;
		state.num_nucls += rsize;
		++v->curr;
	}
	state.finish();
	delete_volume_t(v);
	LOG(stderr, "split \'%s\' (%lld reads, %lld nucls) into %d volumes.", reads, (long long)state.num_reads, (long long)state.num_nucls, state.vol);
	return state.vol;
}

void
split_dataset(const char* reads, const char* wrk_dir, int* num_vols)
{
	std::string name(reads);
	name += ".idx";
	ifstream in_idx_file;
	open_fstream(in_idx_file, name.c_str(), ios::in);
	ifstream pac_file;
	name = reads;
	name += ".pac";
	open_fstream(pac_file, name.c_str(), ios::in | ios::binary);
	u1_t* buffer;
	char* seq;
	safe_malloc(buffer, u1_t, MAX_SEQ_SIZE);
	safe_malloc(seq, char, MAX_SEQ_SIZE);
	volume_t* v = new_volume_t(0, 0);
	int vol = 0;
	int rid = 0;
	char idx_file_name[PATH_MAX], vol_file_name[PATH_MAX];
	generate_idx_file_name(wrk_dir, idx_file_name);
	FILE* idx_file = fopen(idx_file_name, "w");
	idx_t offset, size;
	while (in_idx_file >> offset >> size)
	{
		if (v->curr + size + 1 > MCS)
		{
			v->start_read_id = rid;
			rid += v->num_reads;
			generate_vol_file_name(wrk_dir, vol++, vol_file_name);
			fprintf(idx_file, "%s\n", vol_file_name);
			dump_volume(vol_file_name, v);
			clear_volume_t(v);
		}
		extract_one_seq(pac_file, offset, size, buffer, seq);
		add_one_seq(v, seq, size);
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
	// allow testing to see if file exists
	if (!fvn && num_vols == 0) {
		return vn;
	}
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
		std::cerr << name << "\n";
	}
}
