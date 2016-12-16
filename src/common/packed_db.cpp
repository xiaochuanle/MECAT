#include "packed_db.h"

#include <fstream>
#include <string>

#include "defs.h"
#include "fasta_reader.h"

using namespace std;

void
PackedDB::dump_pac(u1_t* p, const idx_t size, const char* path)
{
    ofstream out;
    open_fstream(out, path, std::ios::out | std::ios::binary);
    streambuf* sb = out.rdbuf();
    sb_write(sb, p, (size + 3)/4);
    sb_write(sb, &size, sizeof(idx_t));
    close_fstream(out);
}

u1_t*
PackedDB::load_pac(const char* path, idx_t& size)
{
    ifstream in;
    open_fstream(in, path, std::ios::in | std::ios::binary);
    streambuf* sb = in.rdbuf();
    in.seekg(-sizeof(idx_t), ios::end);
    sb_read(sb, &size, sizeof(idx_t));
    in.seekg(0, ios::beg);
	u1_t* p;
    safe_calloc(p, u1_t, (size+3)/4);
    sb_read(sb, p, (size + 3)/4);
    close_fstream(in);
	return p;
}

void
PackedDB::dump_idx(PODArray<SeqIndex>& idx_list, const char* path)
{
    ofstream out;
    open_fstream(out, path, ios::out);
    idx_t i = 0, n = idx_list.size();
    for (i = 0; i < n; ++i)
        out << idx_list[i].id << "\t" << idx_list[i].offset << "\t" << idx_list[i].size << "\n";
    close_fstream(out);
}

void
PackedDB::load_idx(const char* path, PODArray<SeqIndex>& idx_list)
{
    ifstream in;
    open_fstream(in, path, ios::in);
    SeqIndex si;
    idx_list.clear();
    while(in >> si.id >> si.offset >> si.size) idx_list.push_back(si);
    close_fstream(in);
}

void
PackedDB::dump_packed_db(const char* path)
{
    string n;
    generate_pac_name(path, n);
    dump_pac(pac, db_size, n.data());
    generate_idx_name(path, n);
    dump_idx(seq_idx, n.data());
}

void
PackedDB::load_packed_db(const char* path)
{
    string n;
    generate_pac_name(path, n);
    pac = load_pac(n.data(), db_size);
	max_db_size = db_size;
    generate_idx_name(path, n);
    load_idx(n.data(), seq_idx);
}

void
PackedDB::pack_fasta_db(const char* path, const char* output_prefix, const idx_t min_size)
{
	u1_t* buffer;
	safe_malloc(buffer, u1_t, MAX_SEQ_SIZE);
	Sequence read;
	FastaReader fr(path);
	string n;
	generate_pac_name(output_prefix, n);
	ofstream pout;
	open_fstream(pout, n.data(), ios::out | ios::binary);
	streambuf* psb = pout.rdbuf();
	generate_idx_name(output_prefix, n);
	ofstream iout;
	open_fstream(iout, n.data(), ios::out);
	
	const u1_t* et = get_dna_encode_table();
	idx_t id = 0, tsize = 0, rsize = 0;
	SeqIndex si;
	while(1)
	{
		rsize = fr.read_one_seq(read);
		if (rsize == -1) break;
		if (rsize < min_size) continue;
		Sequence::str_t& s = read.sequence();
		memset(buffer, 0, MAX_SEQ_SIZE);
		for(idx_t i = 0; i < rsize; ++i)
		{
			u1_t c = s[i];
			c = et[c];
			if (c > 3) c = 0;
			set_char(buffer, i, c);
		}
		
		si.offset = tsize;
		si.size = rsize;
		si.id = id++;
		iout << si.id << "\t" << si.offset << "\t" << si.size << "\n";
		
		rsize = (rsize + 3) / 4;
		sb_write(psb, buffer, rsize);
		rsize *= 4;
		tsize += rsize;
	}
	
	sb_write(psb, &tsize, sizeof(idx_t));
	close_fstream(pout);
	close_fstream(iout);
	safe_free(buffer);
	
	LOG(stdout, "pack %lld reads, totally %lld residues", (long long)id, (long long)tsize);
}

void PackedDB::add_one_seq(const Sequence& seq)
{
	SeqIndex si;
	si.size = seq.size();
	si.offset = db_size;
	seq_idx.push_back(si);
	
	if (db_size + si.size > max_db_size)
	{
		idx_t new_size = (max_db_size) ? max_db_size : 1024;
		while (db_size + si.size > new_size) new_size *= 2;
		u1_t* new_pac = NULL;
		safe_calloc(new_pac, u1_t, (new_size + 3)/4);
		memcpy(new_pac, pac, (db_size + 3)/4);
		safe_free(pac);
		pac = new_pac;
		max_db_size = new_size;
	}
	const Sequence::str_t& org_seq = seq.sequence();
	const u1_t* table = get_dna_encode_table();
	for (idx_t i = 0; i < si.size; ++i)
	{
		u1_t c = org_seq[i];
		c = table[c];
		if (c > 3) c = rand() & 3;
		set_char(db_size, c);
		++db_size;
	}
}

void PackedDB::add_one_seq(const char* seq, const idx_t size)
{
	SeqIndex si;
	si.size = size;
	si.offset = db_size;
	seq_idx.push_back(si);
	
	if (db_size + si.size > max_db_size)
	{
		idx_t new_size = (max_db_size) ? max_db_size : 1024;
		while (db_size + si.size > new_size) new_size *= 2;
		u1_t* new_pac = NULL;
		safe_calloc(new_pac, u1_t, (new_size + 3)/4);
		memcpy(new_pac, pac, (db_size + 3)/4);
		safe_free(pac);
		pac = new_pac;
		max_db_size = new_size;
	}
	
	const u1_t* table = get_dna_encode_table();
	for (idx_t i = 0; i < si.size; ++i)
	{
		u1_t c = seq[i];
		c = table[c];
		if (c > 3) c = rand() & 3;
		set_char(db_size, c);
		++db_size;
	}
}

void PackedDB::load_fasta_db(const char* dbname)
{
	DynamicTimer dtimer(__func__);
	FastaReader freader(dbname);
	Sequence seq;
	while (1)
	{
		idx_t size = freader.read_one_seq(seq);
		if (size == -1) break;
		add_one_seq(seq);
	}
}

idx_t PackedDB::offset_to_rid(const idx_t offset) const
{
	if (offset >= db_size) return -1;
	idx_t left, mid, right;
	left = 0, mid = 0, right = seq_idx.size();
	while (left < right)
	{
		mid = (left + right) >> 1;
		if (offset >= seq_idx[mid].offset)
		{
			if (mid == seq_idx.size() - 1) break;
			if (offset < seq_idx[mid + 1].offset) break;
			left = mid + 1;
		}
		else
		{
			right = mid;
		}
	}
	return mid;
}
