#ifndef PACKED_DB_H
#define PACKED_DB_H

#include "defs.h"
#include "sequence.h"

class PackedDB
{
public:
    struct SeqIndex
    {
        idx_t id;
        idx_t offset;
        idx_t size;
    };

public:
    PackedDB() : pac(NULL), db_size(0), max_db_size(0) {}
    ~PackedDB() { destroy(); }
	void reserve(const idx_t& size)
	{
		destroy();
		seq_idx.clear();
		max_db_size = size;
		idx_t bytes = (max_db_size / 4);
		safe_calloc(pac, u1_t, bytes);
	}
	
	void GetSequence(const index_t id, const bool fwd, char* seq, const index_t size_in_ovlp)
	{
		const index_t offset = seq_idx[id].offset;
		const index_t size = seq_idx[id].size;
		r_assert(size == size_in_ovlp);
		index_t idx = 0;
		if (fwd)
		{
			for (index_t i = 0; i < size; ++i)
			{
				uint1 c = get_char(offset + i);
				seq[idx++] = c;
			}
		}
		else
		{
			for (index_t i = size - 1; i >= 0; --i)
			{
				uint1 c = get_char(offset + i);
				c = 3 - c;
				seq[idx++] = c;
			}
		}
	}

    void get_sequence(const idx_t from, const idx_t to, const bool forward, char* seq) const
    {
        idx_t idx = 0;
        if (forward)
            for(idx_t i = from; i < to; ++i)
            {
                u1_t c = get_char(pac, i);
                seq[idx++] = c;
            }
        else
            for(idx_t i = to - 1; i >= from; --i)
            {
                u1_t c = get_char(pac, i);
                c = 3 - c;
                seq[idx++] = c;
            }
    }

    void get_sequence(const idx_t rid, const bool forward, char* seq) const
    {
        const idx_t s = seq_idx[rid].offset;
        const idx_t e = s + seq_idx[rid].size;
        get_sequence(s, e, forward, seq);
    }

    void get_sequence(const idx_t rid, const idx_t from, const idx_t to, const bool forward, char* seq) const
    {
        const idx_t s = seq_idx[rid].offset + from;
        const idx_t e = seq_idx[rid].offset + to;
        get_sequence(s, e, forward, seq);
    }
	
	void get_decode_sequence(const idx_t rid, const idx_t from, const idx_t to, const bool forward, char* seq) const
	{
		get_sequence(rid, from, to, forward, seq);
		for(idx_t i = 0; i < to - from; ++i)
		{
			u1_t c = seq[i];
			r_assert(c >= 0 && c < 4);
			c = "ACGT"[c];
			seq[i] = c;
		}
	}

    static void set_char(u1_t* p, const idx_t idx, const u1_t c)
    {
        p[idx >> 2] |= c << ((~idx&3)<<1);
    }

    static u1_t get_char(const u1_t* p, const idx_t idx)
    {
        u1_t c = p[idx >> 2] >> ((~idx&3)<<1)&3;
        return c;
    }
	
	void set_char(const idx_t idx, const u1_t c)
	{
		set_char(pac, idx, c);
	}
	u1_t get_char(const idx_t idx) const
	{
		return get_char(pac, idx);
	}

    idx_t size() const { return db_size; }
    idx_t num_seqs() const { return seq_idx.size(); }
    idx_t seq_offset(const idx_t rid) const { return seq_idx[rid].offset; }
    idx_t seq_size(const idx_t rid) const { return seq_idx[rid].size; }
    idx_t offset_to_rid(const idx_t offset) const;
    void add_one_seq(const Sequence& seq);
	void add_one_seq(const char* seq, const idx_t size);
    void destroy() { if (pac) safe_free(pac); db_size = max_db_size = 0; }
	void clear() { seq_idx.clear(); db_size = 0; memset(pac, 0, (max_db_size + 3)/4); }

    static void generate_pac_name(const char* prefix, std::string& ret)
    {
        ret = prefix;
        ret += ".pac";
    }
    static void generate_idx_name(const char* prefix, std::string& ret)
    {
        ret = prefix;
        ret += ".idx";
    }
    static void dump_pac(u1_t* p, const idx_t size, const char* path);
    static u1_t* load_pac(const char* path, idx_t& size);
    static void dump_idx(PODArray<SeqIndex>& idx_list, const char* path);
    static void load_idx(const char* path, PODArray<SeqIndex>& idx_list);
	
    void dump_packed_db(const char* path);
    void load_packed_db(const char* path);
	
	static void pack_fasta_db(const char* fasta, const char* output_prefix, const idx_t min_size);
	void load_fasta_db(const char* fasta);

private:
    u1_t*   pac;
    idx_t   db_size;
    idx_t   max_db_size;
    PODArray<SeqIndex> seq_idx;
};

#endif // PACKED_DB_H
