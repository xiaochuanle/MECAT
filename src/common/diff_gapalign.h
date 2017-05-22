#ifndef DIFF_GAPALIGN_H
#define DIFF_GAPALIGN_H

#include "defs.h"
#include "alignment.h"
#include "gapalign.h"

#include <string>
#include <vector>

typedef long idx;

struct DiffAlignParameters
{
    idx segment_size;
    idx row_size;
    idx column_size;
    idx segment_aln_size;
    idx max_seq_size;
    idx max_aln_size;
    idx d_path_size;
    idx aln_path_size;
	
	void init(const int large_block = 0) {
		if (large_block) {
            segment_size = 1000;
            row_size = 4096;
            column_size = 4096;
            segment_aln_size = 4096;
            max_seq_size = MAX_SEQ_SIZE;
            max_aln_size = MAX_SEQ_SIZE;
            d_path_size = 5000000;
            aln_path_size = 5000000;
        } else {
            segment_size = 500;
            row_size = 4096;
            column_size = 4096;
            segment_aln_size = 4096;
            max_seq_size = MAX_SEQ_SIZE;
            max_aln_size = MAX_SEQ_SIZE;
            d_path_size = 5000000;
            aln_path_size = 5000000;
        }
	}
};

struct Alignment
{
    int aln_str_size;
    int dist;
    int aln_q_s;
    int aln_q_e;
    int aln_t_s;
    int aln_t_e;
    char* q_aln_str;
    char* t_aln_str;

    void init() {
        aln_str_size = 0;
        aln_q_s = aln_q_e = 0;
        aln_t_s = aln_t_e = 0;
    }
	
	Alignment(const idx max_aln_size) {
		snew(q_aln_str, char, max_aln_size);
		snew(t_aln_str, char, max_aln_size);
	}
	
	~Alignment() {
		sfree(q_aln_str);
		sfree(t_aln_str);
	}
};

struct OutputStore
{
	char* left_store1;
	char* left_store2;
	char* right_store1;
	char* right_store2;
	char* out_store1;
	char* out_store2;
	
	int left_store_size;
	int right_store_size;
	int out_store_size;
	int query_start, query_end;
	int target_start, target_end;
	
	double calc_ident() const {
		if (out_store_size == 0) return 0.0;
		int n = 0;
		for (int i = 0; i < out_store_size; ++i) {
			if (out_store1[i] == out_store2[i]) ++n;
		}
		return 100.0 * n / out_store_size;
	}
	
	OutputStore(const idx max_aln_size) {
		snew(left_store1, char, max_aln_size);
		snew(left_store2, char, max_aln_size);
		snew(right_store1, char, max_aln_size);
		snew(right_store2, char, max_aln_size);
		snew(out_store1, char, max_aln_size);
		snew(out_store2, char, max_aln_size);
	}
	
	~OutputStore() {
		sfree(left_store1);
		sfree(left_store2);
		sfree(right_store1);
		sfree(right_store2);
		sfree(out_store1);
		sfree(out_store2);
	}
	
	void init() {
		left_store_size = right_store_size = out_store_size = 0;
	}
};

struct DPathData
{
    int pre_k, x1, y1, x2, y2;
};

struct DPathData2
{
    int d, k, pre_k, x1, y1, x2, y2;
};

struct PathPoint
{
    int x, y;
};

class DiffAligner : public GapAligner
{
public:
    DiffAligner(const int large_block) {
        param.init(large_block);
		align = new Alignment(param.segment_aln_size);
		result = new OutputStore(param.max_aln_size);
        snew(dynq, int, param.row_size);
        snew(dynt, int, param.column_size);
        snew(d_path, DPathData2, param.d_path_size);
        snew(aln_path, PathPoint, param.aln_path_size);
    }
    
    virtual ~DiffAligner() {
		delete align;
		delete result;
        sfree(dynq);
        sfree(dynt);
        sfree(d_path);
        sfree(aln_path);
    }

	virtual bool go(const char* query, const int qstart, const int qsize, 
			const char* target, const int sstart, const int tsize,
			const int min_aln_size);
	
	virtual int query_start() const {
		return result->query_start;
	}
	
	virtual int query_end() const {
		return result->query_end;
	}
	
	virtual int target_start() const {
		return result->target_start;
	}
	
	virtual int target_end() const {
		return result->target_end;
	}
	
	virtual double calc_ident() const {
		return result->calc_ident();
	}
	
	virtual char* query_mapped_string() {
		return result->out_store1;
	}
	
	virtual char* target_mapped_string() {
		return result->out_store2;
	}

public:
    DiffAlignParameters    	param;
    int*           		 	dynq;
    int*            		dynt;
    Alignment*       		align;
    OutputStore*     		result;
    DPathData2*     		d_path;
    PathPoint*      		aln_path;
};

#endif // DIFF_GAPALIGN_H
