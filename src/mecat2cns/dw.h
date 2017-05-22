#ifndef DW_H
#define DW_H

#include <algorithm>

#include "../common/alignment.h"
#include "../common/defs.h"
#include "../common/packed_db.h"

namespace ns_banded_sw {

struct SW_Parameters
{
    idx_t segment_size;
    idx_t row_size;
    idx_t column_size;
    idx_t segment_aln_size;
    idx_t max_seq_size;
    idx_t max_aln_size;
    idx_t d_path_size;
    idx_t aln_path_size;
};

SW_Parameters
get_sw_parameters_small();

SW_Parameters
get_sw_parameters_large();

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

    void init()
    {
        aln_str_size = 0;
        aln_q_s = aln_q_e = 0;
        aln_t_s = aln_t_e = 0;
    }

    Alignment(const idx_t max_aln_size)
    {
        safe_malloc(q_aln_str, char, max_aln_size);
        safe_malloc(t_aln_str, char, max_aln_size);
    }
    ~Alignment()
    {
        safe_free(q_aln_str);
        safe_free(t_aln_str);
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
    char* out_match_pattern;

    int left_store_size;
    int right_store_size;
    int out_store_size;
    int query_start, query_end;
    int target_start, target_end;
    int mat, mis, ins, del;
	double ident;
    
    OutputStore(const idx_t max_aln_size)
    {
        safe_malloc(left_store1, char, max_aln_size);
        safe_malloc(left_store2, char, max_aln_size);
        safe_malloc(right_store1, char, max_aln_size);
        safe_malloc(right_store2, char, max_aln_size);
        safe_malloc(out_store1, char, max_aln_size);
        safe_malloc(out_store2, char, max_aln_size);
        safe_malloc(out_match_pattern, char, max_aln_size);
    }

    ~OutputStore()
    {
        safe_free(left_store1);
        safe_free(left_store2);
        safe_free(right_store1);
        safe_free(right_store2);
        safe_free(out_store1);
        safe_free(out_store2);
        safe_free(out_match_pattern);
    }

    void init()
    {
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

struct DiffRunningData
{
    SW_Parameters   swp;
    char*           query;
    char*           target;
    int*            DynQ;
    int*            DynT;
    Alignment*      align;
    OutputStore*    result;
    DPathData2*     d_path;
    PathPoint*      aln_path;
	
	DiffRunningData(const SW_Parameters& swp_in);
	~DiffRunningData();
};

void fill_m4record_from_output_store(const OutputStore& result, 
									 const idx_t qid, 
									 const idx_t sid,
									 const char qstrand,
									 const char sstrand,
									 const idx_t qsize,
									 const idx_t ssize,
									 const idx_t q_off_in_aln, 
									 const idx_t s_off_in_aln, 
									 const idx_t q_ext,
									 const idx_t s_ext,
									 M4Record& m4);

struct CandidateStartPosition
{
	idx_t qoff;
	idx_t toff;
	idx_t tstart;
	idx_t tsize;
	idx_t tid;
	int left_q, left_t;
	int right_q, right_t;
	int num1, num2;
	int score;
	idx_t toff_in_aln;
	char chain;
};

void print_candidate(const CandidateStartPosition& csp);

int Align(const char* query, const int q_len, const char* target, const int t_len, 
          const int band_tolerance, const int get_aln_str, Alignment* align, 
		  int* V, int* U, DPathData2* d_path, PathPoint* aln_path, const int right_extend);

void dw_in_one_direction(const char* query, const int query_size, const char* target, const int target_size,
						 int* U, int* V, Alignment* align, DPathData2* d_path, PathPoint* aln_path, 
						 SW_Parameters* swp, OutputStore* result, const int right_extend);

int  dw(const char* query, const int query_size, const int query_start,
        const char* target, const int target_size, const int target_start,
        int* U, int* V, Alignment* align, DPathData2* d_path, 
        PathPoint* aln_path, OutputStore* result, SW_Parameters* swp,
	    double error_rate, const int min_aln_size);

bool GetAlignment(const char* query, const int query_start, const int query_size,
				  const char* target, const int target_start, const int target_size,
				  DiffRunningData* drd, M5Record& m5, double error_rate, 
				  const int min_aln_size);

} // end of namespace ns_banded_sw

#endif  // DW_H
