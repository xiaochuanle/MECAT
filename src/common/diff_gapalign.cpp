#include "diff_gapalign.h"
//#include "smart_assert.h"

#include <algorithm>
#include <set>

#include <cstdlib>
#include <cstring>

using namespace std;

struct SCompareDPathData2
{
    bool operator()(const DPathData2& a, const DPathData2& b) {
        return (a.d == b.d) ? (a.k < b.k) : (a.d < b.d);
    }
};

struct DPathDataExtractor
{
	DPathData2* operator()(const int d, const int k)
	{
		target.d = d;
		target.k = k;
		return std::lower_bound(d_path_list, d_path_list + d_path_list_size, target, scmp);
	}
	
	DPathDataExtractor(DPathData2* d_path, const int n)
		: d_path_list(d_path),
		  d_path_list_size(n) {
		  }
	
    DPathData2 target;
	SCompareDPathData2 scmp;
	DPathData2* d_path_list;
	const int d_path_list_size;
};

void
GetAlignString(const char* query, const int q_len, 
			   const char* target, const int t_len,
			   DPathData2* d_path, const int d_path_size,
			   PathPoint* aln_path, Alignment* align,
			   const int qe, const int te, 
			   int d, int k,
			   const int right_extend)
{
	int cd = d;
	int ck = k;
	int aln_path_idx = 0;
	int i;
	DPathDataExtractor dpath_extrator(d_path, d_path_size);
	while (cd >= 0 && aln_path_idx < q_len + t_len + 1)
	{
		DPathData2* d_path_aux = dpath_extrator(cd, ck);
		aln_path[aln_path_idx].x = d_path_aux->x2;
		aln_path[aln_path_idx].y = d_path_aux->y2;
		++aln_path_idx;
		aln_path[aln_path_idx].x = d_path_aux->x1;
		aln_path[aln_path_idx].y = d_path_aux->y1;
		++aln_path_idx;
		ck = d_path_aux->pre_k;
		cd -= 1;
	}
	--aln_path_idx;
	int cx = aln_path[aln_path_idx].x;
	int cy = aln_path[aln_path_idx].y;
	align->aln_q_s = cx;
	align->aln_t_s = cy;
	int aln_pos = 0;
	while (aln_path_idx > 0)
	{
		--aln_path_idx;
		int nx = aln_path[aln_path_idx].x;
		int ny = aln_path[aln_path_idx].y;
		if (cx == nx && cy == ny) continue;
		
		if (cx == nx && cy != ny) {
			for (i = 0; i < ny - cy; ++i) {
				align->q_aln_str[aln_pos + i] = GAP_CODE;
				align->t_aln_str[aln_pos + i] = extract_char<char>(target, cy + i, right_extend);
			}
			aln_pos += ny - cy;
		} else if (cx != nx && cy == ny) {
			for (i = 0; i < nx - cx; ++i) {
				align->q_aln_str[aln_pos + i] = extract_char<char>(query, cx + i, right_extend);
				align->t_aln_str[aln_pos + i] = GAP_CODE;
			}
			aln_pos += nx - cx;
		} else {
			for (i = 0; i < nx - cx; ++i) {
				align->q_aln_str[aln_pos + i] = extract_char<char>(query, cx + i, right_extend);
			}
			for (i = 0; i < ny - cy; ++i) {
				align->t_aln_str[aln_pos + i] = extract_char<char>(target, cy + i, right_extend);
			}
			aln_pos += ny - cy;
		}
		
		cx = nx;
		cy = ny;
	}
	align->aln_str_size = aln_pos;
}
			   

int Align(const char* query, const int q_len, const char* target, const int t_len,
		  const int band_tolerance, const int get_aln_str, Alignment* align,
		  int* V, int* U, DPathData2* d_path, PathPoint* aln_path, const int right_extend)
{
	int k_offset;
	int  d;
	int  k, k2;
	int best_m;
	int min_k, new_min_k, max_k, new_max_k, pre_k;
	int x = -1, y = -1;
	int max_d, band_size;
	unsigned long d_path_idx = 0, max_idx = 0;
	int aligned = 0;
	int best_x = -1, best_y = -1, best_d = q_len + t_len + 100, best_k = 0, best_d_path_idx = -1;

	max_d = (int)(.3 * (q_len + t_len));
	k_offset = max_d;
	band_size = band_tolerance * 2;
	align->init();
	best_m = -1;
	min_k = 0;
	max_k = 0;
	d_path_idx = 0;
	max_idx = 0;

	for (d = 0; d < max_d; ++d)
	{
		if (max_k - min_k > band_size) break;
		
		for (k = min_k; k <= max_k; k += 2)
		{
			if( k == min_k || (k != max_k && V[k - 1 + k_offset] < V[k + 1 + k_offset]) )
			{ pre_k = k + 1; x = V[k + 1 + k_offset]; }
			else
			{ pre_k = k - 1; x = V[k - 1 + k_offset] + 1; }
			y = x - k;
			d_path[d_path_idx].d = d;
			d_path[d_path_idx].k = k;
			d_path[d_path_idx].x1 = x;
			d_path[d_path_idx].y1 = y;

			if (right_extend)
				while( x < q_len && y < t_len && query[x] == target[y]) { ++x; ++y; }
			else
				while( x < q_len && y < t_len && query[-x] == target[-y]) { ++x; ++y; }

			d_path[d_path_idx].x2 = x;
			d_path[d_path_idx].y2 = y;
			d_path[d_path_idx].pre_k = pre_k;
			++d_path_idx;

			V[k + k_offset] = x;
			U[k + k_offset] = x + y;
			if (x + y > best_m) {
				best_m = x + y;
				best_x = x;
				best_y = y;
				best_d = d;
				best_k = k;
				best_d_path_idx = d_path_idx;
			}
			if (x >= q_len || y >= t_len)
			{ aligned = 1; max_idx = d_path_idx; break; }
		}

		// for banding
		new_min_k = max_k;
		new_max_k = min_k;
		for (k2 = min_k; k2 <= max_k; k2 += 2)
			if (U[k2 + k_offset] >= best_m - band_tolerance)
			{ new_min_k = std::min(new_min_k, k2); new_max_k = std::max(new_max_k, k2); }
		max_k = new_max_k + 1;
		min_k = new_min_k - 1;

		if (aligned)
		{
			align->aln_q_e = x;
			align->aln_t_e = y;
			align->dist = d;
			align->aln_str_size = (x + y + d) / 2;
			align->aln_q_s = 0;
			align->aln_t_s = 0;

			if (get_aln_str) {
				GetAlignString(query, q_len, target, t_len, d_path, max_idx, aln_path, align, x, y, d, k, right_extend);
			}
			break;
		} 
	}
	
	if (!aligned) {
		if (best_x > 0) {
			align->aln_q_e = best_x;
			align->aln_t_e = best_y;
			align->dist = best_d;
			align->aln_str_size = (best_x + best_y + best_d) / 2;
			align->aln_q_s = 0;
			align->aln_t_s = 0;
			if (get_aln_str) {
				GetAlignString(query, q_len, target, t_len, d_path, best_d_path_idx, aln_path, align, best_x, best_y, best_d, best_k, right_extend);
			}
		} else {
			align->aln_q_e = 0;
			align->aln_t_e = 0;
			align->dist = 0;
			align->aln_str_size = 0;
			align->aln_q_s = 0;
			align->aln_t_s = 0;
		}
	}
	
	return (align->aln_q_e == q_len || align->aln_t_e == t_len);
}

void dw_in_one_direction(const char* query, const int query_size, const char* target, const int target_size,
						 int* U, int* V, Alignment* align, DPathData2* d_path, PathPoint* aln_path,
						 DiffAlignParameters* swp, OutputStore* result, const int right_extend)
{
	const int kTailMatchBP = 4;
	const int kBlkSize = swp->segment_size;
	int qidx = 0, tidx = 0;
	int qblk, tblk;
	const char* seq1;
	const char* seq2;
	while (1) {
		fill(U, U + swp->row_size, 0);
		fill(V, V + swp->column_size, 0);
		bool last_block = retrieve_next_aln_block(query,
												  qidx,
												  query_size,
												  target,
												  tidx,
												  target_size,
												  kBlkSize,
												  right_extend,
												  seq1,
												  seq2,
												  qblk,
												  tblk);	  
		Align(seq1,
			  qblk,
			  seq2,
			  tblk,
			  0.3 * max(qblk, tblk),
			  400,
			  align,
			  U, 
			  V,
			  d_path,
			  aln_path,
			  right_extend);
		
		int qcnt = 0, tcnt = 0, acnt = 0;
		const bool trim = trim_mismatch_end(align->q_aln_str, 
											align->t_aln_str, 
											align->aln_str_size, 
											kTailMatchBP, 
											qcnt, 
											tcnt, 
											acnt);
		if (!trim) break;
		
		bool full_map = false;
		if (qblk - align->aln_q_e <= 20 || tblk - align->aln_t_e <= 20) full_map = true;

		if (last_block || (!full_map)) {
			qcnt -= kTailMatchBP;
			tcnt -= kTailMatchBP;
			acnt -=kTailMatchBP;
		}
		align->aln_str_size -= acnt;
		if (right_extend) {
			memcpy(result->right_store1 + result->right_store_size, align->q_aln_str, align->aln_str_size);
			memcpy(result->right_store2 + result->right_store_size, align->t_aln_str, align->aln_str_size);
			result->right_store_size += align->aln_str_size;
		} else {
			memcpy(result->left_store1 + result->left_store_size, align->q_aln_str, align->aln_str_size);
			memcpy(result->left_store2 + result->left_store_size, align->t_aln_str, align->aln_str_size);
			result->left_store_size += align->aln_str_size;
		}
		
		if (last_block || (!full_map)) break;
		qidx += (align->aln_q_e - qcnt);
		tidx += (align->aln_t_e - tcnt);
	}
}

bool
DiffAligner::go(const char* query, const int qstart, const int qsize, 
				const char* target, const int tstart, const int tsize,
				const int min_aln_size)
{
	result->init();
	align->init();
	dw_in_one_direction(query + qstart - 1, qstart, 
						target + tstart - 1, tstart,
						dynq, dynt, align, d_path, 
						aln_path, &param, result, 0);
	dw_in_one_direction(query + qstart, qsize - qstart,
						target + tstart, tsize - tstart,
						dynq, dynt, align, d_path,
						aln_path, &param, result, 1);
	int i, j, k, idx = 0;
	const char* dt = "ACGT-";
	for (k = result->left_store_size - 1, i = 0, j = 0; k >= 0; --k, ++idx) {
		unsigned char ch = result->left_store1[k];
		d_assert(ch >= 0 && ch <= 4);
		ch = dt[ch];
		result->out_store1[idx] = ch;
		if (ch != '-') ++i;
		
		ch = result->left_store2[k];
		d_assert(ch >= 0 && ch <= 4);
		ch = dt[ch];
		result->out_store2[idx] = ch;
		if (ch != '-') ++j;
	}
	result->query_start = qstart - i;
	result->target_start = tstart - j;
	d_assert(result->query_start >= 0);
	d_assert(result->target_start >= 0);
	
	for (k = 0, i = 0, j = 0; k < result->right_store_size; ++k, ++idx) {
		unsigned char ch = result->right_store1[k];
		d_assert(ch >= 0 && ch <= 4);
		ch = dt[ch];
		result->out_store1[idx] = ch;
		if (ch != '-') ++i;
		
		ch = result->right_store2[k];
		d_assert(ch >= 0 && ch <= 4);
		ch = dt[ch];
		result->out_store2[idx] = ch;
		if (ch != '-') ++j;
	}
	result->out_store_size = idx;
	result->query_end = qstart + i;
	result->target_end = tstart + j;
	result->out_store1[result->out_store_size] = '\0';
	result->out_store2[result->out_store_size] = '\0';
	
	return result->out_store_size >= min_aln_size;
}
