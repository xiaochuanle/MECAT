#ifndef AUX_H
#define AUX_H

#include "../common/gapalign.h"
#include "output.h"
#include "mecat2ref_defs.h"
#include <vector>

struct AlignInfo
{
	int qid, qoff, qend;
	int parent_id, id, prev_id, next_id;
	char valid;
	char qdir;
	long soff, send;
	
	bool operator < (const AlignInfo& rhs) const {
		return (qend - qoff) > (rhs.qend - rhs.qoff);
	}
};

inline bool
is_full_align(const AlignInfo& ai, const int qsize)
{
	return ai.qend - ai.qoff >= qsize * 0.9;
}

inline bool
AlignInfoContained(const AlignInfo& a, const AlignInfo& b)
{
	const int extra = 100;
	bool r = (a.qdir == b.qdir)
			 &&
			 (b.qoff + extra >= a.qoff)
			 &&
			 (b.qend <= a.qend + extra)
			 &&
			 (b.soff + extra >= a.soff)
			 &&
			 (b.send <= a.send + extra);
	return r;
}

int 
find_location(int *t_loc,int *t_seedn,int *t_score,long *loc,int k,int *rep_loc,float len,int read_len1, double ddfs_cutoff);

bool extend_candidate(candidate_save& can,
					  GapAligner* aligner,
					  const char* raw_ref,
					  const long ref_size,
					  const char* fwd_raw_read,
					  const char* rev_raw_read,
					  std::vector<char>& qstr,
					  std::vector<char>& tstr,
					  int read_name,
					  int read_len,
					  AlignInfo* alns,
					  int* naln,
					  TempResult* trv,
					  int& ntr);

void
rescue_clipped_align(AlignInfo* alnv,
					 int& naln,
					 TempResult* results,
					 int& nresults,
					 GapAligner* aligner,
					 const char* raw_ref,
					 const long ref_size,
					 const char* fwd_raw_read,
					 const char* rev_raw_read,
					 std::vector<char>& qstr,
					 std::vector<char>& tstr,
					 int read_name,
					 int read_len,
					 int block_size,
					 int BC,
					 Back_List* fwd_database, 
					 Back_List* rev_database,
					 double ddfs_cutoff);

void
output_results(AlignInfo* alnv,
			   int& naln,
			   TempResult* results,
			   int& nresults,
			   int num_output,
			   FILE* out);

#define CLIPPED 2000

#endif // AUX_H
