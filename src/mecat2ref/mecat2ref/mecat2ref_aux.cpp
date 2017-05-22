#include "mecat2ref_aux.h"
#include "../common/defs.h"
#include <algorithm>
using namespace std;

int find_location(int *t_loc,int *t_seedn,int *t_score,long *loc,int k,int *rep_loc,float len,int read_len1, double ddfs_cutoff)
{
    int i,j,maxval=0,maxi,rep=0,lasti=0;
    for(i=0; i<k; i++)t_score[i]=0;
    for(i=0; i<k-1; i++)for(j=i+1; j<k; j++)if(t_seedn[j]-t_seedn[i]>0&&t_loc[j]-t_loc[i]>0&&t_loc[j]-t_loc[i]<read_len1&&fabs((t_loc[j]-t_loc[i])/((t_seedn[j]-t_seedn[i])*len)-1)<ddfs_cutoff)
            {
                t_score[i]++;
                t_score[j]++;
            }

    for(i=0; i<k; i++)
    {
        if(maxval<t_score[i])
        {
            maxval=t_score[i];
            maxi=i;
            rep=0;
        }
        else if(maxval==t_score[i])
        {
            rep++;
            lasti=i;
        }
    }
    for(i=0; i<4; i++)loc[i]=0;
    if(maxval>=5&&rep==maxval)
    {
        loc[0]=t_loc[maxi],loc[1]=t_seedn[maxi];
        *rep_loc=maxi;
        loc[2]=t_loc[lasti],loc[3]=t_seedn[lasti];
        return(1);
    }
    else if(maxval>=5&&rep!=maxval)
    {
        for(j=0; j<maxi; j++)if(t_seedn[maxi]-t_seedn[j]>0&&t_loc[maxi]-t_loc[j]>0&&t_loc[maxi]-t_loc[j]<read_len1&&fabs((t_loc[maxi]-t_loc[j])/((t_seedn[maxi]-t_seedn[j])*len)-1)<ddfs_cutoff)
            {
                if(loc[0]==0)
                {
                    loc[0]=t_loc[j];
                    loc[1]=t_seedn[j];
                    *rep_loc=j;
                }
                else
                {
                    loc[2]=t_loc[j];
                    loc[3]=t_seedn[j];
                }
            }
        j=maxi;
        if(loc[0]==0)
        {
            loc[0]=t_loc[j];
            loc[1]=t_seedn[j];
            *rep_loc=j;
        }
        else
        {
            loc[2]=t_loc[j];
            loc[3]=t_seedn[j];
        }
        for(j=maxi+1; j<k; j++)if(t_seedn[j]-t_seedn[maxi]>0&&t_loc[j]-t_loc[maxi]>0&&t_loc[j]-t_loc[maxi]<=read_len1&&fabs((t_loc[j]-t_loc[maxi])/((t_seedn[j]-t_seedn[maxi])*len)-1)<ddfs_cutoff)
            {
                if(loc[0]==0)
                {
                    loc[0]=t_loc[j];
                    loc[1]=t_seedn[j];
                    *rep_loc=j;
                }
                else
                {
                    loc[2]=t_loc[j];
                    loc[3]=t_seedn[j];
                }
            }
        return(1);
    }
    else return(0);
}

void
extract_sequences(const char* raw_ref,
				  const char* raw_read,
				  const long ref_start,
				  const int read_start,
				  const int read_size,
				  const long ref_size,
				  long& left_ref_size,
				  long& right_ref_size,
				  vector<char>& qstr,
				  vector<char>& tstr)
{
	long L1 = read_start, R1 = read_size - read_start;
	long L2 = ref_start, R2= ref_size - ref_start;
	long L = min(L1, L2), R = min(R1, R2);
	left_ref_size = min(L2, (long)(L * 1.2));
	right_ref_size = min(R2, (long)(R * 1.2));
	const u1_t* et = get_dna_encode_table();
	
	const char* rs = raw_ref + ref_start - left_ref_size;
	const char* re = raw_ref + ref_start + right_ref_size;
	tstr.clear();
	for (const char* r = rs; r != re; ++r) {
		u1_t c = *r;
		c = et[c];
		if (c > 3) c = 0;
		tstr.push_back((char)c);
	}
	
	qstr.clear();
	for (int i = 0; i < read_size; ++i) {
		u1_t c = raw_read[i];
		c = et[c];
		if (c > 3) c = 0;
		qstr.push_back((char)c);
	}
}

bool extend_candidate(candidate_save& can,
					  GapAligner* aligner,
					  const char* raw_ref,
					  const long ref_size,
					  const char* fwd_raw_read,
					  const char* rev_raw_read,
					  vector<char>& qstr,
					  vector<char>& tstr,
					  int read_name,
					  int read_len,
					  AlignInfo* alns,
					  int* naln,
					  TempResult* trv,
					  int& ntr)
{
	const char* raw_read = (can.chain == 'F') ? fwd_raw_read : rev_raw_read;
	int read_start = can.loc2;
	long ref_start = can.loc1 - 1;
	long left_ref_size, right_ref_size;
	extract_sequences(raw_ref, 
					  raw_read, 
					  ref_start, 
					  read_start,
					  read_len, 
					  ref_size, 
					  left_ref_size, 
					  right_ref_size, 
					  qstr, 
					  tstr);
	if (aligner->go(qstr.data(), read_start, read_len, tstr.data(), left_ref_size, tstr.size(), 1000)) {

		TempResult& r = trv[ntr++];
		r.read_id = read_name;
		r.read_dir = can.chain;
		r.vscore = can.score;
		r.qb = aligner->query_start();
		r.qe = aligner->query_end();
		r.qs = read_len;
		r.sb = ref_start - left_ref_size + aligner->target_start();
		r.se = ref_start - left_ref_size + aligner->target_end();
		strcpy(r.qmap, aligner->query_mapped_string());
		strcpy(r.smap, aligner->target_mapped_string());

		if (alns) {
			alns[*naln].qoff = r.qb;
			alns[*naln].qend = r.qe;
			alns[*naln].qdir = r.read_dir;
			alns[*naln].soff = r.sb;
			alns[*naln].send = r.se;
			alns[*naln].valid = 1;
			alns[*naln].id = ntr - 1;
			alns[*naln].prev_id = -1;
			alns[*naln].next_id = -1;
			alns[*naln].parent_id = -1;
			++(*naln);
		} 
		
		return true;
	}
	return false;
}

bool
fill_clipped_candidate(Back_List* block,
					   long bid,
					   candidate_save& can,
					   char chain,
					   int read_size,
					   int BC,
					   int block_size,
					   double ddfs_cutoff)
{
	int seedn[SM], boff[SM], score[SM], rep_loc;
	long locations[4];
	int n = min((int)block->score2, SM);
	for (int i = 0; i < n; ++i) {
		seedn[i] = block->seedno[i];
		boff[i] = block->loczhi[i];
		score[i] = 0;
	}
	int r = find_location(boff, seedn, score, locations, n, &rep_loc, BC, read_size, ddfs_cutoff);
	if (r) {
		can.score = score[rep_loc];
		can.chain = chain;
		can.loc1 = bid * block_size + locations[0];
		can.loc2 = (locations[1] - 1) * BC;
		return true;
	}
	return false;
}

bool
find_left_clipped_candidate(AlignInfo& aln,
							candidate_save& can,
						    Back_List* database,
						    int block_size,
						    int read_size,
						    int BC,
						    double ddfs_cutoff)
{
	if (aln.qoff <= CLIPPED || aln.soff <= CLIPPED) return false;
	int n1 = aln.qoff / block_size;
	int n2 = aln.soff / block_size;
	int n = min(n1, n2);
	int max_score = 0;
	Back_List* block = NULL;
	long bid = -1;
	for (--n2; n >= 0 && n2 >= 0; --n, --n2) {
		if (database[n2].score2 > max_score) {
			max_score = database[n2].score2;
			block = database + n2;
			bid = n2;
		}
	}
	bool ret = false;
	if (block && block->score2 > 4) {
		ret = fill_clipped_candidate(block, bid, can, aln.qdir, read_size, BC, block_size, ddfs_cutoff);
	}
	return ret;
}

bool
find_right_clipped_candidate(AlignInfo& aln,
							 candidate_save& can,
							 Back_List* database,
							 int block_size,
							 int read_size,
							 long ref_size,
							 int BC,
							 double ddfs_cutoff)
{
	if (read_size - aln.qend <= CLIPPED || ref_size - aln.send <= CLIPPED) return false;
	int n1 = (read_size - aln.qend) / block_size;
	int n2 = (ref_size - aln.send) / block_size;
	int n = min(n1, n2);
	int max_score = 0;
	long bid = -1;
	Back_List* block = NULL;
	int k = aln.send / block_size + 1;
	for (; n >= 0; --n, ++k) {
		if (database[k].score2 > max_score) {
			max_score = database[k].score2;
			block = database + k;
			bid = k;
		}
	}
	bool ret = false;
	if (block && block->score2 > 4) {
		ret = fill_clipped_candidate(block, bid, can, aln.qdir, read_size, BC, block_size, ddfs_cutoff);
	}
	return ret;
}

inline bool
is_left_clipped_align(const AlignInfo& a, const AlignInfo& b)
{
	if (a.qdir != b.qdir) return false;
	
	if (abs(b.qend - a.qoff) <= 200) {
		if (a.soff - b.send > -200 && a.soff - b.send < 10000) return true;
	}
	
	if (abs(b.send - a.soff) <= 200) {
		if (a.qoff - b.qend > -200 && a.qoff - b.qend < 10000) return true;
	}
	
	return false;
}

inline bool
is_right_clipped_align(const AlignInfo& a, const AlignInfo& b)
{
	if (a.qdir != b.qdir) return false;
	
	if (abs(a.qend - b.qoff) <= 200) {
		if (b.soff - a.send > -200 && b.soff - a.send < 10000) return true;
	}
	
	if (abs(a.send - b.soff) <= 200) {
		if (b.qoff - a.qend > -200 && b.qoff - a.qend < 10000) return true;
	}
	
	return false;
}

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
					  vector<char>& qstr,
					  vector<char>& tstr,
					  int read_name,
					  int read_len,
					  int block_size,
					  int BC,
					  Back_List* fwd_database, 
					  Back_List* rev_database,
					  double ddfs_cutoff)
{
	sort(alnv, alnv + naln);
	for (int i = 0; i < naln - 1; ++i) {
		if (!alnv[i].valid) continue;
		for (int j = i + 1; j < naln; ++j) {
			if (!alnv[j].valid) continue;
			if (AlignInfoContained(alnv[i], alnv[j])) alnv[j].valid = 0;
		}
	}
	int k = 0;
	for (int i = 0; i < naln; ++i) {
		if (alnv[i].valid) alnv[k++] = alnv[i];
	}
	naln = k;
	if (is_full_align(alnv[0], read_len)) return;
			
	for (int i = 0; i < naln - 1; ++i) {
		if (alnv[i].parent_id != -1) continue;
		for (int j = i + 1; j < naln; ++j) {
			if (alnv[j].parent_id != -1) continue;
			
			if (alnv[i].prev_id != -1 && is_left_clipped_align(alnv[i], alnv[j])) {
				alnv[i].prev_id = alnv[j].id;
				alnv[j].parent_id = alnv[i].id;
			}
			if (alnv[i].next_id != -1 && is_right_clipped_align(alnv[i], alnv[j])) {
				alnv[i].next_id = alnv[j].id;
				alnv[j].parent_id = alnv[i].id;
			}
		}
	}
	
	int n = min(naln, 3);
	k = 0;
	candidate_save can;
	for (int i = 0; i < n; ++i) {
		if (alnv[i].parent_id != -1) continue;
		Back_List* database = (alnv[i].qdir == 'F') ? fwd_database : rev_database;
		if (alnv[i].prev_id == -1 && find_left_clipped_candidate(alnv[i], can, database, block_size, read_len, BC, ddfs_cutoff)) {
			bool r = extend_candidate(can, 
							 aligner, 
							 raw_ref, 
							 ref_size,
							 fwd_raw_read, 
							 rev_raw_read, 
							 qstr, 
							 tstr, 
							 read_name, 
							 read_len, 
							 NULL, 
							 NULL, 
							 results, 
							 nresults);
			if (r) {
				TempResult& rslt = results[nresults - 1];
				AlignInfo& ai = alnv[naln + k];
				ai.qoff = rslt.qb;
				ai.qend = rslt.qe;
				ai.qdir = rslt.read_dir;
				ai.soff = rslt.sb;
				ai.send = rslt.se;
				ai.valid = 1;
				ai.id = nresults - 1;
				ai.prev_id = -1;
				ai.next_id = -1;
				ai.parent_id = -1;
				
				if (is_left_clipped_align(alnv[i], ai)) {
					ai.parent_id = alnv[i].id;
					alnv[i].prev_id = ai.id;
					++k;
				}
			}
		}
		if (alnv[i].next_id == -1 && find_right_clipped_candidate(alnv[i], can, database, block_size, read_len, ref_size, BC, ddfs_cutoff)) {
			bool r = extend_candidate(can, 
									  aligner, 
									  raw_ref, 
									  ref_size,
									  fwd_raw_read, 
									  rev_raw_read, 
									  qstr, 
									  tstr, 
									  read_name, 
									  read_len, 
									  NULL, 
									  NULL, 
									  results, 
									  nresults);
			if (r) {
				TempResult& rslt = results[nresults - 1];
				AlignInfo& ai = alnv[naln + k];
				ai.qoff = rslt.qb;
				ai.qend = rslt.qe;
				ai.qdir = rslt.read_dir;
				ai.soff = rslt.sb;
				ai.send = rslt.se;
				ai.valid = 1;
				ai.id = nresults - 1;
				ai.prev_id = -1;
				ai.next_id = -1;
				ai.parent_id = -1;

				if (is_right_clipped_align(alnv[i], ai)) {
					ai.parent_id = alnv[i].id;
					alnv[i].next_id = ai.id;
					++k;
				}
			}
		}
	}
	if (!k) return;
	
	naln += k;
	sort(alnv, alnv + naln);
	k = 0;
	for (int i = 0; i < naln; ++i) {
		if (is_full_align(alnv[i], read_len)) {
			alnv[i].parent_id = -1;
			alnv[i].prev_id = -1;
			alnv[i].next_id = -1;
			++k;
		}
	}
	if (k) naln = k;
}

void
output_results(AlignInfo* alnv,
			   int& naln,
			   TempResult* results,
			   int& nresults,
			   int num_output,
			   FILE* out)
{
	int n = 0;
	for (int i = 0; i < naln && n < num_output; ++i) {
		if (alnv[i].parent_id != -1) continue;
		int id = alnv[i].id;
		output_temp_result(results + id, out);
		id = alnv[i].prev_id;
		if (id != -1) output_temp_result(results + id, out);
		id = alnv[i].next_id;
		if (id != -1) output_temp_result(results + id, out);
		++n;
	}
}
