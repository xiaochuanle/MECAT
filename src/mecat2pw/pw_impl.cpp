#include "../common/split_database.h"
#include "pw_options.h"
#include "../common/diff_gapalign.h"
#include "../common/xdrop_gapalign.h"
#include "../common/packed_db.h"
#include "../common/lookup_table.h"
#include "pw_impl.h"

#include <algorithm>

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#define MSS MAX_SEQ_SIZE

static int MAXC = 100;
static int output_gapped_start_point = 1;
static int kmer_size = 13;
static const double ddfs_cutoff_pacbio = 0.25;
static const double ddfs_cutoff_nanopore = 0.25;
static double ddfs_cutoff = ddfs_cutoff_pacbio;
static int min_align_size = 0;
static int min_kmer_match = 0;
static int min_kmer_dist = 0;

using namespace std;

PWThreadData::PWThreadData(options_t* opt, volume_t* ref, volume_t* rd, ref_index* idx, std::ostream* o)
	: options(opt), used_thread_id(0), reference(ref), reads(rd), ridx(idx), out(o), m4_results(NULL), ec_results(NULL), next_processed_id(0)
{
	pthread_mutex_init(&id_lock, NULL);
	if (options->task == TASK_SEED)
	{
		safe_malloc(ec_results, ExtensionCandidate*, options->num_threads);
		for (int i = 0; i < options->num_threads; ++i)
			safe_malloc(ec_results[i], ExtensionCandidate, kResultListSize);
	}
	else if (options->task == TASK_ALN)
	{
		safe_malloc(m4_results, M4Record*, options->num_threads);
		for (int i = 0; i < options->num_threads; ++i)
			safe_malloc(m4_results[i], M4Record, kResultListSize);
	}
	else
	{
		LOG(stderr, "Task must be either %d or %d, not %d!", TASK_SEED, TASK_ALN, options->task);
		abort();
	}
	pthread_mutex_init(&result_write_lock, NULL);
	pthread_mutex_init(&read_retrieve_lock, NULL);
}

PWThreadData::~PWThreadData()
{
	if (ec_results)
	{
		for (int i = 0; i < options->num_threads; ++i) safe_free(ec_results[i]);
		safe_free(ec_results);
	}
	if (m4_results)
	{
		for (int i = 0; i < options->num_threads; ++i) safe_free(m4_results[i]);
		safe_free(m4_results);
	}
}

void
reverse_complement(char* dst, const char* src, const int size)
{
	const uint8_t* rc_table = get_dna_complement_table();
	int i;
	for (i = 0; i < size; ++i)
	{
		uint8_t c = src[i];
		assert(c >= 0 && c < 4);
		c = rc_table[c];
		dst[size - 1 - i] = c;
	}
}

int
extract_kmers(const char* s, const int ssize, int* kmer_ids)
{
	int max_id = 1 << (2 * kmer_size);
	int num_kmers = (ssize - kmer_size) / BC + 1;
	int i, j;
	for (i = 0; i < num_kmers; ++i)
	{
		int eit = 0, start = i * BC;
		for (j = 0; j < kmer_size; ++j) eit = (eit << 2) | s[start + j];
		kmer_ids[i] = eit;
		assert(eit < max_id);
	}
	return num_kmers;
}

SeedingBK::SeedingBK(const int ref_size)
{
	const int num_segs = ref_size / ZV + 5;
	safe_malloc(index_list, int, num_segs);
	safe_malloc(index_score, short, num_segs);
	safe_malloc(database, Back_List, num_segs);
	safe_malloc(kmer_ids, int, MAX_SEQ_SIZE);
	for (int i = 0; i < num_segs; ++i) 
	{
		database[i].score = 0;
		database[i].index = -1;
	}
}

SeedingBK::~SeedingBK()
{
	safe_free(index_list);
	safe_free(index_score);
	safe_free(database);
	safe_free(kmer_ids);
}

void insert_loc(Back_List *spr,int loc,int seedn,float len)
{
    int list_loc[SI],list_score[SI],list_seed[SI],i,j,minval,mini;
    for(i=0; i<SM; i++)
    {
        list_loc[i]=spr->loczhi[i];
        list_seed[i]=spr->seedno[i];
        list_score[i]=0;
    }
    list_loc[SM]=loc;
    list_seed[SM]=seedn;
    list_score[SM]=0;
    mini=-1;
    minval=10000;
    for(i=0; i<SM; i++)for(j=i+1; j<SI; j++)if(list_seed[j]-list_seed[i]>0&&list_loc[j]-list_loc[i]>0&&fabs((list_loc[j]-list_loc[i])/((list_seed[j]-list_seed[i])*len)-1.0)<ddfs_cutoff)
            {
                list_score[i]++;
                list_score[j]++;
            }
    for(i=0; i<SI; i++)if(minval>list_score[i])
        {
            minval=list_score[i];
            mini=i;
        }
    if(minval==SM)
    {
        spr->loczhi[SM-1]=loc;
        spr->seedno[SM-1]=seedn;
    }
    else if(minval<SM&&mini<SM)
    {
        for(i=mini; i<SM; i++)
        {
            spr->loczhi[i]=list_loc[i+1];
            spr->seedno[i]=list_seed[i+1];
        }
        spr->score--;
    }
}

int find_location(int *t_loc,int *t_seedn,int *t_score,int *loc,int k,int *rep_loc,float len,int read_len1)
{
    int i,j,maxval=0,maxi,rep=0,lasti = 0,tempi;
    for(i=0; i<k; i++)t_score[i]=0;
    for(i=0; i<k-1; i++)for(j=i+1,tempi=t_seedn[i]; j<k; j++)if(tempi!=t_seedn[j]&&t_seedn[j]-t_seedn[i]>0&&t_loc[j]-t_loc[i]>0&&t_loc[j]-t_loc[i]<read_len1&&fabs((t_loc[j]-t_loc[i])/((t_seedn[j]-t_seedn[i])*len)-1)<ddfs_cutoff)
            {
                t_score[i]++;
                t_score[j]++;
                tempi=t_seedn[j];
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

int
seeding(const char* read, const int read_size, ref_index* ridx, SeedingBK* sbk)
{
	int* kmer_ids = sbk->kmer_ids;
	int* index_list = sbk->index_list;
	int* index_spr = index_list;
	short* index_score = sbk->index_score;
	short* index_ss = index_score;
	Back_List* database = sbk->database;
	
	int num_kmers = extract_kmers(read, read_size, kmer_ids);
	int km;
	int used_segs = 0;
	for (km = 0; km < num_kmers; ++km)
	{
		int num_seeds = ridx->kmer_counts[kmer_ids[km]];
		int* seed_arr = ridx->kmer_starts[kmer_ids[km]];
		int sid;
		int endnum = 0;
		for (sid = 0; sid < num_seeds; ++sid)
		{
			int seg_id = seed_arr[sid] / ZV;
			int seg_off = seed_arr[sid] % ZV;
			Back_List* spr = database + seg_id;
			if (spr->score == 0 || spr->seednum < km + 1)
			{
				int loc = ++spr->score;
				if (loc <= SM) { spr->loczhi[loc - 1] = seg_off; spr->seedno[loc - 1] = km + 1; }
				else insert_loc(spr, seg_off, km + 1, BC);
				int s_k;
				if (seg_id > 0) s_k = spr->score + (spr - 1)->score;
				else s_k = spr->score;
				if (endnum < s_k) endnum = s_k;
				if (spr->index == -1)
				{
					*(index_spr++) = seg_id;
					*(index_ss++) = s_k;
					spr->index = used_segs++;
				}
				else index_score[spr->index] = s_k;
			}
			spr->seednum = km + 1;
		}
	}
	return used_segs;
}

int
get_candidates(volume_t* ref, 
			   SeedingBK* sbk, 
			   const int num_segs, 
			   const int read_id, 
			   const int read_size, 
			   const char chain,
			   candidate_save* candidates, 
			   int candidatenum)
{
	int* index_list = sbk->index_list;
	int* index_spr = index_list;
	short* index_score = sbk->index_score;
	short* index_ss = index_score;
	Back_List* database = sbk->database;
	const int temp_arr_size = 2 * SM + 10;
	int temp_list[temp_arr_size],temp_seedn[temp_arr_size],temp_score[temp_arr_size];
	candidate_save *candidate_loc = candidates, candidate_temp;
	int location_loc[4],repeat_loc;
	int i, j, k, u_k;
	for (i = 0; i < num_segs; ++i, ++index_spr, ++index_ss) 
		if (*index_ss >= 2 * min_kmer_match)
		{
			Back_List *spr = database + (*index_spr), *spr1;
			if (spr->score == 0) continue;
			int s_k = spr->score;
			int start_loc = *index_spr;
			start_loc = MUL_ZV(start_loc);
			int loc;
			if ((*index_spr) > 0)
			{
				loc = (spr - 1)->score;
				if (loc > 0) 
				{
					start_loc = (*index_spr - 1);
					start_loc = MUL_ZV(start_loc);
				}
			}
			else loc = 0;
			
			if (!loc)
				for (j = 0, u_k = 0; j < s_k && j < SM; ++j)
				{
					temp_list[u_k] = spr->loczhi[j];
					temp_seedn[u_k] = spr->seedno[j];
					++u_k;
				}
			else
			{
				k = loc;
				u_k = 0;
				spr1 = spr - 1;
				for (j = 0; j < k && j < SM; ++j)
				{
					temp_list[u_k] = spr1->loczhi[j];
					temp_seedn[u_k] = spr1->seedno[j];
					++u_k;
				}
				for (j = 0; j < s_k && j < SM; ++j)
				{
					temp_list[u_k] = spr->loczhi[j] + ZV;
					temp_seedn[u_k] = spr->seedno[j];
					++u_k;
				}
			}
			
			{
				int f = find_location(temp_list, temp_seedn, temp_score, location_loc, u_k, &repeat_loc, BC, read_size);
				if (!f) continue;
				if (temp_score[repeat_loc] < 2 * min_kmer_match + 2) continue;
			}
			
			candidate_temp.score = temp_score[repeat_loc];
			candidate_temp.chain = chain;
			int loc_seed = temp_seedn[repeat_loc];
			location_loc[0] = start_loc + location_loc[0];
			int loc_list = location_loc[0];
			int sid = get_read_id_from_offset_list(ref->offset_list, location_loc[0]);
			int sstart = ref->offset_list->offset_list[sid].offset;
			int ssize = ref->offset_list->offset_list[sid].size;
			int send = sstart + ssize + 1;
			sid += ref->start_read_id;
			if (sid > read_id) continue;
			if (sid == read_id)
			{
				u_k = DIV_ZV(sstart);
				spr = database + u_k;
				s_k = MOD_ZV(sstart);
				for (j = 0, k = 0; j < spr->score && j < SM; ++j)
					if (spr->loczhi[j] < s_k) { spr->loczhi[k] = spr->loczhi[j]; ++k; }
				spr->score = k;
				for (++spr, ++u_k, k = DIV_ZV(send); u_k < k; ++u_k, ++spr) spr->score = 0;
				for (j = 0, k = 0, s_k = MOD_ZV(send); j < spr->score && j < SM; ++j)
					if (spr->loczhi[j] > s_k) { spr->loczhi[k] = spr->loczhi[j]; ++k; }
				spr->score = k;
			}
			else
			{
				candidate_temp.readno = sid;
				candidate_temp.readstart = sstart;
				location_loc[1] = (location_loc[1] - 1) * BC;
				int left_length1 = location_loc[0] - sstart + kmer_size - 1;
				int right_length1 = send - location_loc[0];
				int left_length2 = location_loc[1] + kmer_size - 1;
				int right_length2 = read_size - location_loc[1];
				int num1 = (left_length1 > left_length2) ? left_length2 : left_length1;
				int num2 = (right_length1 > right_length2) ? right_length2 : right_length1;
				if (num1 + num2 < min_kmer_dist) continue;
				candidate_temp.loc1=location_loc[0] - sstart;
				candidate_temp.num1=num1;
				candidate_temp.loc2=location_loc[1];
				candidate_temp.num2=num2;
				candidate_temp.left1=left_length1;
				candidate_temp.left2=left_length2;
				candidate_temp.right1=right_length1;
				candidate_temp.right2=right_length2;
				
				int seedcount = 0;
				int nlb = (num1 + ZV - 1);
				nlb = DIV_ZV(nlb);
				for(u_k=*index_spr-1, spr1=spr-1; u_k>=0&&nlb>0; spr1--,--nlb,u_k--)if(spr1->score>0)
					{
						start_loc = MUL_ZV(u_k);
						int scnt = min((int)spr1->score, SM);
						for(j=0,s_k=0; j < scnt; j++)if(fabs((loc_list-start_loc-spr1->loczhi[j])/((loc_seed-spr1->seedno[j])*BC*1.0)-1.0)<ddfs_cutoff)
							{
								seedcount++;
								s_k++;
							}
						if(s_k*1.0 / scnt > 0.4)
						{
							spr1->score=0;
						}
					}
				//find all right seed
				int nrb = (num2 + ZV - 1);
				nrb = DIV_ZV(nrb);
				for(u_k=*index_spr+1,spr1=spr+1; nrb; spr1++,--nrb,u_k++)if(spr1->score>0)
					{
						start_loc = MUL_ZV(u_k);
						int scnt = min((int)spr1->score, SM);
						for(j=0,s_k=0; j < scnt; j++)if(fabs((start_loc+spr1->loczhi[j]-loc_list)/((spr1->seedno[j]-loc_seed)*BC*1.0)-1.0)<ddfs_cutoff)
							{
								seedcount++;
								s_k++;
							}
						if(s_k*1.0 / scnt > 0.4)
						{
							spr1->score=0;
						}
					}
				
				candidate_temp.score=candidate_temp.score+seedcount;
				//insert canidate position or delete this position
				int low=0;
				int high=candidatenum-1;
				int mid;
				while(low<=high)
				{
					mid=(low+high)/2;
					if(mid>=candidatenum||candidate_loc[mid].score<candidate_temp.score)high=mid-1;
					else low=mid+1;
				}
				if(candidatenum<MAXC)for(u_k=candidatenum-1; u_k>high; u_k--)candidate_loc[u_k+1]=candidate_loc[u_k];
				else for(u_k=candidatenum-2; u_k>high; u_k--)candidate_loc[u_k+1]=candidate_loc[u_k];
				if(high+1<MAXC)candidate_loc[high+1]=candidate_temp;
				if(candidatenum<MAXC)candidatenum++;
				else candidatenum=MAXC;
			}
		}
	
	for(i=0,index_spr=index_list; i<num_segs; i++,index_spr++)
	{
		database[*index_spr].score=0;
		database[*index_spr].index=-1;
	}
	return candidatenum;
}

void
fill_m4record(GapAligner* aligner, const int qid, const int sid,
			  const char qchain, int qsize, int ssize,
			  int qstart, int sstart, int vscore, M4Record* m)
{
	if (qchain == 'F')
	{
		m->qid = sid;
		m->sid = qid;
		m->ident = aligner->calc_ident();
		m->vscore = vscore;
		m->qdir = 0;
		m->qoff = aligner->target_start();
		m->qend = aligner->target_end();
		m->qsize = ssize;
		m->sdir = 0;
		m->soff = aligner->query_start();
		m->send = aligner->query_end();
		m->ssize = qsize;
		m->qext = sstart;
		m->sext = qstart;
	}
	else
	{
		m->qid = sid;
		m->sid = qid;
		m->ident = aligner->calc_ident();
		m->vscore = vscore;
		m->qdir = 0;
		m->qoff = aligner->target_start();
		m->qend = aligner->target_end();
		m->qsize = ssize;
		m->sdir = 1;
		m->soff = qsize - aligner->query_end();
		m->send = qsize - aligner->query_start();
		m->ssize = qsize;
		m->qext = sstart;
		m->sext = qsize - 1 - qstart;
	}
}


void output_m4record(ostream& out, const M4Record& m4)
{
	const char sep = '\t';
	
	out << m4qid(m4)    << sep
	    << m4sid(m4)    << sep
	    << m4ident(m4) << sep
	    << m4vscore(m4)   << sep
	    << m4qdir(m4)   << sep
	    << m4qoff(m4)   << sep
	    << m4qend(m4)   << sep
	    << m4qsize(m4)  << sep
	    << m4sdir(m4)   << sep
	    << m4soff(m4)   << sep
	    << m4send(m4)   << sep
	    << m4ssize(m4);
	
	if (output_gapped_start_point)
		out << sep
			<< m4qext(m4)	<< sep
			<< m4sext(m4);
	out << "\n";
}

void
print_m4record_list(ostream* out, M4Record* m4_list, int num_m4)
{
	for (int i = 0; i < num_m4; ++i) output_m4record(*out, m4_list[i]);
}

struct CmpM4RecordByQidAndOvlpSize
{
	bool operator()(const M4Record& a, const M4Record& b)
	{
		if (m4qid(a) != m4qid(b)) return m4qid(a) < m4qid(b);
		int o1 = M4RecordOverlapSize(a);
		int o2 = M4RecordOverlapSize(b);
		return o1 > o2;
	}
};

void
check_records_containment(M4Record* m4v, int s, int e, int* valid)
{
	const int soft = 100;
	
	for (int i = s; i < e; ++i)
	{
		if (!valid[i]) continue;
		int qb1 = m4qoff(m4v[i]);
		int qe1 = m4qend(m4v[i]);
		int sb1 = m4soff(m4v[i]);
		int se1 = m4send(m4v[i]);
		for (int j = i + 1; j < e; ++j)
		{
			if (!valid[j]) continue;
			if (m4sdir(m4v[i]) != m4sdir(m4v[j])) continue;
			int qb2 = m4qoff(m4v[j]);
			int qe2 = m4qend(m4v[j]);
			int sb2 = m4soff(m4v[j]);
			int se2 = m4send(m4v[j]);
			
			if (qb2 + soft >= qb1 && qe2 - soft <= qe1 && sb2 + soft >= sb1 && se2 - soft <= se1) valid[j] = 0;
		}
	}
}

void
append_m4v(M4Record* glist, int* glist_size,
		   M4Record* llist, int* llist_size,
		   ostream* out, pthread_mutex_t* results_write_lock)
{
	sort(llist, llist + *llist_size, CmpM4RecordByQidAndOvlpSize());
	int i = 0, j;
	int valid[*llist_size];
	fill(valid, valid + *llist_size, 1);
	while (i < *llist_size)
	{
		idx_t qid = m4qid(llist[i]);
		j = i + 1;
		while (j < *llist_size && m4qid(llist[j]) == qid) ++j;
		if (j - i > 1) check_records_containment(llist, i, j, valid);
		i = j;
	}
	
	if ((*glist_size) + (*llist_size) > PWThreadData::kResultListSize)
	{
		pthread_mutex_lock(results_write_lock);
		print_m4record_list(out, glist, *glist_size);
		*glist_size = 0;
		pthread_mutex_unlock(results_write_lock);
	}
	
	for (i = 0; i < *llist_size; ++i)
		if (valid[i])
		{
			glist[*glist_size] = llist[i];
			++(*glist_size);
		}
	
	*llist_size = 0;
}

inline void
get_next_chunk_reads(PWThreadData* data, int& Lid, int& Rid)
{
	pthread_mutex_lock(&data->read_retrieve_lock);
	Lid = data->next_processed_id;
	Rid = Lid + CHUNK_SIZE;
	if (Rid > data->reads->num_reads) Rid = data->reads->num_reads;
	data->next_processed_id += CHUNK_SIZE;
	pthread_mutex_unlock(&data->read_retrieve_lock);
}

void
pairwise_mapping(PWThreadData* data, int tid)
{
	char *read, *read1, *read2, *subject;
	safe_malloc(read1, char, MSS);
	safe_malloc(read2, char, MSS);
	safe_malloc(subject, char, MSS);
	SeedingBK* sbk = new SeedingBK(data->reference->curr);
	candidate_save candidates[MAXC];
	int num_candidates = 0;
	M4Record* m4_list = data->m4_results[tid];
	int m4_list_size = 0;
	M4Record* m4v = new M4Record[MAXC];
	int num_m4 = 0;
	GapAligner* aligner = NULL;
	if (data->options->tech == TECH_PACBIO) {
		aligner = new DiffAligner(0);
	} else if (data->options->tech == TECH_NANOPORE) {
		aligner = new XdropAligner(0);
	} else {
		ERROR("TECH must be either %d or %d", TECH_PACBIO, TECH_NANOPORE);
	}

	int rid, Lid, Rid;
	while (1)
	{
		get_next_chunk_reads(data, Lid, Rid);
		if (Lid >= data->reads->num_reads) break;
		for (rid = Lid; rid < Rid; ++rid)
		{
			int rsize = data->reads->offset_list->offset_list[rid].size;
			extract_one_seq(data->reads, rid, read1);
			reverse_complement(read2, read1, rsize);
			int s;
			char chain;
			num_candidates = 0;
			for (s = 0; s < 2; ++s)
			{
				if (s%2) { chain = 'R'; read = read2; }
				else { chain = 'F'; read = read1; }
				int num_segs = seeding(read, rsize, data->ridx, sbk);
				num_candidates = get_candidates(data->reference, 
												sbk, 
												num_segs, 
												rid + data->reads->start_read_id, 
												rsize, 
												chain, 
												candidates, 
												num_candidates); 
			}

			for (s = 0; s < num_candidates; ++s)
			{
				if (candidates[s].chain == 'F') read = read1;
				else read = read2;
				extract_one_seq(data->reference, candidates[s].readno - data->reference->start_read_id, subject);
				int sstart = candidates[s].loc1;
				int qstart = candidates[s].loc2;
				if (qstart && sstart)
				{
					qstart += kmer_size / 2;
					sstart += kmer_size / 2;
				}
				int ssize = data->reference->offset_list->offset_list[candidates[s].readno - data->reference->start_read_id].size;
				
				int flag = aligner->go(read, qstart, rsize, subject, sstart, ssize, min_align_size);
				
				if (flag)
				{
					fill_m4record(aligner, rid + data->reads->start_read_id, 
								  candidates[s].readno, candidates[s].chain, 
								  rsize, ssize, qstart, sstart, candidates[s].score,
								  m4v + num_m4);
					++num_m4;
				}
			}
			
			append_m4v(m4_list, &m4_list_size, m4v, &num_m4, data->out, &data->result_write_lock);
		}
	}
		
		if (m4_list_size)
		{
			pthread_mutex_lock(&data->result_write_lock);
			print_m4record_list(data->out, m4_list, m4_list_size);
			m4_list_size = 0;
			pthread_mutex_unlock(&data->result_write_lock);
		}
		
		safe_free(read1);
		safe_free(read2);
		safe_free(subject);
		delete sbk;
		delete aligner;
		delete[] m4v;
}

void
candidate_detect(PWThreadData* data, int tid)
{
	char *read, *read1, *read2, *subject;
	safe_malloc(read1, char, MAX_SEQ_SIZE);
	safe_malloc(read2, char, MAX_SEQ_SIZE);
	safe_malloc(subject, char, MAX_SEQ_SIZE);
	SeedingBK* sbk = new SeedingBK(data->reference->curr);
	Candidate candidates[MAXC];
	int num_candidates = 0;
	r_assert(data->ec_results);
	ExtensionCandidate* eclist = data->ec_results[tid];
	int nec = 0;
	ExtensionCandidate ec;

	int rid, Lid, Rid;
	while (1)
	{
		get_next_chunk_reads(data, Lid, Rid);
		if (Lid >= data->reads->num_reads) break;
	for (rid = Lid; rid < Rid; ++rid)
	{
		int rsize = data->reads->offset_list->offset_list[rid].size;
        if (rsize >= MAX_SEQ_SIZE) {
            cout << "rsize = " << rsize << "\t" << MAX_SEQ_SIZE << endl;
            abort();
        }
		extract_one_seq(data->reads, rid, read1);
		reverse_complement(read2, read1, rsize);
		int s;
		int chain;
		num_candidates = 0;
		for (s = 0; s < 2; ++s)
		{
			if (s%2) { chain = REV; read = read2; }
			else { chain = FWD; read = read1; }
			int num_segs = seeding(read, rsize, data->ridx, sbk);
			num_candidates = get_candidates(data->reference, 
											sbk, 
											num_segs, 
											rid + data->reads->start_read_id, 
											rsize, 
											chain, 
											candidates, 
											num_candidates); 
		}
		
		for (s = 0; s < num_candidates; ++s)
		{
			int qstart = candidates[s].loc2;
			int sstart = candidates[s].loc1;
			if (qstart && sstart)
			{
				qstart += kmer_size / 2;
				sstart += kmer_size / 2;
			}
			int qdir = candidates[s].chain;
			int sdir = FWD;
			int qid = rid + data->reads->start_read_id;
			int sid = candidates[s].readno;
			int score = candidates[s].score;
			
			ec.qid = qid;
			ec.qdir = qdir;
			ec.qext = qstart;
			ec.sid = sid;
			ec.sdir = sdir;
			ec.sext = sstart;
			ec.score = score;
			ec.qsize = rsize;
			ec.ssize = data->reference->offset_list->offset_list[sid - data->reference->start_read_id].size;
            if (ec.qdir == REV) ec.qext = ec.qsize - 1 - ec.qext;
            if (ec.sdir == REV) ec.sext = ec.ssize - 1 - ec.sext;
			eclist[nec] = ec;
			++nec;
			if (nec == PWThreadData::kResultListSize)
			{
				pthread_mutex_lock(&data->result_write_lock);
				for (int i = 0; i < nec; ++i) (*data->out) << eclist[i];
				nec = 0;
				pthread_mutex_unlock(&data->result_write_lock);
			}
		}
	}
	}
	
	if (nec)
	{
		pthread_mutex_lock(&data->result_write_lock);
		for (int i = 0; i < nec; ++i) (*data->out) << eclist[i];
		nec = 0;
		pthread_mutex_unlock(&data->result_write_lock);
	}
	
	safe_free(read1);
	safe_free(read2);
	safe_free(subject);
	delete sbk;
}

void*
multi_thread_func(void* p)
{
	PWThreadData* data = (PWThreadData*)p;
	int t = 0;
	pthread_mutex_lock(&data->id_lock);
	t = data->used_thread_id;
	++data->used_thread_id;
	pthread_mutex_unlock(&data->id_lock);
	r_assert(data->options->task == TASK_SEED || data->options->task == TASK_ALN);
	if (data->options->task == TASK_SEED) candidate_detect(data, t);
	else pairwise_mapping(data, t);
	return NULL;
}

void
process_one_volume(options_t* options, const int svid, const int evid, volume_names_t* vn, ostream* out)
{
	MAXC = options->num_candidates;
	output_gapped_start_point = options->output_gapped_start_point;
	min_align_size = options->min_align_size;
	min_kmer_match = options->min_kmer_match;
	
	if (options->tech == TECH_PACBIO) {
		ddfs_cutoff = ddfs_cutoff_pacbio;
		min_kmer_dist = 1800;
	} else if (options->tech == TECH_NANOPORE) {
		ddfs_cutoff = ddfs_cutoff_nanopore;
		min_kmer_dist = 400;
	} else {
		ERROR("TECH must be either %d or %d", TECH_PACBIO, TECH_NANOPORE);
	}
	
	const char* ref_name = get_vol_name(vn, svid);
	volume_t* ref = load_volume(ref_name);
	ref_index* ridx = create_ref_index(ref, kmer_size, options->num_threads);
	pthread_t tids[options->num_threads];
	char volume_process_info[1024];;
	int vid, tid;
	for(vid = svid; vid < evid; ++vid)
	{
		sprintf(volume_process_info, "process volume %d", vid);
		DynamicTimer dtimer(volume_process_info);
		const char* read_name = get_vol_name(vn, vid);
		LOG(stderr, "processing %s\n", read_name);
		volume_t* read = load_volume(read_name);
		PWThreadData* data = new PWThreadData(options, ref, read, ridx, out);
		for (tid = 0; tid < options->num_threads; ++tid)
		{
			int err_code = pthread_create(tids + tid, NULL, multi_thread_func, (void*)data);
			if (err_code)
			{
				LOG(stderr, "Error: return code is %d\n", err_code);
				abort();
			}
		}
		for (tid = 0; tid < options->num_threads; ++tid) pthread_join(tids[tid], NULL);
		read = delete_volume_t(read);
		delete data;
	}
	ref = delete_volume_t(ref);
	ridx = destroy_ref_index(ridx);
}
