#ifndef PW_IMPL_H
#define PW_IMPL_H

#include <iostream>

#include "../common/alignment.h"
#include "../common/packed_db.h"
#include "../common/lookup_table.h"

#define RM 			100000
#define DN 			500
#define BC 			10
#define SM 			40
#define SI 			41
#define CHUNK_SIZE 	500
#define ZV 			2000
#define MUL_ZV(a) 	((a)*ZV)
#define DIV_ZV(a) 	((a)/ZV)
#define MOD_ZV(a) 	((a)%ZV)

struct candidate_save
{
    int loc1,loc2,left1,left2,right1,right2,score,num1,num2,readno,readstart;
    char chain;
};

typedef candidate_save Candidate;

struct Back_List
{
    short score,loczhi[SM],seedno[SM],seednum;
    int index;
};

struct PWThreadData
{
	options_t*				options;
	int 					used_thread_id;
	pthread_mutex_t 		id_lock;
	volume_t* 				reference;
	volume_t* 				reads;
	ref_index* 				ridx;
	std::ostream*			out;	
	M4Record** 				m4_results;
	ExtensionCandidate**	ec_results;
	static const int		kResultListSize = 10000;
	pthread_mutex_t			result_write_lock;
	int						next_processed_id;
	pthread_mutex_t			read_retrieve_lock;
	
	PWThreadData(options_t* opt, volume_t* ref, volume_t* rd, ref_index* idx, std::ostream* o);
	~PWThreadData();
};

struct SeedingBK
{
	int* index_list;
	short* index_score;
	Back_List* database;
	int* kmer_ids;
	
	SeedingBK(const int ref_size);
	~SeedingBK();
};

void
process_one_volume(options_t* options, const int svid, const int evid, volume_names_t* vn, std::ostream* out);

#endif // PW_IMPL_H
