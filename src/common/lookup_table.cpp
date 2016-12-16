#include "lookup_table.h"
#include "packed_db.h"

#include <cstdio>

ref_index*
destroy_ref_index(ref_index* ridx)
{
	safe_free(ridx->kmer_counts);
	safe_free(ridx->kmer_starts);
	safe_free(ridx->kmer_offsets);
	safe_free(ridx);
	return NULL;
}

typedef struct
{
	ref_index* ridx;
	uint32_t min_key;
	uint32_t max_key;
	volume_t* v;
	int kmer_size;
} ref_index_thread_info;

void*
fill_ref_index_offsets_func(void* arg)
{
	ref_index_thread_info* riti = (ref_index_thread_info*)(arg);
	volume_t* v = riti->v;
	int num_reads = v->num_reads;
	ref_index* index = riti->ridx;
	int kmer_size = riti->kmer_size;
	uint32_t index_count = 1 << (kmer_size * 2);
	uint32_t leftnum = 34 - 2 * kmer_size;
	int i, j;
	for (i = 0; i < num_reads; ++i)
	{
		int read_start = v->offset_list->offset_list[i].offset;
		int read_size = v->offset_list->offset_list[i].size;
		uint32_t eit = 0;
		for (j = 0; j < read_size; ++j)
		{
			int k = read_start + j;
			uint8_t c = PackedDB::get_char(v->data, k);
			eit = (eit << 2) | c;
			assert(eit < index_count);
			if (j >= kmer_size - 1)
			{
				if (index->kmer_starts[eit] && eit >= riti->min_key && eit <= riti->max_key)
				{
					index->kmer_starts[eit][index->kmer_counts[eit]] = k + 1 - kmer_size;
					++index->kmer_counts[eit];
				}
				eit <<= leftnum;
				eit >>= leftnum;
			}
		}
	}
	
	return NULL;
}

ref_index*
create_ref_index(volume_t* v, int kmer_size, int num_threads)
{
	DynamicTimer dtimer(__func__);
	uint32_t index_count = 1 << (kmer_size * 2);
	uint32_t leftnum = 34 - 2 * kmer_size;
	ref_index* index = (ref_index*)malloc(sizeof(ref_index));
	safe_calloc(index->kmer_counts, int, index_count);
	int num_reads = v->num_reads;
	for (uint32_t i = 0; i != index_count; ++i) assert(index->kmer_counts[i] == 0);
	for (int i = 0; i != num_reads; ++i)
	{
		int read_start = v->offset_list->offset_list[i].offset;
		int read_size = v->offset_list->offset_list[i].size;
		uint32_t eit = 0;
		for (int j = 0; j < read_size; ++j)
		{
			int k = read_start + j;
			uint8_t c = PackedDB::get_char(v->data, k);
			assert(c>= 0 && c < 4);
			eit = (eit << 2) | c;
			if (j >= kmer_size - 1)
			{
				assert(eit < index_count);
				++index->kmer_counts[eit];
				eit = eit << leftnum;
				eit = eit >> leftnum;
			}
		}
	}
	
	int num_kmers = 0;
	for (uint32_t i = 0; i != index_count; ++i) 
	{
		if (index->kmer_counts[i] > 128) index->kmer_counts[i] = 0;
		num_kmers += index->kmer_counts[i];
	}
	printf("number of kmers: %d\n", num_kmers);
	safe_malloc(index->kmer_offsets, int, num_kmers);
	safe_malloc(index->kmer_starts, int*, index_count);
	
	if (v->curr < 10 * 1000000) num_threads = 1;
	int kmers_per_thread = (num_kmers + num_threads - 1) / num_threads;
	fprintf(stderr, "%d threads are used for filling offset lists.\n", num_threads);
	uint32_t hash_boundaries[2 * num_threads];
	uint32_t L = 0;
	num_kmers = 0;
	int kmer_cnt = 0;
	int tid = 0;
	for (uint32_t i = 0; i != index_count; ++i)
	{
		if (index->kmer_counts[i])
		{
			index->kmer_starts[i] = index->kmer_offsets + num_kmers;
			num_kmers += index->kmer_counts[i];
			kmer_cnt += index->kmer_counts[i];
			index->kmer_counts[i] = 0;
			
			if (kmer_cnt >= kmers_per_thread)
			{
				printf("thread %d: %d\t%d\n", tid, L, i);
				hash_boundaries[2 * tid] = L;
				hash_boundaries[2 * tid + 1] = i;
				++tid;
				L = i + 1;
				kmer_cnt = 0;
			}
		}
		else
		{
			index->kmer_starts[i] = NULL;
		}
	}
	if (kmer_cnt)
	{
		printf("thread %d: %d\t%d\n", tid, L, index_count - 1);
		hash_boundaries[2 * tid] = L;
		hash_boundaries[2 * tid + 1] = index_count - 1;
	}
	
	ref_index_thread_info ritis[num_threads];
	for (int i = 0; i != num_threads; ++i)
	{
		ritis[i].ridx = index;
		ritis[i].min_key = hash_boundaries[2 * i];
		ritis[i].max_key = hash_boundaries[2 * i + 1];
		ritis[i].v = v;
		ritis[i].kmer_size = kmer_size;
	}
	
	pthread_t tids[num_threads];
	for (int j = 0; j < num_threads; ++j)
		pthread_create(tids + j, NULL, fill_ref_index_offsets_func, (void*)(ritis + j));
	for (int j = 0; j < num_threads; ++j)
		pthread_join(tids[j], NULL);
	
	return index;
}
