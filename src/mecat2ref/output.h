#ifndef _OUTPUT_H
#define _OUTPUT_H

#include <stdio.h>

#define FMT_REF 0
#define FMT_M4 1
#define FMT_SAM 2

typedef struct 
{
	long chrstart;
	long chrsize;
	char chrname[64];
}fastaindexinfo;

void
print_m4_result(const int read_id,
				const char* chr_name,
				const char qdir,
				const int qstart,
				const int qend,
				const int qsize,
				const int vscore,
				const long sstart,
				const long send,
				const long ssize,
				const char* qmap,
				const char* smap,
				FILE* out);

void
print_ref_result(const int read_id,
				 const char* chr_name,
				 const char qdir,
				 const int qstart,
				 const int qend,
				 const int qsize,
				 const int vscore,
				 const long sstart,
				 const long send,
				 const long ssize,
				 const char* qmap,
				 const char* smap,
				 FILE* out);

void print_sam_header(FILE* out);

void print_sam_references(fastaindexinfo* fii, const int num_chr, FILE* out);

void print_sam_program(int argc, char* argv[], FILE* out);

void
output_cigar(const int qstart,
			const int qend,
			const int qsize,
			const char* qmap,
			const char* smap,
			FILE* out);

void
output_sam(const int read_id,
		   const char* chr_name,
		   const char qdir,
		   const int qstart,
		   const int qend,
		   const int qsize,
		   const int vscore,
		   const long sstart,
		   const long send,
		   const long ssize,
		   const char* qmap,
		   const char* smap,
		   FILE* out);

void
output_one_result(const int read_id,
				  const char* chr_name,
				  const char qdir,
				  const int qstart,
				  const int qend,
				  const int qsize,
				  const int vscore,
				  const long sstart,
				  const long send,
				  const long ssize,
				  const char* qmap,
				  const char* smap,
				  const int format,
				  FILE* out);

typedef struct 
{
	int read_id;
	char read_dir;
	int vscore;
	int qb;
	int qe;
	int qs;
	long sb;
	long se;
	int aln_size;
	char* qmap;
	char* smap;
} TempResult;

TempResult*
create_temp_result();

TempResult*
destroy_temp_result(TempResult* result);

void
copy_temp_result(TempResult* src, TempResult* dst);

void
set_temp_result(int read_id,
				char read_dir,
				int vscore,
				int qb,
				int qe,
				int qs,
				long s,
				long se,
				char* qmap,
				char* smap,
				TempResult* result);

void
output_temp_result(TempResult* result, FILE* out);

int
load_temp_result(TempResult* result, FILE* in);

#endif // _OUTPUT_H
