#include "output.h"

#include <assert.h>
#include <string.h>
#include <stdlib.h>

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
		FILE* out)
{
	char sdir = 'F';
	int qb = qstart, qe = qend;
	if (qdir == 'R') 
	{
		sdir = 'R';
		qb = qsize - qend;
		qe = qsize - qstart;
	}
	
	fprintf(out, "%d\t", read_id); // read_id
	fprintf(out, "%s\t", chr_name); // ref name
	fprintf(out, "%c\t", sdir); // ref strand
	fprintf(out, "%d\t", vscore); // voting score
	fprintf(out, "%d\t", qb); // read start
	fprintf(out, "%d\t", qe); // read end
	fprintf(out, "%d\t", qsize); // read size
	fprintf(out, "%ld\t", sstart); // ref start
	fprintf(out, "%ld\t", send); // ref end
	fprintf(out, "%ld\n", ssize); // ref size
	
	fprintf(out, "%s\n", qmap); // mapped read
	fprintf(out, "%s\n", smap); // mapped ref
}

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
				 FILE* out)
{
	int qb = qstart, qe = qend;
	if (qdir == 'R')
	{
		qb = qsize - qend;
		qe = qsize - qstart;
	}
	double ident = 0.0;
	int n = strlen(qmap);
	assert(n > 0);
	int i;
	for (i = 0; i < n; ++i) if (qmap[i] == smap[i]) ident += 1.0;
	ident = ident / n;
	ident *= 100.0;
	
	fprintf(out, "%d\t", read_id); // read id
	fprintf(out, "%s\t", chr_name); // ref name
	fprintf(out, "%.4f\t", ident); // % identity
	fprintf(out, "%d\t", vscore); // voting score
	fprintf(out, "%d\t", (qdir == 'F') ? 0 : 1); // query strand
	fprintf(out, "%d\t", qb); // query start
	fprintf(out, "%d\t", qe); // query end
	fprintf(out, "%d\t", qsize); // query size
	fprintf(out, "0\t"); // ref strand
	fprintf(out, "%ld\t", sstart); // ref start
	fprintf(out, "%ld\t", send); // ref end
	fprintf(out, "%ld", ssize); // ref size
	fprintf(out, "\n");
}

static const char sam_tab = '\t';

void
print_sam_header(FILE* out)
{
	fprintf(out, "@HD\tVN:1.4\tSO:unknown\tGO:query\n");
}

void print_sam_references(fastaindexinfo* fii, const int num_chr, FILE* out)
{
	int i;
	for (i = 0; i < num_chr; ++i)
	{
		fprintf(out, "@SQ\tSN:%s\tLN:%ld\n", fii[i].chrname, fii[i].chrsize);
	}
}

void print_sam_program(int argc, char* argv[], FILE* out)
{
	fprintf(out, "@PG\tID:0\tVN:0.0.1\tCL:");
	int i;
	for (i = 0; i < argc; ++i) fprintf(out, "%s ",argv[i]);
	fprintf(out, "\t");
	fprintf(out, "PN:mecat2ref\n");
}

void
output_cigar(const int qstart,
			 const int qend,
			 const int qsize,
			 const char* qmap,
			 const char* smap,
			 FILE* out)
{
	int qlen = strlen(qmap);
	int slen = strlen(smap);
	assert(qlen == slen);
	if (qstart) fprintf(out, "%dH", qstart);
	int i = 0, j, n = qlen;
	while (i < n)
	{
		if (qmap[i] == '-') // delete from reference
		{
			j = i + 1;
			while (j < n && qmap[j] == '-') ++j;
			fprintf(out, "%dD", j - i);
		}
		else if (smap[i] == '-') // insertion into reference
		{
			j = i + 1;
			while (j < n && smap[j] == '-') ++j;
			fprintf(out, "%dI", j - i);
		}
		else // match or mismatch
		{
			j = i + 1;
			while (j < n && qmap[j] != '-' && smap[j] != '-') ++j;
			fprintf(out, "%dM", j - i);
		}
		i = j;
	}
	if (qend != qsize) fprintf(out, "%dH", qsize - qend);
}

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
		   FILE* out)
{
	int flag = 0;
	if (qdir == 'R') flag = 0x10; // SEQ being reverse complemented
	fprintf(out, "%d\t", read_id); /// 1) qname
	fprintf(out, "%d\t", flag); /// 2) flag
	fprintf(out, "%s\t", chr_name); /// 3) rname
	fprintf(out, "%ld\t", sstart + 1); /// 4) 1-based left most position
	fprintf(out, "255\t"); /// 5) mapq
	output_cigar(qstart, qend, qsize, qmap, smap, out); /// 6) cigar
	fprintf(out, "\t");
	fprintf(out, "*\t"); /// 7) rnext
	fprintf(out, "0\t"); /// 8) pnext
	fprintf(out, "0\t"); /// 9) tlen
	/// 10) seq
	const int n = strlen(qmap);
	int i;
	for (i = 0; i < n; ++i) if (qmap[i] != '-') fprintf(out, "%c", qmap[i]);
	fprintf(out, "\t");
	fprintf(out, "*\n"); /// 11) qual
}
		   
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
				  FILE* out)
{
	if (format == 0)
		print_ref_result(read_id,chr_name,qdir,qstart,qend,qsize,vscore,sstart,send,ssize,qmap,smap,out);
	else if (format == 1)
		print_m4_result(read_id,chr_name,qdir,qstart,qend,qsize,vscore,sstart,send,ssize,qmap,smap,out);
	else if (format == 2)
		output_sam(read_id,chr_name,qdir,qstart,qend,qsize,vscore,sstart,send,ssize,qmap,smap,out);
}

void
set_temp_result(int read_id,
				char read_dir,
				int vscore,
				int qb,
				int qe,
				int qs,
				long sb,
				long se,
				char* qmap,
				char* smap,
				TempResult* result)
{
	result->read_id = read_id;
	result->read_dir = read_dir;
	result->vscore = vscore;
	result->qb = qb;
	result->qe = qe;
	result->qs = qs;
	result->sb = sb;
	result->se = se;
	result->qmap = qmap;
	result->smap = smap;
}

void
output_temp_result(TempResult* result, FILE* out)
{
	fprintf(out, "%d\t%c\t%d\t%d\t%d\t%d\t%ld\t%ld\n%s\n%s\n",
				 result->read_id,
				 result->read_dir,
				 result->vscore,
				 result->qb,
				 result->qe,
				 result->qs,
				 result->sb,
				 result->se,
				 result->qmap,
				 result->smap);
}

int
load_temp_result(TempResult* result, FILE* in)
{
	int r = fscanf(in, "%d\t%c\t%d\t%d\t%d\t%d\t%ld\t%ld\n%s\n%s\n",
				   &result->read_id,
				   &result->read_dir,
				   &result->vscore,
				   &result->qb,
				   &result->qe,
				   &result->qs,
				   &result->sb,
				   &result->se,
				   result->qmap,
				   result->smap);
	return (r == EOF) ? 0 : 1;
}

TempResult*
create_temp_result()
{
	TempResult* result = (TempResult*)calloc(sizeof(TempResult), 1);
	result->qmap = (char*)malloc(100000);
	result->smap = (char*)malloc(100000);
	return result;
}

TempResult*
destroy_temp_result(TempResult* result)
{
	free(result->qmap);
	free(result->smap);
	free(result);
	return NULL;
}

void
copy_temp_result(TempResult* src, TempResult* dst)
{
	dst->read_id = src->read_id;
	dst->read_dir = src->read_dir;
	dst->vscore = src->vscore;
	dst->qb = src->qb;
	dst->qe = src->qe;
	dst->qs = src->qs;
	dst->sb = src->sb;
	dst->se = src->se;
	strcpy(dst->qmap, src->qmap);
	strcpy(dst->smap, src->smap);
}
