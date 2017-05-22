#ifndef MEAP_REF_DEFS_H
#define MEAP_REF_DEFS_H

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <assert.h>

#define RM 100000
#define ZV 1000
#define ZVS 2000
#define ZVL 1000L
#define ZVSL 2000L
#define DN 500
#define SM 20
#define SI 21
#define MK 100
#define MAXSTR 1000000000
#define MAX 10
#define SVM 100000
#define PLL 1000

typedef struct
{
    int aln_str_size,dist,aln_q_s,aln_q_e,aln_t_s,aln_t_e;
    char q_aln_str[2500];
    char t_aln_str[2500];

} alignment;

typedef struct
{
    int x1,y1,x2,y2;
    char *seq1;
    char *seq2;
    float evalue;
} queryresult;

typedef struct
{
    int x,y,k,min,max_index;
} deffpoint;

typedef struct
{
    int slist,llist;
} groupdata;

typedef struct
{
    int pre_k,x1,y1,x2,y2;
} d_path_data;

typedef struct
{
    int d,k,pre_k,x1,y1,x2,y2;
} d_path_data2;

typedef struct
{
    int x,y;
} path_point;

typedef struct
{
    char tempsw1[RM],tempsw2[RM],left_store1[RM],left_store2[RM],right_store1[RM],right_store2[RM],out_store1[RM],out_store2[RM];
} output_store;

typedef struct
{
    int readno,readlen;
    char *seqloc;
} ReadFasta;

struct Back_List
{
    short int score, score2, loczhi[SM],seedno[SM],seednum;
    int index;
};

typedef struct
{
    long loc1,loc2,left1,left2,right1,right2;
    int score, num1,num2;
    char chain;
} candidate_save;

#endif // MEAP_REF_DEFS_H
