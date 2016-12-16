#include "mecat2ref_defs.h"
#include "output.h"

#include <algorithm>
using namespace std;

static int MAXC = 0;

typedef struct
{
    int loc1,loc2,left1,left2,right1,right2,score,num1,num2,readno,readstart;
    char chain;
} candidate_save;

static pthread_t *thread;
static int threadnum=2;
static FILE **outfile;
static pthread_mutex_t mutilock; 
static int runnumber=0,runthreadnum=0, readcount,terminalnum;
static int *countin,**databaseindex,*allloc,sumcount;
static int seed_len,seqcount;
static char *REFSEQ;
static char *savework,workpath[300],fastqfile[300];
static ReadFasta *readinfo;

static int compare_d_path(const void * a, const void * b)
{
    const d_path_data2 * arg1 = (d_path_data2 *)a;
    const d_path_data2 * arg2 = (d_path_data2 *)b;
    if (arg1->d - arg2->d == 0)
    {
        return  arg1->k - arg2->k;
    }
    else
    {
        return arg1->d - arg2->d;
    }
}

static d_path_data2 * get_dpath_idx( int d, int k, unsigned long max_idx, d_path_data2 * base)
{
    d_path_data2 d_tmp;
    d_path_data2 *rtn;
    d_tmp.d = d;
    d_tmp.k = k;
    rtn = (d_path_data2 *)  bsearch( &d_tmp, base, max_idx, sizeof(d_path_data2), compare_d_path);
    //printf("dp %ld %ld %ld %ld %ld %ld %ld\n", (rtn)->d, (rtn)->k, (rtn)->x1, (rtn)->y1, (rtn)->x2, (rtn)->y2, (rtn)->pre_k);
    return rtn;

}

struct SCompareDPathData2
{
    bool operator()(const d_path_data2& a, const d_path_data2& b)
    { return (a.d == b.d) ? (a.k < b.k) : (a.d < b.d); }
};

static int align(char * query_seq, char * target_seq,int band_tolerance,int get_aln_str,alignment * align_rtn,int * V,int * U,d_path_data2 *d_path,path_point * aln_path,deffpoint *deffvalue)
{
    int k_offset,d,k, k2,best_m,min_k, new_min_k,max_k, new_max_k,pre_k,x, y,dl;
    int ck,cd,cx, cy, nx, ny,max_d,band_size,q_len,t_len;
    unsigned long d_path_idx = 0,max_idx = 0;
    int aln_path_idx,aln_pos,i,aligned=0;
    d_path_data2 * d_path_aux;
    q_len=strlen(query_seq);
    t_len=strlen(target_seq);
    max_d = (int) (0.3*(q_len + t_len));
    band_size = band_tolerance * 2;
    k_offset = max_d;
    align_rtn->aln_str_size = 0;
    align_rtn->aln_q_s = 0;
    align_rtn->aln_q_e = 0;
    align_rtn->aln_t_s = 0;
    align_rtn->aln_t_e = 0;
    best_m = -1;
    min_k = 0;
    max_k = 0;
    d_path_idx = 0;
    max_idx = 0;
    for (d = 0; d < max_d; d ++ )
    {
        if(max_k - min_k > band_size)break;
        for (k = min_k; k <= max_k;  k += 2)
        {
            if ( k == min_k || (k != max_k && V[ k - 1 + k_offset ] < V[ k + 1 + k_offset]) )
            {
                pre_k = k + 1;
                x = V[ k + 1 + k_offset];
            }
            else
            {
                pre_k = k - 1;
                x = V[ k - 1 + k_offset] + 1;
            }
            y = x - k;
            //search the same string and save the pathway
            d_path[d_path_idx].d = d;
            d_path[d_path_idx].k = k;
            d_path[d_path_idx].x1 = x;
            d_path[d_path_idx].y1 = y;
            while ( x < q_len && y < t_len && query_seq[x] == target_seq[y] )
            {
                x++;
                y++;
            }
            d_path[d_path_idx].x2 = x;
            d_path[d_path_idx].y2 = y;
            d_path[d_path_idx].pre_k = pre_k;
            d_path_idx ++;
            V[ k + k_offset ] = x;
            U[ k + k_offset ] = x + y;
            if( x + y > best_m)
            {
                best_m = x + y;
                if(d%20==0)
                {
                    dl=d/20;
                    if(dl==0)
                    {
                        deffvalue[dl].x=0;
                        deffvalue[dl].y=0;
                        deffvalue[dl].max_index=0;
                        deffvalue[dl].min=0;
                        deffvalue[dl].k=0;
                    }
                    else
                    {
                        deffvalue[dl].x=x;
                        deffvalue[dl].y=best_m;
                        deffvalue[dl].max_index=d_path_idx;
                        if(x<y)deffvalue[dl].min=x;
                        else deffvalue[dl].min=y;
                        deffvalue[dl].k=k;
                    }
                }
            }
            if ( x >= q_len || y >= t_len)
            {
                aligned = 1;
                max_idx = d_path_idx;
                break;
            }
        }

        // For banding
        new_min_k = max_k;
        new_max_k = min_k;
        for (k2 = min_k; k2 <= max_k;  k2 += 2)
        {
            if (U[ k2 + k_offset] >= best_m - band_tolerance )
            {
                if( k2 < new_min_k )new_min_k = k2;
                if( k2 > new_max_k )new_max_k = k2;
            }
        }
        max_k = new_max_k + 1;
        min_k = new_min_k - 1;
        /* if(d%20==0)
         {
             dl=d/20;
             if(dl>1&&(20.0/(deffvalue[dl].x-deffvalue[dl-1].x)>0.5||20.0/((deffvalue[dl].y-deffvalue[dl].x)-(deffvalue[dl-1].y-deffvalue[dl-1].x))>0.5))
             {
                 aligned=2;
                 dl=dl-1;
             }
         }*/


        if (aligned == 1)
        {
            align_rtn->aln_q_e = x;
            align_rtn->aln_t_e = y;
            align_rtn->dist = d;
            align_rtn->aln_str_size = (x + y + d) / 2;
            align_rtn->aln_q_s = 0;
            align_rtn->aln_t_s = 0;
            std::sort(d_path, d_path + max_idx, SCompareDPathData2());
            if (get_aln_str > 0)
            {
                cd = d;
                ck = k;
                aln_path_idx = 0;
                while (cd >= 0 && aln_path_idx < q_len + t_len + 1)
                {
                    d_path_aux = (d_path_data2 *) get_dpath_idx( cd, ck, max_idx, d_path);
                    aln_path[aln_path_idx].x = d_path_aux -> x2;
                    aln_path[aln_path_idx].y = d_path_aux -> y2;
                    aln_path_idx ++;
                    aln_path[aln_path_idx].x = d_path_aux -> x1;
                    aln_path[aln_path_idx].y = d_path_aux -> y1;
                    aln_path_idx ++;
                    ck = d_path_aux -> pre_k;
                    cd -= 1;
                }
                aln_path_idx --;
                cx = aln_path[aln_path_idx].x;
                cy = aln_path[aln_path_idx].y;
                align_rtn->aln_q_s = cx;
                align_rtn->aln_t_s = cy;
                aln_pos = 0;
                while ( aln_path_idx > 0 )
                {
                    aln_path_idx --;
                    nx = aln_path[aln_path_idx].x;
                    ny = aln_path[aln_path_idx].y;
                    if (cx == nx && cy == ny)
                    {
                        continue;
                    }
                    if (nx == cx && ny != cy)  //advance in y
                    {
                        for (i = 0; i <  ny - cy; i++)align_rtn->q_aln_str[aln_pos + i] = '-';
                        for (i = 0; i <  ny - cy; i++)align_rtn->t_aln_str[aln_pos + i] = target_seq[cy + i];
                        aln_pos += ny - cy;
                    }
                    else if (nx != cx && ny == cy)  //advance in x
                    {
                        for (i = 0; i <  nx - cx; i++)align_rtn->q_aln_str[aln_pos + i] = query_seq[cx + i];
                        for (i = 0; i <  nx - cx; i++)align_rtn->t_aln_str[aln_pos + i] = '-';
                        aln_pos += nx - cx;
                    }
                    else
                    {
                        for (i = 0; i <  nx - cx; i++)align_rtn->q_aln_str[aln_pos + i] = query_seq[cx + i];
                        for (i = 0; i <  ny - cy; i++)align_rtn->t_aln_str[aln_pos + i] = target_seq[cy + i];
                        aln_pos += ny - cy;
                    }
                    cx = nx;
                    cy = ny;
                }
                align_rtn->aln_str_size = aln_pos;
            }
            break;
        }
        else if(aligned == 2)
        {
            //get start value
            align_rtn->aln_q_e = deffvalue[dl].x;
            align_rtn->aln_t_e = deffvalue[dl].y-deffvalue[dl].x;
            d=dl*20;
            max_idx=deffvalue[dl].max_index;
            k=deffvalue[dl].k;
            q_len=deffvalue[dl].x;
            t_len=deffvalue[dl].y-deffvalue[dl].x;

            align_rtn->dist = d;
            align_rtn->aln_str_size = (x + y + d) / 2;
            align_rtn->aln_q_s = 0;
            align_rtn->aln_t_s = 0;
            qsort(d_path, max_idx, sizeof(d_path_data2), compare_d_path);
            if (get_aln_str > 0)
            {
                cd = d;
                ck = k;
                aln_path_idx = 0;
                while (cd >= 0 && aln_path_idx < q_len + t_len + 1)
                {
                    d_path_aux = (d_path_data2 *) get_dpath_idx( cd, ck, max_idx, d_path);
                    aln_path[aln_path_idx].x = d_path_aux -> x2;
                    aln_path[aln_path_idx].y = d_path_aux -> y2;
                    aln_path_idx ++;
                    aln_path[aln_path_idx].x = d_path_aux -> x1;
                    aln_path[aln_path_idx].y = d_path_aux -> y1;
                    aln_path_idx ++;
                    ck = d_path_aux -> pre_k;
                    cd -= 1;
                }
                aln_path_idx --;
                cx = aln_path[aln_path_idx].x;
                cy = aln_path[aln_path_idx].y;
                align_rtn->aln_q_s = cx;
                align_rtn->aln_t_s = cy;
                aln_pos = 0;
                while ( aln_path_idx > 0 )
                {
                    aln_path_idx --;
                    nx = aln_path[aln_path_idx].x;
                    ny = aln_path[aln_path_idx].y;
                    if (cx == nx && cy == ny)
                    {
                        continue;
                    }
                    if (nx == cx && ny != cy)  //advance in y
                    {
                        for (i = 0; i <  ny - cy; i++)align_rtn->q_aln_str[aln_pos + i] = '-';
                        for (i = 0; i <  ny - cy; i++)align_rtn->t_aln_str[aln_pos + i] = target_seq[cy + i];
                        aln_pos += ny - cy;
                    }
                    else if (nx != cx && ny == cy)  //advance in x
                    {
                        for (i = 0; i <  nx - cx; i++)align_rtn->q_aln_str[aln_pos + i] = query_seq[cx + i];
                        for (i = 0; i <  nx - cx; i++)align_rtn->t_aln_str[aln_pos + i] = '-';
                        aln_pos += nx - cx;
                    }
                    else
                    {
                        for (i = 0; i <  nx - cx; i++)align_rtn->q_aln_str[aln_pos + i] = query_seq[cx + i];
                        for (i = 0; i <  ny - cy; i++)align_rtn->t_aln_str[aln_pos + i] = target_seq[cy + i];
                        aln_pos += ny - cy;
                    }
                    cx = nx;
                    cy = ny;
                }
                align_rtn->aln_str_size = aln_pos;
            }
            break;
        }
    }
    if(align_rtn->aln_q_e==q_len||align_rtn->aln_t_e==t_len)return(aligned);
    else return(0);
}

static int filesize(FILE *stream)
{
    int curpos, length;
    curpos = ftell(stream);
    fseek(stream, 0L, SEEK_END);
    length = ftell(stream);
    fseek(stream, curpos, SEEK_SET);
    return length;
}

static unsigned short atcttrans(char c)
{
    if(c=='A')return 0;
    else if(c=='T')return 1;
    else if(c=='C')return 2;
    else if(c=='G')return 3;
    else return 4;
}

static int sumvalue_x(int *intarry,int count)
{
    int i,sumval=0;
    for(i=0; i<count; i++)
    {
        if(intarry[i]>0&&intarry[i]<129)sumval=sumval+intarry[i];
        else if(intarry[i]>128)intarry[i]=0;
    }
    return(sumval);
}

static int transnum_buchang(char *seqm,int *value,int *endn,int len_str,int readnum,int BC)
{
    int eit=0,temp;
    int i,j,start,num;
    num=(len_str-readnum)/BC+1;
    *endn=(len_str-readnum)%BC;
    //if((len_str%readnum)>0)num=num+1;
    for(i=0; i<num; i++)
    {
        eit=0;
        start=i*BC;
        //if(i==num-1){eit=0;start=len_str-readnum;}
        for(j=0; j<readnum; j++)
        {
            temp=atcttrans(seqm[start+j]);
            if(temp==4)
            {
                eit=-1;
                break;
            }
            eit=eit<<2;
            eit=eit+temp;
        }
        value[i]=eit;
    }
    return(num);
}

static void string_check(char *seq1,char *seq2,char *str1,char *str2)
{
    //char seq1[500]="GACCGCCGGACAGCCCACAAACACAACAGCATTTGGCGTATTTCCCGTCAAAGGACTGCGAGTGGGACCGCGCACCGATTTATAGAGTAACGGTGGGACTTACCCCCGACGACTAGAGG";
    //char seq2[500]="GACCGCCGGACAGGCCACAAACACAAATCCGCGAGGCGTATTCCTGTCAAAGGGACTACGGCCCAGTGGGACGCCGCACGACTATATAGTAGTAAGGTGTGCTTTACCCGACGCCCTAGAGG";
    //char str1[700]="GACCGCCGGACAG-CCCACAAACACAACA---GC-ATTTGGCGTATTTCCC-GTCAAAGG-ACTG-CG----AGTGGGACCGC-GCACCGA-T-TTATAG-AGTAACGGTG-GGACTT-ACCCCCGACGAC--TAGAGG";
    //char str2[700]="GACCGCCGGACAGGCC-ACAAACACAA-ATCCGCGA---GGCGTATT-CC-TGTCAAAGGGACT-ACGGCCCAGTGGGAC-GCCGCAC-GACTAT-ATAGTAGTAA-GGTGTG--CTTTACCC--GACG-CCCTAGAGG";
    int len1,len2,slen1,loc1,loc2,k,s,j;
    len1=strlen(seq1)-1;
    len2=strlen(seq2)-1;
    for(slen1=strlen(str1)-1,loc1=0,loc2=0; slen1>-1; slen1--)
    {
        if(str1[slen1]!='-')loc1++;
        else if(loc1<=len1&&loc2<=len2&&seq1[len1-loc1]==seq2[len2-loc2])
        {
            k=1;
            while(loc1+k<=len1&&loc2+k<=len2&&seq1[len1-loc1-k]==seq2[len2-loc2-k])k++;
            s=0;
            j=slen1;
            while(s<k)
            {
                if(str1[j]!='-')
                {
                    str1[j]='-';
                    s++;
                }
                j--;
            }
            s=0;
            j=slen1;
            while(s<k)
            {
                if(str2[j]!='-')
                {
                    str2[j]='-';
                    s++;
                }
                j--;
            }
            for(s=0,j=slen1; s<k; j--)
            {
                str1[j]=seq1[len1-loc1-s];
                str2[j]=seq2[len2-loc2-s];
                s++;
            }
            if(str1[slen1]!='-')loc1++;
        }
        if(str2[slen1]!='-')loc2++;
        else if(loc1-1<=len1&&loc2<=len2&&seq1[len1-loc1+1]==seq2[len2-loc2])
        {
            k=1;
            while(loc1+k-1<=len1&&loc2+k<=len2&&seq1[len1-loc1+1-k]==seq2[len2-loc2-k])k++;
            s=0;
            j=slen1;
            while(s<k)
            {
                if(str1[j]!='-')
                {
                    str1[j]='-';
                    s++;
                }
                j--;
            }
            s=0;
            j=slen1;
            while(s<k)
            {
                if(str2[j]!='-')
                {
                    str2[j]='-';
                    s++;
                }
                j--;
            }
            for(s=0,j=slen1; s<k; j--)
            {
                str1[j]=seq1[len1-loc1+1-s];
                str2[j]=seq2[len2-loc2-s];
                s++;
            }
            if(str2[slen1]!='-')loc2++;
        }


    }
}

static void insert_loc(struct Back_List *spr,int loc,int seedn,float len)
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
    for(i=0; i<SM; i++)for(j=i+1; j<11; j++)if(list_seed[j]-list_seed[i]>0&&list_loc[j]-list_loc[i]>0&&fabs((list_loc[j]-list_loc[i])/((list_seed[j]-list_seed[i])*len)-1.0)<0.3)
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
static int find_location(int *t_loc,int *t_seedn,int *t_score,int *loc,int k,int *rep_loc,float len,int read_len1)
{
    int i,j,maxval=0,maxi,rep=0,lasti = 0;
    for(i=0; i<k; i++)t_score[i]=0;
    for(i=0; i<k-1; i++)for(j=i+1; j<k; j++)if(t_seedn[j]-t_seedn[i]>0&&t_loc[j]-t_loc[i]>0&&t_loc[j]-t_loc[i]<read_len1&&fabs((t_loc[j]-t_loc[i])/((t_seedn[j]-t_seedn[i])*len)-1)<0.25)
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
        for(j=0; j<maxi; j++)if(t_seedn[maxi]-t_seedn[j]>0&&t_loc[maxi]-t_loc[j]>0&&t_loc[maxi]-t_loc[j]<read_len1&&fabs((t_loc[maxi]-t_loc[j])/((t_seedn[maxi]-t_seedn[j])*len)-1)<0.20)
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
        for(j=maxi+1; j<k; j++)if(t_seedn[j]-t_seedn[maxi]>0&&t_loc[j]-t_loc[maxi]>0&&t_loc[j]-t_loc[maxi]<=read_len1&&fabs((t_loc[j]-t_loc[maxi])/((t_seedn[j]-t_seedn[maxi])*len)-1)<0.20)
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

static void creat_ref_index(char *fastafile)
{
    unsigned int eit,temp;
    int i, start, indexcount=0,leftnum=0;
    long length,count, rsize = 0;
    FILE *fasta,*fastaindex;
    char *seq,ch,nameall[200];
    if(seed_len==14)indexcount=268435456;
    else if(seed_len==13)indexcount=67108864;
    else if(seed_len==12)indexcount=16777216;
    else if(seed_len==11)indexcount=4194304;
    else if(seed_len==10)indexcount=1048576;
    else if(seed_len==9)indexcount=262144;
    else if(seed_len==8)indexcount=65536;
    else if(seed_len==7)indexcount=16384;
    else if(seed_len==6)indexcount=4096;
    leftnum=34-2*seed_len;
    //read reference seq
    fasta=fopen(fastafile, "r");
    length=filesize(fasta);
    sprintf(nameall,"%s/chrindex.txt",workpath);
    fastaindex=fopen(nameall,"w");

    REFSEQ=(char *)malloc((length+1000)*sizeof(char));
    seq=REFSEQ;
    for (ch=getc(fasta),count=0; ch!=EOF; ch=getc(fasta))
    {
        if(ch=='>')
        {
            assert(fscanf(fasta,"%[^\n]s",nameall) == 1);
			if (rsize) fprintf(fastaindex, "%ld\n", rsize);
			rsize = 0;
			for(i=0;i<strlen(nameall);i++)if(nameall[i]==' '||nameall[i]=='\t')break;
			nameall[i]='\0';
            fprintf(fastaindex,"%ld\t%s\t",count,nameall);
        }
        else if(ch!='\n'&&ch!='\r')
        {
            if(ch>'Z')ch=toupper(ch);
            seq[count]=ch;
            count=count+1;
			++rsize;
        }
    }
    fclose(fasta);
	fprintf(fastaindex, "%ld\n", rsize);
    fprintf(fastaindex,"%ld\t%s\n",count,"FileEnd");
    seq[count]='\0';
    fclose(fastaindex);
    seqcount=count;

    countin=(int *)malloc((indexcount)*sizeof(int));
    for(i=0; i<indexcount; i++)countin[i]=0;

    eit=0;
    start=0;
    for(i=0; i<seqcount; i++)
    {
        if(seq[i]=='N'||(temp=atcttrans(seq[i]))==4)
        {
            eit=0;
            start=0;
            continue;
        }
        temp=atcttrans(seq[i]);
        if(start<seed_len-1)
        {
            eit=eit<<2;
            eit=eit+temp;
            start=start+1;
        }
        else if(start>=seed_len-1)
        {
            eit=eit<<2;
            eit=eit+temp;
            start=start+1;
            countin[eit]=countin[eit]+1;
            eit=eit<<leftnum;
            eit=eit>>leftnum;
        }
    }

    sumcount=sumvalue_x(countin,indexcount);
    allloc=(int *)malloc(sumcount*sizeof(int));
    databaseindex=(int **)malloc((indexcount)*sizeof(int *));
    sumcount=0;
    for(i=0; i<indexcount; i++)
    {
        if(countin[i]>0)
        {
            databaseindex[i]=allloc+sumcount;
            sumcount=sumcount+countin[i];
            countin[i]=0;
        }
        else databaseindex[i]=NULL;
    }

    eit=0;
    start=0;
    for(i=0; i<seqcount; i++)
    {
        if(seq[i]=='N'||(temp=atcttrans(seq[i]))==4)
        {
            eit=0;
            start=0;
            continue;
        }
        temp=atcttrans(seq[i]);
        if(start<seed_len-1)
        {
            eit=eit<<2;
            eit=eit+temp;
            start=start+1;
        }
        else if(start>=seed_len-1)
        {
            eit=eit<<2;
            eit=eit+temp;
            start=start+1;

            if(databaseindex[eit]!=NULL)
            {
                countin[eit]=countin[eit]+1;
                databaseindex[eit][countin[eit]-1]=i+2-seed_len;
            }
            eit=eit<<leftnum;
            eit=eit>>leftnum;
        }
    }
}

static void reference_mapping(int threadint)
{
    int eit;
    int cleave_num,read_len,s_k,loc;
    int mvalue[10000],*leadarray,loc_flag,flag_end,u_k;
    int count1=0,i,j,k,templong,read_name;
    struct Back_List *database,*temp_spr,*temp_spr1;
    int location_loc[4],repeat_loc,*index_list,*index_spr;
    short int *index_score,*index_ss;
    int temp_list[200],temp_seedn[200],temp_score[200],start_loc;
    int sci=0,localnum,read_i,read_end,fileid;
    int endnum,ii;
    char *seq,*onedata,onedata1[RM],onedata2[RM],seq1[2500],seq2[2500],*seq_pr1,*seq_pr2,FR;
    int left_loc,left_loc1,right_loc=0,right_loc1,cc1,canidatenum,loc_seed,loc_list;
    int num1,num2,BC;
    int low,high,mid,seedcount;
    alignment strvalue;
    int longstr1[2000],longstr2[2000],left_length1,right_length1,left_length2,right_length2,align_flag;
    d_path_data2 *d_path;
    candidate_save canidate_loc[MAXC],canidate_temp;
    path_point aln_path[5000];
    output_store *resultstore;
    deffpoint deffvalue[100];
    seq=REFSEQ;
    j=seqcount/ZV+5;
    index_list=(int *)malloc(j*sizeof(int));
    index_score=(short int *)malloc(j*sizeof(short int));
    database=(struct Back_List *)malloc(j*sizeof(struct Back_List));
    for(i=0,temp_spr=database; i<j; temp_spr++,i++)
    {
        temp_spr->score=0;
        temp_spr->index=-1;
    }
    d_path=(d_path_data2*)malloc(600*700*2*sizeof(d_path_data2));
    resultstore=(output_store *)malloc(1*sizeof(output_store));
	TempResult result;

    fileid=1;
    while(fileid)
    {
        pthread_mutex_lock(&mutilock);
        localnum=runnumber;
        runnumber++;
        pthread_mutex_unlock(&mutilock);
        if(localnum>=terminalnum)
        {
            fileid=0;
            break;
        }
        if(localnum==terminalnum-1)read_end=readcount;
        else read_end=(localnum+1)*PLL;
        for(read_i=localnum*PLL; read_i<read_end; read_i++)
        {
            read_name=readinfo[read_i].readno;
            read_len=readinfo[read_i].readlen;
            strcpy(onedata1,readinfo[read_i].seqloc);
            canidatenum=0;
            loc_flag=0;
            for(ii=1; ii<=2; ii++)
            {
                read_len=strlen(onedata1);
                BC=5+(read_len/1000);
                if(BC>20)BC=20;
                if(ii==1)onedata=onedata1;
                else if(ii==2)
                {
                    strcpy(onedata2,onedata1);
                    onedata=onedata2;
                    for(j=read_len-1,i=0; j>i; j--,i++)
                    {
                        FR=onedata[i];
                        onedata[i]=onedata[j];
                        onedata[j]=FR;
                    }
                    for(i=0; i<read_len; i++)
                    {
                        FR=onedata[i];
                        switch(FR)
                        {
                        case 'A':
                        {
                            onedata[i]='T';
                            break;
                        }
                        case 'T':
                        {
                            onedata[i]='A';
                            break;
                        }
                        case 'C':
                        {
                            onedata[i]='G';
                            break;
                        }
                        case 'G':
                        {
                            onedata[i]='C';
                            break;
                        }
                        }
                    }
                }
                endnum=0;
                read_len=strlen(onedata);
                cleave_num=transnum_buchang(onedata,mvalue,&endnum,read_len,seed_len,BC);
                j=0;
                index_spr=index_list;
                index_ss=index_score;
                endnum=0;
                for(k=0; k<cleave_num; k++)if(mvalue[k]>=0)
                    {
                        count1=countin[mvalue[k]];
                        //if(count1>20)continue;
                        leadarray=databaseindex[mvalue[k]];
                        for(i=0; i<count1; i++,leadarray++)
                        {
                            templong=(*leadarray)/ZV;
                            u_k=(*leadarray)%ZV;
                            if(templong>=0)
                            {
                                temp_spr=database+templong;
                                if(temp_spr->score==0||temp_spr->seednum<k+1)
                                {
                                    loc=++(temp_spr->score);
                                    if(loc<=SM)
                                    {
                                        temp_spr->loczhi[loc-1]=u_k;
                                        temp_spr->seedno[loc-1]=k+1;
                                    }
                                    else insert_loc(temp_spr,u_k,k+1,BC);
                                    if(templong>0)s_k=temp_spr->score+(temp_spr-1)->score;
                                    else s_k=temp_spr->score;
                                    if(endnum<s_k)endnum=s_k;
                                    if(temp_spr->index==-1)
                                    {
                                        *(index_spr++)=templong;
                                        *(index_ss++)=s_k;
                                        temp_spr->index=j;
                                        j++;
                                    }
                                    else index_score[temp_spr->index]=s_k;
                                }
                                temp_spr->seednum=k+1;
                            }
                        }
                    }
                cc1=j;
                for(i=0,index_spr=index_list,index_ss=index_score; i<cc1; i++,index_spr++,index_ss++)if(*index_ss>6)
                    {
                        temp_spr=database+*index_spr;
                        if(temp_spr->score==0)continue;
                        s_k=temp_spr->score;
						if(*index_spr>0)loc=(temp_spr-1)->score;
						else loc=0;
                        start_loc=(*index_spr)*ZV;
                        if(*index_spr>0)
                        {
                            loc=(temp_spr-1)->score;
                            if(loc>0)start_loc=(*index_spr-1)*ZV;
                        }
                        else loc=0;
                        if(loc==0)for(j=0,u_k=0; j<s_k&&j<SM; j++)
                            {
                                temp_list[u_k]=temp_spr->loczhi[j];
                                temp_seedn[u_k]=temp_spr->seedno[j];
                                u_k++;
                            }
                        else
                        {
                            k=loc;
                            u_k=0;
                            temp_spr1=temp_spr-1;
                            for(j=0; j<k&&j<SM; j++)
                            {
                                temp_list[u_k]=temp_spr1->loczhi[j];
                                temp_seedn[u_k]=temp_spr1->seedno[j];
                                u_k++;
                            }
                            for(j=0; j<s_k&&j<SM; j++)
                            {
                                temp_list[u_k]=temp_spr->loczhi[j]+ZV;
                                temp_seedn[u_k]=temp_spr->seedno[j];
                                u_k++;
                            }
                        }
                        flag_end=find_location(temp_list,temp_seedn,temp_score,location_loc,u_k,&repeat_loc,BC,read_len);
                        if(flag_end==0)continue;
                        if(temp_score[repeat_loc]<6)continue;
                        canidate_temp.score=temp_score[repeat_loc];
                        loc_seed=temp_seedn[repeat_loc];
                        //loc_list=temp_list[repeat_loc];
                        location_loc[0]=start_loc+location_loc[0];
                        location_loc[1]=(location_loc[1]-1)*BC;
                        loc_list=location_loc[0];
                        left_length1=location_loc[0]+seed_len-1;
                        right_length1=seqcount-location_loc[0];
                        left_length2=location_loc[1]+seed_len-1;
                        right_length2=read_len-location_loc[1];
                        if(left_length1>=left_length2)num1=left_length2;
                        else num1=left_length1;
                        if(right_length1>=right_length2)num2=right_length2;
                        else num2=right_length1;
                        seedcount=0;
                        canidate_temp.loc1=location_loc[0];
                        canidate_temp.num1=num1;
                        canidate_temp.loc2=location_loc[1];
                        canidate_temp.num2=num2;
                        canidate_temp.left1=left_length1;
                        canidate_temp.left2=left_length2;
                        canidate_temp.right1=right_length1;
                        canidate_temp.right2=right_length2;
                        //find all left seed
                        for(u_k=*index_spr-2,k=num1/ZV,temp_spr1=temp_spr-2; u_k>=0&&k>=0; temp_spr1--,k--,u_k--)if(temp_spr1->score>0)
                            {
                                start_loc=u_k*ZV;
                                for(j=0,s_k=0; j<temp_spr1->score; j++)if(fabs((loc_list-start_loc-temp_spr1->loczhi[j])/((loc_seed-temp_spr1->seedno[j])*BC*1.0)-1.0)<0.2)
                                    {
                                        seedcount++;
                                        s_k++;
                                    }
                                if(s_k*1.0/temp_spr1->score>0.4)temp_spr1->score=0;
                            }
                        //find all right seed
                        for(u_k=*index_spr+1,k=num2/ZV,temp_spr1=temp_spr+1; k>0; temp_spr1++,k--,u_k++)if(temp_spr1->score>0)
                            {
                                start_loc=u_k*ZV;
                                for(j=0,s_k=0; j<temp_spr1->score; j++)if(fabs((start_loc+temp_spr1->loczhi[j]-loc_list)/((temp_spr1->seedno[j]-loc_seed)*BC*1.0)-1.0)<0.2)
                                    {
                                        seedcount++;
                                        s_k++;
                                    }
                                if(s_k*1.0/temp_spr1->score>0.4)temp_spr1->score=0;
                            }
                        canidate_temp.score=canidate_temp.score+seedcount;
                        if(ii==1)canidate_temp.chain='F';
                        else canidate_temp.chain='R';
                        //insert canidate position or delete this position
                        low=0;
                        high=canidatenum-1;
                        while(low<=high)
                        {
                            mid=(low+high)/2;
                            if(mid>=canidatenum||canidate_loc[mid].score<canidate_temp.score)high=mid-1;
                            else low=mid+1;
                        }
                        if(canidatenum<MAXC)for(u_k=canidatenum-1; u_k>high; u_k--)canidate_loc[u_k+1]=canidate_loc[u_k];
                        else for(u_k=canidatenum-2; u_k>high; u_k--)canidate_loc[u_k+1]=canidate_loc[u_k];
                        if(high+1<MAXC)canidate_loc[high+1]=canidate_temp;
                        if(canidatenum<MAXC)canidatenum++;
                        else canidatenum=MAXC;
                    }
                for(i=0,index_spr=index_list; i<cc1; i++,index_spr++)
                {
                    database[*index_spr].score=0;
                    database[*index_spr].index=-1;
                }
            }

            loc_flag=0;
            for(i=0; i<canidatenum; i++)
            {

                location_loc[0]=canidate_loc[i].loc1;
                location_loc[1]=canidate_loc[i].loc2;
                //readno=canidate_loc[i].readno;
                num1=canidate_loc[i].num1;
                num2=canidate_loc[i].num2;
                left_length1=canidate_loc[i].left1;
                left_length2=canidate_loc[i].left2;
                right_length1=canidate_loc[i].right1;
                right_length2=canidate_loc[i].right2;
                if(canidate_loc[i].chain=='F')onedata=onedata1;
                else if(canidate_loc[i].chain=='R')onedata=onedata2;
                //left alignment search
                seq_pr1=seq+location_loc[0]+seed_len-2;
                seq_pr2=onedata+location_loc[1]+seed_len-1;
                left_loc1=0;
                left_loc=0;
                resultstore->left_store1[0]='\0';
                resultstore->left_store2[0]='\0';
                flag_end=1;
                while(flag_end)
                {
                    if(num1>600)
                    {
                        for(s_k=0; s_k<DN; s_k++,seq_pr1--,seq_pr2--)
                        {
                            seq1[s_k]=*seq_pr1;
                            seq2[s_k]=*seq_pr2;
                        }
                        seq1[s_k]='\0';
                        seq2[s_k]='\0';
                    }
                    else
                    {
                        flag_end=0;
                        for(s_k=0; s_k<num1; s_k++,seq_pr1--,seq_pr2--)
                        {
                            seq1[s_k]=*seq_pr1;
                            seq2[s_k]=*seq_pr2;
                        }
                        seq1[s_k]='\0';
                        seq2[s_k]='\0';
                    }
                    for(loc=0; loc<2000; loc++)
                    {
                        longstr1[loc]=0;
                        longstr2[loc]=0;
                    }
                    align_flag=align(seq1,seq2,0.3*s_k,400,&strvalue,longstr1,longstr2,d_path,aln_path,deffvalue);
                    if(align_flag==1)
                    {
                        strvalue.q_aln_str[strvalue.aln_str_size]='\0';
                        strvalue.t_aln_str[strvalue.aln_str_size]='\0';
                        for(k=strvalue.aln_str_size-1,loc=0,sci=0,eit=0; k>-1&&eit<6; k--)
                        {
                            if(strvalue.q_aln_str[k]!='-')loc++;
                            if(strvalue.t_aln_str[k]!='-')sci++;
                            if(strvalue.q_aln_str[k]==strvalue.t_aln_str[k])eit++;
                            else eit=0;
                        }
                        if(flag_end==1)
                        {
                            loc=DN-strvalue.aln_q_e+loc;
                            sci=DN-strvalue.aln_t_e+sci;
                            if(loc==DN)align_flag=0;
                            seq_pr1=seq_pr1+loc;
                            seq_pr2=seq_pr2+sci;
                            strvalue.q_aln_str[k+1]='\0';
                            strvalue.t_aln_str[k+1]='\0';
                            left_loc1=left_loc1+DN-loc;
                            left_loc=left_loc+DN-sci;
                        }
                        else
                        {
                            loc=num1-strvalue.aln_q_e;
                            sci=num1-strvalue.aln_t_e;
                            if(loc==num1)align_flag=0;
                            left_loc1=left_loc1+num1-loc;
                            left_loc=left_loc+num1-sci;
                            seq_pr1=seq_pr1+loc+1;
                            seq_pr2=seq_pr2+sci+1;
                            //strvalue.q_aln_str[k+7]='\0';
                            //strvalue.t_aln_str[k+7]='\0';
                        }
                        strcat(resultstore->left_store1,strvalue.q_aln_str);
                        strcat(resultstore->left_store2,strvalue.t_aln_str);
                        if(left_length1-left_loc1>=left_length2-left_loc)num1=left_length2-left_loc;
                        else num1=left_length1-left_loc1;
                    }
                    else if(align_flag==2)
                    {
                        strvalue.q_aln_str[strvalue.aln_str_size]='\0';
                        strvalue.t_aln_str[strvalue.aln_str_size]='\0';
                        loc=DN-strvalue.aln_q_e;
                        sci=DN-strvalue.aln_t_e;
                        if(loc==DN)align_flag=0;
                        seq_pr1=seq_pr1+loc;
                        seq_pr2=seq_pr2+sci;
                        //strvalue.q_aln_str[k+1]='\0';
                        // strvalue.t_aln_str[k+1]='\0';
                        left_loc1=left_loc1+DN-loc;
                        left_loc=left_loc+DN-sci;
                        strcat(resultstore->left_store1,strvalue.q_aln_str);
                        strcat(resultstore->left_store2,strvalue.t_aln_str);
                    }
                    if(align_flag!=1)break;
                }

                //right alignment search
                right_loc1=0;
                right_loc=0;
                seq_pr1=seq+location_loc[0]-1;
                seq_pr2=onedata+location_loc[1];
                resultstore->right_store1[0]='\0';
                resultstore->right_store2[0]='\0';
                flag_end=1;
                while(flag_end)
                {
                    if(num2>600)
                    {
                        for(s_k=0; s_k<DN; s_k++,seq_pr1++,seq_pr2++)
                        {
                            seq1[s_k]=*seq_pr1;
                            seq2[s_k]=*seq_pr2;
                        }
                        seq1[s_k]='\0';
                        seq2[s_k]='\0';
                    }
                    else
                    {
                        flag_end=0;
                        for(s_k=0; s_k<num2; s_k++,seq_pr1++,seq_pr2++)
                        {
                            seq1[s_k]=*seq_pr1;
                            seq2[s_k]=*seq_pr2;
                        }
                        seq1[s_k]='\0';
                        seq2[s_k]='\0';
                    }
                    for(loc=0; loc<2000; loc++)
                    {
                        longstr1[loc]=0;
                        longstr2[loc]=0;
                    }
                    align_flag=align(seq1,seq2,0.3*s_k,400,&strvalue,longstr1,longstr2,d_path,aln_path,deffvalue);
                    if(align_flag==1)
                    {
                        strvalue.q_aln_str[strvalue.aln_str_size]='\0';
                        strvalue.t_aln_str[strvalue.aln_str_size]='\0';
                        for(k=strvalue.aln_str_size-1,loc=0,sci=0,eit=0; k>-1&&eit<6; k--)
                        {
                            if(strvalue.q_aln_str[k]!='-')loc++;
                            if(strvalue.t_aln_str[k]!='-')sci++;
                            if(strvalue.q_aln_str[k]==strvalue.t_aln_str[k])eit++;
                            else eit=0;
                        }
                        if(flag_end==1)
                        {
                            loc=DN-strvalue.aln_q_e+loc;
                            sci=DN-strvalue.aln_t_e+sci;
                            if(loc==DN)align_flag=0;
                            seq_pr1=seq_pr1-loc;
                            seq_pr2=seq_pr2-sci;
                            strvalue.q_aln_str[k+1]='\0';
                            strvalue.t_aln_str[k+1]='\0';
                            right_loc1=right_loc1+DN-loc;
                            right_loc=right_loc+DN-sci;
                        }
                        else
                        {
                            loc=num2-strvalue.aln_q_e;
                            sci=num2-strvalue.aln_t_e;
                            if(loc==num2)align_flag=0;
                            right_loc1=right_loc1+num2-loc;
                            right_loc=right_loc+num2-sci;
                            seq_pr1=seq_pr1-loc;
                            seq_pr2=seq_pr2-sci;
                            //strvalue.q_aln_str[k+7]='\0';
                            //strvalue.t_aln_str[k+7]='\0';
                        }
                        strcat(resultstore->right_store1,strvalue.q_aln_str);
                        strcat(resultstore->right_store2,strvalue.t_aln_str);
                        if(right_length1-right_loc1>=right_length2-right_loc)num2=right_length2-right_loc;
                        else num2=right_length1-right_loc1;
                    }
                    else if(align_flag==2)
                    {
                        strvalue.q_aln_str[strvalue.aln_str_size]='\0';
                        strvalue.t_aln_str[strvalue.aln_str_size]='\0';
                        loc=DN-strvalue.aln_q_e;
                        sci=DN-strvalue.aln_t_e;
                        if(loc==DN)align_flag=0;
                        seq_pr1=seq_pr1-loc;
                        seq_pr2=seq_pr2-sci;;
                        right_loc1=right_loc1+DN-loc;
                        right_loc=right_loc+DN-sci;
                        strcat(resultstore->right_store1,strvalue.q_aln_str);
                        strcat(resultstore->right_store2,strvalue.t_aln_str);
                    }
                    if(align_flag!=1)break;
                }

                s_k=strlen(resultstore->left_store1);
                for(j=0,loc=0,k=0; j<s_k; j++)
                {
                    if(resultstore->left_store1[j]!='-')
                    {
                        resultstore->out_store1[loc]=resultstore->left_store1[j];
                        loc++;
                    }
                    if(resultstore->left_store2[j]!='-')
                    {
                        resultstore->out_store2[k]=resultstore->left_store2[j];
                        k++;
                    }
                }
                resultstore->out_store1[loc]='\0';
                resultstore->out_store2[k]='\0';
                //printf("%s\n%s\n%s\n%s\n",resultstore->out_store1,resultstore->out_store2,resultstore->left_store1,resultstore->left_store2);
                string_check(resultstore->out_store1,resultstore->out_store2,resultstore->left_store1,resultstore->left_store2);
                //printf("%s\n%s\n",resultstore->left_store1,resultstore->left_store2);

                //output result
                u_k=strlen(resultstore->left_store1);
                for(j=u_k-1,loc=0,eit=0,k=0; j>-1; j--,k++)
                {
                    FR=resultstore->left_store1[j];
                    resultstore->out_store1[k]=FR;
                    if(FR!='-')loc++;
                    FR=resultstore->left_store2[j];
                    resultstore->out_store2[k]=FR;
                    if(FR!='-')eit++;
                }
                resultstore->out_store1[k]='\0';
                resultstore->out_store2[k]='\0';
                //resultstore->out_store1[k-seed_len]='\0';resultstore->out_store2[k-seed_len]='\0';
                if(u_k>0)
                {
                    left_loc1=location_loc[0]+seed_len-loc;
                    left_loc=location_loc[1]+seed_len-eit+1;
                }
                else
                {
                    left_loc1=location_loc[0];
                    left_loc=location_loc[1]+1;
                }

                s_k=strlen(resultstore->right_store1);
                for(k=0,loc=0,eit=0; k<s_k; k++)
                {
                    if(resultstore->right_store1[k]!='-')loc++;
                    if(resultstore->right_store2[k]!='-')eit++;
                }
                if(s_k>0)
                {
                    right_loc1=location_loc[0]+loc-1;
                    right_loc=location_loc[1]+eit;
                }
                else
                {
                    right_loc1=location_loc[0]+seed_len-1;
                    right_loc=location_loc[1]+seed_len;
                }
                if(s_k>=seed_len&&u_k>=seed_len)
                {
                    strcat(resultstore->out_store1,resultstore->right_store1+seed_len);
                    strcat(resultstore->out_store2,resultstore->right_store2+seed_len);
                }
                else if(u_k<seed_len)
                {
                    strcpy(resultstore->out_store1,resultstore->right_store1);
                    strcpy(resultstore->out_store2,resultstore->right_store2);
                }
                FR=canidate_loc[i].chain;
                if(strlen(resultstore->out_store1)>=1000)
                {
					set_temp_result(read_name, FR, canidate_loc[i].score, left_loc, 
									right_loc, read_len, left_loc1, right_loc1, 
									resultstore->out_store2,resultstore->out_store1, &result);
					output_temp_result(&result, outfile[threadint]);
                    loc_flag=1;
                }
            }//search result
            //when the above search fail,the program should be the following accuracy search
            if(loc_flag==0)
            {
                canidatenum=0;
                for(ii=1; ii<=2; ii++)
                {
                    read_len=strlen(onedata1);
                    BC=5;
                    if(ii==1)onedata=onedata1;
                    else if(ii==2)
                    {
                        strcpy(onedata2,onedata1);
                        onedata=onedata2;
                        for(j=read_len-1,i=0; j>i; j--,i++)
                        {
                            FR=onedata[i];
                            onedata[i]=onedata[j];
                            onedata[j]=FR;
                        }
                        for(i=0; i<read_len; i++)
                        {
                            FR=onedata[i];
                            switch(FR)
                            {
                            case 'A':
                            {
                                onedata[i]='T';
                                break;
                            }
                            case 'T':
                            {
                                onedata[i]='A';
                                break;
                            }
                            case 'C':
                            {
                                onedata[i]='G';
                                break;
                            }
                            case 'G':
                            {
                                onedata[i]='C';
                                break;
                            }
                            }
                        }
                    }

                    loc_flag=0;
                    endnum=0;
                    read_len=strlen(onedata);
                    cleave_num=transnum_buchang(onedata,mvalue,&endnum,read_len,seed_len,BC);
                    j=0;
                    index_spr=index_list;
                    index_ss=index_score;
                    endnum=0;
                    for(k=0; k<cleave_num; k++)if(mvalue[k]>=0)
                        {
                            count1=countin[mvalue[k]];
                            //if(count1>20)continue;
                            leadarray=databaseindex[mvalue[k]];
                            for(i=0; i<count1; i++,leadarray++)
                            {
                                templong=(*leadarray)/ZVS;
                                u_k=(*leadarray)%ZVS;
                                if(templong>=0)
                                {
                                    temp_spr=database+templong;
                                    if(temp_spr->score==0||temp_spr->seednum<k+1)
                                    {
                                        loc=++(temp_spr->score);
                                        if(loc<=SM)
                                        {
                                            temp_spr->loczhi[loc-1]=u_k;
                                            temp_spr->seedno[loc-1]=k+1;
                                        }
                                        else insert_loc(temp_spr,u_k,k+1,BC);
                                        if(templong>0)s_k=temp_spr->score+(temp_spr-1)->score;
                                        else s_k=temp_spr->score;
                                        if(endnum<s_k)endnum=s_k;
                                        if(temp_spr->index==-1)
                                        {
                                            *(index_spr++)=templong;
                                            *(index_ss++)=s_k;
                                            temp_spr->index=j;
                                            j++;
                                        }
                                        else index_score[temp_spr->index]=s_k;
                                    }
                                    temp_spr->seednum=k+1;
                                }
                            }
                        }
                    cc1=j;
                    for(i=0,index_spr=index_list,index_ss=index_score; i<cc1; i++,index_spr++,index_ss++)if(*index_ss>4)
                        {
                            temp_spr=database+*index_spr;
                            if(temp_spr->score==0)continue;
                            s_k=temp_spr->score;
							if(*index_spr>0)loc=(temp_spr-1)->score;
                            else loc=0;
                            start_loc=(*index_spr)*ZVS;
                            if(*index_spr>0)
                            {
                                loc=(temp_spr-1)->score;
                                if(loc>0)start_loc=(*index_spr-1)*ZVS;
                            }
                            else loc=0;
                            if(loc==0)for(j=0,u_k=0; j<s_k&&j<SM; j++)
                                {
                                    temp_list[u_k]=temp_spr->loczhi[j];
                                    temp_seedn[u_k]=temp_spr->seedno[j];
                                    u_k++;
                                }
                            else
                            {
                                k=loc;
                                u_k=0;
                                temp_spr1=temp_spr-1;
                                for(j=0; j<k&&j<SM; j++)
                                {
                                    temp_list[u_k]=temp_spr1->loczhi[j];
                                    temp_seedn[u_k]=temp_spr1->seedno[j];
                                    u_k++;
                                }
                                for(j=0; j<s_k&&j<SM; j++)
                                {
                                    temp_list[u_k]=temp_spr->loczhi[j]+ZVS;
                                    temp_seedn[u_k]=temp_spr->seedno[j];
                                    u_k++;
                                }
                            }
                            flag_end=find_location(temp_list,temp_seedn,temp_score,location_loc,u_k,&repeat_loc,BC,read_len);
                            if(flag_end==0)continue;
                            if(temp_score[repeat_loc]<6)continue;
                            canidate_temp.score=temp_score[repeat_loc];
                            loc_seed=temp_seedn[repeat_loc];
                            //loc_list=temp_list[repeat_loc];
                            location_loc[0]=start_loc+location_loc[0];
                            location_loc[1]=(location_loc[1]-1)*BC;
                            loc_list=location_loc[0];
                            left_length1=location_loc[0]+seed_len-1;
                            right_length1=seqcount-location_loc[0];
                            left_length2=location_loc[1]+seed_len-1;
                            right_length2=read_len-location_loc[1];
                            if(left_length1>=left_length2)num1=left_length2;
                            else num1=left_length1;
                            if(right_length1>=right_length2)num2=right_length2;
                            else num2=right_length1;
                            seedcount=0;
                            canidate_temp.loc1=location_loc[0];
                            canidate_temp.num1=num1;
                            canidate_temp.loc2=location_loc[1];
                            canidate_temp.num2=num2;
                            canidate_temp.left1=left_length1;
                            canidate_temp.left2=left_length2;
                            canidate_temp.right1=right_length1;
                            canidate_temp.right2=right_length2;
                            //find all left seed
                            for(u_k=*index_spr-2,k=num1/ZVS,temp_spr1=temp_spr-2; u_k>=0&&k>=0; temp_spr1--,k--,u_k--)if(temp_spr1->score>0)
                                {
                                    start_loc=u_k*ZVS;
                                    for(j=0,s_k=0; j<temp_spr1->score; j++)if(fabs((loc_list-start_loc-temp_spr1->loczhi[j])/((loc_seed-temp_spr1->seedno[j])*BC*1.0)-1.0)<0.2)
                                        {
                                            seedcount++;
                                            s_k++;
                                        }
                                    if(s_k*1.0/temp_spr1->score>0.4)temp_spr1->score=0;
                                }
                            //find all right seed
                            for(u_k=*index_spr+1,k=num2/ZVS,temp_spr1=temp_spr+1; k>0; temp_spr1++,k--,u_k++)if(temp_spr1->score>0)
                                {
                                    start_loc=u_k*ZVS;
                                    for(j=0,s_k=0; j<temp_spr1->score; j++)if(fabs((start_loc+temp_spr1->loczhi[j]-loc_list)/((temp_spr1->seedno[j]-loc_seed)*BC*1.0)-1.0)<0.2)
                                        {
                                            seedcount++;
                                            s_k++;
                                        }
                                    if(s_k*1.0/temp_spr1->score>0.4)temp_spr1->score=0;
                                }
                            canidate_temp.score=canidate_temp.score+seedcount;
                            if(ii==1)canidate_temp.chain='F';
                            else canidate_temp.chain='R';
                            //insert canidate position or delete this position
                            low=0;
                            high=canidatenum-1;
                            while(low<=high)
                            {
                                mid=(low+high)/2;
                                if(mid>=canidatenum||canidate_loc[mid].score<canidate_temp.score)high=mid-1;
                                else low=mid+1;
                            }
                            if(canidatenum<MAXC)for(u_k=canidatenum-1; u_k>high; u_k--)canidate_loc[u_k+1]=canidate_loc[u_k];
                            else for(u_k=canidatenum-2; u_k>high; u_k--)canidate_loc[u_k+1]=canidate_loc[u_k];
                            if(high+1<MAXC)canidate_loc[high+1]=canidate_temp;
                            if(canidatenum<MAXC)canidatenum++;
                            else canidatenum=MAXC;
                        }
                    for(i=0,index_spr=index_list; i<cc1; i++,index_spr++)
                    {
                        database[*index_spr].score=0;
                        database[*index_spr].index=-1;
                    }
                }

                for(i=0; i<canidatenum; i++)
                {

                    location_loc[0]=canidate_loc[i].loc1;
                    location_loc[1]=canidate_loc[i].loc2;
                    num1=canidate_loc[i].num1;
                    num2=canidate_loc[i].num2;
                    left_length1=canidate_loc[i].left1;
                    left_length2=canidate_loc[i].left2;
                    right_length1=canidate_loc[i].right1;
                    right_length2=canidate_loc[i].right2;
                    if(canidate_loc[i].chain=='F')onedata=onedata1;
                    else if(canidate_loc[i].chain=='R')onedata=onedata2;
                    //left alignment search
                    seq_pr1=seq+location_loc[0]+seed_len-2;
                    seq_pr2=onedata+location_loc[1]+seed_len-1;
                    left_loc1=0;
                    left_loc=0;
                    resultstore->left_store1[0]='\0';
                    resultstore->left_store2[0]='\0';
                    flag_end=1;
                    while(flag_end)
                    {
                        if(num1>600)
                        {
                            for(s_k=0; s_k<DN; s_k++,seq_pr1--,seq_pr2--)
                            {
                                seq1[s_k]=*seq_pr1;
                                seq2[s_k]=*seq_pr2;
                            }
                            seq1[s_k]='\0';
                            seq2[s_k]='\0';
                        }
                        else
                        {
                            flag_end=0;
                            for(s_k=0; s_k<num1; s_k++,seq_pr1--,seq_pr2--)
                            {
                                seq1[s_k]=*seq_pr1;
                                seq2[s_k]=*seq_pr2;
                            }
                            seq1[s_k]='\0';
                            seq2[s_k]='\0';
                        }
                        for(loc=0; loc<2000; loc++)
                        {
                            longstr1[loc]=0;
                            longstr2[loc]=0;
                        }
                        align_flag=align(seq1,seq2,0.3*s_k,400,&strvalue,longstr1,longstr2,d_path,aln_path,deffvalue);
                        if(align_flag==1)
                        {
                            strvalue.q_aln_str[strvalue.aln_str_size]='\0';
                            strvalue.t_aln_str[strvalue.aln_str_size]='\0';
                            for(k=strvalue.aln_str_size-1,loc=0,sci=0,eit=0; k>-1&&eit<6; k--)
                            {
                                if(strvalue.q_aln_str[k]!='-')loc++;
                                if(strvalue.t_aln_str[k]!='-')sci++;
                                if(strvalue.q_aln_str[k]==strvalue.t_aln_str[k])eit++;
                                else eit=0;
                            }
                            if(flag_end==1)
                            {
                                loc=DN-strvalue.aln_q_e+loc;
                                sci=DN-strvalue.aln_t_e+sci;
                                if(loc==DN)align_flag=0;
                                seq_pr1=seq_pr1+loc;
                                seq_pr2=seq_pr2+sci;
                                strvalue.q_aln_str[k+1]='\0';
                                strvalue.t_aln_str[k+1]='\0';
                                left_loc1=left_loc1+DN-loc;
                                left_loc=left_loc+DN-sci;
                            }
                            else
                            {
                                loc=num1-strvalue.aln_q_e;
                                sci=num1-strvalue.aln_t_e;
                                if(loc==num1)align_flag=0;
                                left_loc1=left_loc1+num1-loc;
                                left_loc=left_loc+num1-sci;
                                seq_pr1=seq_pr1+loc+1;
                                seq_pr2=seq_pr2+sci+1;
                                //strvalue.q_aln_str[k+7]='\0';
                                //strvalue.t_aln_str[k+7]='\0';
                            }
                            strcat(resultstore->left_store1,strvalue.q_aln_str);
                            strcat(resultstore->left_store2,strvalue.t_aln_str);
                            if(left_length1-left_loc1>=left_length2-left_loc)num1=left_length2-left_loc;
                            else num1=left_length1-left_loc1;
                        }
                        else if(align_flag==2)
                        {
                            strvalue.q_aln_str[strvalue.aln_str_size]='\0';
                            strvalue.t_aln_str[strvalue.aln_str_size]='\0';
                            loc=DN-strvalue.aln_q_e;
                            sci=DN-strvalue.aln_t_e;
                            if(loc==DN)align_flag=0;
                            seq_pr1=seq_pr1+loc;
                            seq_pr2=seq_pr2+sci;
                            //strvalue.q_aln_str[k+1]='\0';
                            // strvalue.t_aln_str[k+1]='\0';
                            left_loc1=left_loc1+DN-loc;
                            left_loc=left_loc+DN-sci;
                            strcat(resultstore->left_store1,strvalue.q_aln_str);
                            strcat(resultstore->left_store2,strvalue.t_aln_str);
                        }
                        if(align_flag!=1)break;
                    }

                    //right alignment search
                    right_loc1=0;
                    right_loc=0;
                    seq_pr1=seq+location_loc[0]-1;
                    seq_pr2=onedata+location_loc[1];
                    resultstore->right_store1[0]='\0';
                    resultstore->right_store2[0]='\0';
                    flag_end=1;
                    while(flag_end)
                    {
                        if(num2>600)
                        {
                            for(s_k=0; s_k<DN; s_k++,seq_pr1++,seq_pr2++)
                            {
                                seq1[s_k]=*seq_pr1;
                                seq2[s_k]=*seq_pr2;
                            }
                            seq1[s_k]='\0';
                            seq2[s_k]='\0';
                        }
                        else
                        {
                            flag_end=0;
                            for(s_k=0; s_k<num2; s_k++,seq_pr1++,seq_pr2++)
                            {
                                seq1[s_k]=*seq_pr1;
                                seq2[s_k]=*seq_pr2;
                            }
                            seq1[s_k]='\0';
                            seq2[s_k]='\0';
                        }
                        for(loc=0; loc<2000; loc++)
                        {
                            longstr1[loc]=0;
                            longstr2[loc]=0;
                        }
                        align_flag=align(seq1,seq2,0.3*s_k,400,&strvalue,longstr1,longstr2,d_path,aln_path,deffvalue);
                        if(align_flag==1)
                        {
                            strvalue.q_aln_str[strvalue.aln_str_size]='\0';
                            strvalue.t_aln_str[strvalue.aln_str_size]='\0';
                            for(k=strvalue.aln_str_size-1,loc=0,sci=0,eit=0; k>-1&&eit<6; k--)
                            {
                                if(strvalue.q_aln_str[k]!='-')loc++;
                                if(strvalue.t_aln_str[k]!='-')sci++;
                                if(strvalue.q_aln_str[k]==strvalue.t_aln_str[k])eit++;
                                else eit=0;
                            }
                            if(flag_end==1)
                            {
                                loc=DN-strvalue.aln_q_e+loc;
                                sci=DN-strvalue.aln_t_e+sci;
                                if(loc==DN)align_flag=0;
                                seq_pr1=seq_pr1-loc;
                                seq_pr2=seq_pr2-sci;
                                strvalue.q_aln_str[k+1]='\0';
                                strvalue.t_aln_str[k+1]='\0';
                                right_loc1=right_loc1+DN-loc;
                                right_loc=right_loc+DN-sci;
                            }
                            else
                            {
                                loc=num2-strvalue.aln_q_e;
                                sci=num2-strvalue.aln_t_e;
                                if(loc==num2)align_flag=0;
                                right_loc1=right_loc1+num2-loc;
                                right_loc=right_loc+num2-sci;
                                seq_pr1=seq_pr1-loc;
                                seq_pr2=seq_pr2-sci;
                                //strvalue.q_aln_str[k+7]='\0';
                                //strvalue.t_aln_str[k+7]='\0';
                            }
                            strcat(resultstore->right_store1,strvalue.q_aln_str);
                            strcat(resultstore->right_store2,strvalue.t_aln_str);
                            if(right_length1-right_loc1>=right_length2-right_loc)num2=right_length2-right_loc;
                            else num2=right_length1-right_loc1;
                        }
                        else if(align_flag==2)
                        {
                            strvalue.q_aln_str[strvalue.aln_str_size]='\0';
                            strvalue.t_aln_str[strvalue.aln_str_size]='\0';
                            loc=DN-strvalue.aln_q_e;
                            sci=DN-strvalue.aln_t_e;
                            if(loc==DN)align_flag=0;
                            seq_pr1=seq_pr1-loc;
                            seq_pr2=seq_pr2-sci;;
                            right_loc1=right_loc1+DN-loc;
                            right_loc=right_loc+DN-sci;
                            strcat(resultstore->right_store1,strvalue.q_aln_str);
                            strcat(resultstore->right_store2,strvalue.t_aln_str);
                        }
                        if(align_flag!=1)break;
                    }

                    s_k=strlen(resultstore->left_store1);
                    for(j=0,loc=0,k=0; j<s_k; j++)
                    {
                        if(resultstore->left_store1[j]!='-')
                        {
                            resultstore->out_store1[loc]=resultstore->left_store1[j];
                            loc++;
                        }
                        if(resultstore->left_store2[j]!='-')
                        {
                            resultstore->out_store2[k]=resultstore->left_store2[j];
                            k++;
                        }
                    }
                    resultstore->out_store1[loc]='\0';
                    resultstore->out_store2[k]='\0';
                    //printf("%s\n%s\n%s\n%s\n",resultstore->out_store1,resultstore->out_store2,resultstore->left_store1,resultstore->left_store2);
                    string_check(resultstore->out_store1,resultstore->out_store2,resultstore->left_store1,resultstore->left_store2);
                    //printf("%s\n%s\n",resultstore->left_store1,resultstore->left_store2);

                    //output result
                    u_k=strlen(resultstore->left_store1);
                    for(j=u_k-1,loc=0,eit=0,k=0; j>-1; j--,k++)
                    {
                        FR=resultstore->left_store1[j];
                        resultstore->out_store1[k]=FR;
                        if(FR!='-')loc++;
                        FR=resultstore->left_store2[j];
                        resultstore->out_store2[k]=FR;
                        if(FR!='-')eit++;
                    }
                    resultstore->out_store1[k]='\0';
                    resultstore->out_store2[k]='\0';
                    //resultstore->out_store1[k-seed_len]='\0';resultstore->out_store2[k-seed_len]='\0';
                    if(u_k>0)
                    {
                        left_loc1=location_loc[0]+seed_len-loc;
                        left_loc=location_loc[1]+seed_len-eit+1;
                    }
                    else
                    {
                        left_loc1=location_loc[0];
                        left_loc=location_loc[1]+1;
                    }

                    s_k=strlen(resultstore->right_store1);
                    for(k=0,loc=0,eit=0; k<s_k; k++)
                    {
                        if(resultstore->right_store1[k]!='-')loc++;
                        if(resultstore->right_store2[k]!='-')eit++;
                    }
                    if(s_k>0)
                    {
                        right_loc1=location_loc[0]+loc-1;
                        right_loc=location_loc[1]+eit;
                    }
                    else
                    {
                        right_loc1=location_loc[0]+seed_len-1;
                        right_loc=location_loc[1]+seed_len;
                    }
                    if(s_k>=seed_len&&u_k>=seed_len)
                    {
                        strcat(resultstore->out_store1,resultstore->right_store1+seed_len);
                        strcat(resultstore->out_store2,resultstore->right_store2+seed_len);
                    }
                    else if(u_k<seed_len)
                    {
                        strcpy(resultstore->out_store1,resultstore->right_store1);
                        strcpy(resultstore->out_store2,resultstore->right_store2);
                    }
                    FR=canidate_loc[i].chain;
                    if(strlen(resultstore->out_store1)>1000)
                    {
						set_temp_result(read_name, FR, canidate_loc[i].score, left_loc, 
										right_loc, read_len, left_loc1, right_loc1, 
										resultstore->out_store2,resultstore->out_store1, &result);
						output_temp_result(&result, outfile[threadint]);
                        loc_flag=1;
                    }
                }//search result
            }
        }
    }
    free(database);
    free(index_list);
    free(index_score);
    free(d_path);
    free(resultstore);
}


static void* multithread(void* arg)
{
    int localthreadno;
    pthread_mutex_lock(&mutilock);
    localthreadno=runthreadnum;
    runthreadnum++;
    pthread_mutex_unlock(&mutilock);
    reference_mapping(localthreadno);
	
	return NULL;
}

static int load_fastq(FILE *fq)
{
    int readlen,readno,sum=0,flag;
    char *pre;
    readcount=0;
    pre=savework;
    while((flag=fscanf(fq,"%d\t%d\t%s\n",&readno,&readlen,pre))!=EOF&&readcount<SVM&&sum<MAXSTR)
    {
        readinfo[readcount].seqloc=pre;
        readinfo[readcount].readno=readno;
        readlen=strlen(pre);
        readinfo[readcount].readlen=readlen;
        sum=sum+readlen+1;
        pre=pre+readlen+1;
        readcount++;
    }
    if(flag!=EOF)
    {
        readinfo[readcount].seqloc=pre;
        readinfo[readcount].readno=readno;
        readlen=strlen(pre);
        readinfo[readcount].readlen=readlen;
        readcount++;
        return(1);
    }
    else return(0);
}

int meap_ref_impl_small(int maxc)
{
	MAXC = maxc;
    char tempstr[300],fastafile[300];
    int corenum,threadflag,readall;
    int fileflag,threadno;
    FILE *fp,*fastq;
    struct timeval tpstart, tpend;
    float timeuse;
    fp=fopen("config.txt","r");
    assert(fscanf(fp,"%s\n%s\n%s\n%s\n%d %d\n",workpath,fastafile,fastqfile,tempstr,&corenum,&readall) == 6);
    fclose(fp);
    threadnum=corenum;
    //building reference index
    gettimeofday(&tpstart, NULL);
    seed_len=13;
    creat_ref_index(fastafile);
    gettimeofday(&tpend, NULL);
    timeuse = 1000000 * (tpend.tv_sec - tpstart.tv_sec) + tpend.tv_usec - tpstart.tv_usec;
    timeuse /= 1000000;
    fp = fopen("config.txt", "a");
    fprintf(fp, "The Building Reference Index Time: %f sec\n", timeuse);
    fclose(fp);
    gettimeofday(&tpstart, NULL);

    savework=(char *)malloc((MAXSTR+RM)*sizeof(char));
    readinfo=(ReadFasta*)malloc((SVM+2)*sizeof(ReadFasta));
    thread=(pthread_t*)malloc(threadnum*sizeof(pthread_t));
    outfile=(FILE **)malloc(threadnum*sizeof(FILE *));
    for(threadno=0; threadno<threadnum; threadno++)
    {
        sprintf(tempstr,"%s/%d.r",workpath,threadno+1);
        outfile[threadno]=fopen(tempstr,"w");
    }
    sprintf(tempstr,"%s/0.fq",workpath);
    fastq=fopen(tempstr,"r");
    //multi process thread
    fileflag=1;
    while(fileflag)
    {
        fileflag=load_fastq(fastq);
        if(readcount%PLL==0)terminalnum=readcount/PLL;
        else terminalnum=readcount/PLL+1;
        if(readcount<=0)break;
        runnumber=0;
        runthreadnum=0;
        pthread_mutex_init(&mutilock,NULL);
        //creat thread
        if(readcount>0)
        {
            for(threadno=0; threadno<threadnum; threadno++)
            {
                threadflag= pthread_create(&thread[threadno], NULL, multithread, NULL);
                if(threadflag)
                {
                    printf("ERROR; return code is %d\n", threadflag);
                    return EXIT_FAILURE;
                }
            }
            //waiting thread
            for(threadno=0; threadno<threadnum; threadno++)pthread_join(thread[threadno],NULL);
        }
    }
    fclose(fastq);
    //clear creat index memory
    free(countin);
    free(databaseindex);
    free(allloc);
    free(REFSEQ);

    gettimeofday(&tpend, NULL);
    timeuse = 1000000 * (tpend.tv_sec - tpstart.tv_sec) + tpend.tv_usec - tpstart.tv_usec;
    timeuse /= 1000000;
    fp = fopen("config.txt", "a");
    fprintf(fp, "The Mapping Time: %f sec\n", timeuse);
    fclose(fp);

    for(threadno=0; threadno<threadnum; threadno++)fclose(outfile[threadno]);
    free(outfile);
    free(savework);
    free(readinfo);
    free(thread);
    return 0;
}
