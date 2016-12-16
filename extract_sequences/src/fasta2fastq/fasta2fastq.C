#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <time.h>

#define RM 100000


int main(int argc, char* argv[]) {
	int i,j,k,readcount,readall;
	char tempstr[300];
	const char* dataset = argv[1];
	const char* refset = argv[2];
	FILE *fp,*ot;
	char read_name[300],onedata[RM],buff2[RM],*fq,*oq;	
	fp=fopen(dataset,"r");
	ot=fopen(refset,"w");
	fq=(char *)malloc(100000000); 
	setvbuf(fp,fq,_IOFBF, 100000000);
	oq=(char *)malloc(100000000); 
	setvbuf(ot,oq,_IOFBF, 100000000);
        readcount=1,readall=1;
	while(fscanf(fp,">%[^\n]s",read_name)!=EOF&&fscanf(fp,"%s\n",onedata)!=EOF){
		k=strlen(onedata);
		for(i=0;i<k;i++)buff2[i]=']';
		buff2[k]='\0';
        fprintf(ot,"@%s\n%s\n+\n%s\n",read_name,onedata,buff2);
               //  printf("%d\n",readcount);
	}
	fclose(fp);
	fclose(ot);
	free(fq);free(oq);
    return 0;
}
