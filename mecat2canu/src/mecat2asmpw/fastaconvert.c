#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<string.h>
#define RM 150000
void main(int argc, char*argv[]){
	FILE *fp, *ot;
    int kk = 0,i,read_len,j,k;
    char read_name[1000], onedata[RM], *fq, tempstr[200], *oq,ch;
	char dataset[300],output[300];
	for(i=1;i<argc;i++){
		k=strlen(argv[i]);
		for(j=2;j<k;j++)tempstr[j-2]=argv[i][j];tempstr[j-2]='\0';
		//printf("%s\n",tempstr);
		if(argv[i][0]=='-'){
			switch(argv[i][1]){
			case 'D':{strcpy(dataset,tempstr); break;}
			case 'O':{strcpy(output,tempstr);  break;}
			}
		}
	}
    fp = fopen(dataset, "r");
    ot = fopen(output, "w");
    fq = (char *)malloc(100000000);
    setvbuf(fp, fq, _IOFBF, 100000000);
    oq = (char *)malloc(100000000);
    setvbuf(ot, oq, _IOFBF, 100000000);
	ch=getc(fp);
	if(ch=='>'){
		kk=0;
       for (; ch!=EOF; ch=getc(fp)){
        if(ch=='>'){
			if(kk>0){
				onedata[read_len]='\0';
				fprintf(ot, ">%s\n%s\n",read_name,onedata);
				fscanf(fp,"%[^\n]s",read_name);
				read_len=0;
				kk++;
			}
			else {
				 fscanf(fp,"%[^\n]s",read_name);
				  read_len=0;
				  kk++;
			}
		}
        else if(ch!='\n'&&ch!='\r')onedata[read_len++]=ch;
       }
	   onedata[read_len]='\0';
	   fprintf(ot, ">%s\n%s\n",read_name,onedata);
	}
   fclose(fp);
   fclose(ot);
   free(fq);
   free(oq);
}
