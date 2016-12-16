#include <iostream>
#include <cstring>
#include <string>

using namespace std;


void
print_usage(const char* prog)
{
	const char delim = ' ';
	cerr << "usage:"
		 << "\n"
		 << prog
		 << delim
		 << "inputReads"
		 << delim
		 << "outputReads-prefix"
		 << delim
		 << "genomeSize"
		 << delim
		 << "coverage"
		 << "\n";
}

#define __run_system(cmd) \
	do { \
	int __rc_status = system(cmd); \
	if (__rc_status != 0) { fprintf(stderr, "[%s, %u] system() error. Error code is %d.\n", __func__, __LINE__, __rc_status); exit(1); } \
} while (0);


int main(int argc, char* argv[])
{
	if (argc != 5)
	{
		cerr << "Input Error!\n";
		print_usage(argv[0]);
		return 1;
	}
	
	const char* input_reads = argv[1];
	const char* output_reads = argv[2];
	long long genome_size = atoll(argv[3]);
	long long coverage = atoll(argv[4]);
	char cmd[2048];
	
	// fasta to fastq
	cerr << "step 1: convert fasta to fastq\n";
	string fastq = input_reads;
	fastq += ".fastq";
	sprintf(cmd, "es_fasta2fastq %s %s", input_reads, fastq.c_str());
	__run_system(cmd);
	
	// fastq to ca
	cerr << "step 2: convert fastq to CA\n";
	string libraryname = fastq + ".libname";
	string librarynamefrg = libraryname + ".frg";
	sprintf(cmd, "es_fastq2ca -libraryname %s -technology pacbio-corrected -type sanger -reads %s > %s",
			libraryname.c_str(), fastq.c_str(), librarynamefrg.c_str());
	__run_system(cmd);
	
	// gatekeeper 1
	cerr << "step 3\n";
	string librarynamegkpstore = libraryname + ".gkpStore";
	sprintf(cmd, "es_gatekeeper -T -F -o %s %s", librarynamegkpstore.c_str(), librarynamefrg.c_str());
	__run_system(cmd);
	
	// gatekeeper 2
	cerr << "step 4\n";
	long long totalbase = genome_size * coverage;
	string outputfrg = output_reads;
	outputfrg += ".frg";
	sprintf(cmd, "es_gatekeeper -dumpfrg -longestlength 0 %lld %s > %s", totalbase, librarynamegkpstore.c_str(), outputfrg.c_str());
	__run_system(cmd);
	
	// gatekeeper 3
	cerr << "step 5\n";
	sprintf(cmd, "es_gatekeeper -dumpfasta %s -longestlength 0 %lld %s",
		   output_reads, totalbase, librarynamegkpstore.c_str());
	__run_system(cmd);
	
	sprintf(cmd, "rm -rf %s.*", input_reads);
	__run_system(cmd);
	
    return 0;
}

