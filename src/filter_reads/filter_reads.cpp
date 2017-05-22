#include "../common/fasta_reader.h"

#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

typedef index_t idx;

void print_usage(const char* prog)
{
	const char sep = ' ';
	cerr << "USAGE:\n"
		 << prog << sep
		 << "input" << sep
		 << "output" << sep
		 << "min-length" << sep
		 << "max-length" << endl;
}

idx parse_int(const char* arg)
{
	istringstream in(arg);
	idx n;
	in >> n;
	if (!in) ERROR("failed to transfer '%s' to integer.", arg);
	return n;
}

void
output_one_read(Sequence& read, int& id, const idx min_size, const idx max_size, ostream& out)
{
	idx n = read.size();
	idx l = 0, r, left = n;
	Sequence::str_t& seq = read.sequence();
	while (left) {
		r = min(n, l + max_size);
		idx s = r - l;
		if (s < min_size) break;
		out << ">" << id++ << "\n";
		for (idx i = l; i < r; ++i) {
			out << seq[i];
		}
		out << "\n";
		l = r;
		left -= s;
	}
}

int main(int argc, char* argv[])
{
	if (argc != 5) {
		print_usage(argv[0]);
		exit(1);
	}
	const char* input = argv[1];
	const char* output = argv[2];
	const idx min_size = parse_int(argv[3]);
	const idx max_size = parse_int(argv[4]);
	
	FastaReader reader(input);
	Sequence read;
	ofstream out;
	open_fstream(out, output, ios::out);
	int id = 0;
	while (1) {
		idx s = reader.read_one_seq(read);
		if (s == -1) break;
		output_one_read(read, id, min_size, max_size, out);
	}
	close_fstream(out);
}
