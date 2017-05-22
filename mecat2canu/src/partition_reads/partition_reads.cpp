#include <stdint.h>
#include <iostream>
#include <fstream>
#include <string>

#include <cstdlib>

using namespace std;

void
print_usage(const char* prog)
{
    cerr << "usage:\n"
         << prog << " " << "input" << " " << "output" << "\n";
}

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        cerr << "input error!\n";
        print_usage(argv[0]);
        return 1;
    }

    const char* input = argv[1];
    const char* output = argv[2];
    const int64_t vs = 2000000000;
    cerr << "vs = " << vs << "\n";
    int b = 1, n = 0;
    int64_t s = 0;

    ifstream in(input);
    if (!in)
    {
        cerr << "fail to open file \'" << input << "\' for reading.\n";
        return 1;
    }
    ofstream out(output);
    if (!out)
    {
        cerr << "fail to open file \'" << output << "\' for writing.\n";
        return 1;
    }

    string line;
    while (getline(in, line))
    {
        if (line[0] != '>' && line[0] != '@')
        {
            cerr << line.substr(0, 20) << "\n";
            cerr << "file \'" << input << "\' is not a FASTA or FASTQ format file!\n";
            abort();
        }
        const char first_header_char = line[0];
        getline(in, line);
        int64_t ss = line.size();
        if (s + ss > vs)
        {
            out << b << "\t" << b + n - 1 << "\n";
            s = 0;
            b = b + n;
            n = 0;
        }

        s += ss;
        ++n;

        if (first_header_char == '@')
        {
            getline(in, line);
            getline(in, line);
        }
    }

    if (s > 0)
    {
        out << b << "\t" << b + n - 1 << "\n";
        s = 0;
        b = b + n;
        n = 0;
    }

    in.close();
    out.close();

    return 0;
}
