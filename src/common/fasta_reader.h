#ifndef FASTA_READER_H
#define FASTA_READER_H

#include "buffer_line_iterator.h"
#include "sequence.h"

class FastaReader
{
public:
    typedef Sequence::str_t                 str_t;
    typedef BufferLineReader::OneDataLine   OneDataLine;

public:
    FastaReader(const char* fasta_file_name) : m_Reader(fasta_file_name) { encode_table = get_dna_encode_table(); }
    idx_t read_one_seq(Sequence& seq);

private:
    void x_parse_defline(const OneDataLine& line, str_t& header);
    void x_parse_data_line(const OneDataLine& line, str_t& seq);
    void x_check_data_line(const OneDataLine& line);
    bool is_header_line(const OneDataLine& line)
    { return line.size() > 0 && (line.front() == '>' || line.front() == '@'); }
    bool is_nucl(const unsigned char ch)
    {
        int r = encode_table[ch];
        return r < 16;
    }
    bool is_ambig_nucl(const unsigned char ch)
    {
        int r = encode_table[ch];
        return (r < 16) && (r > 3);
    }
    bool is_upper_case_letter(const char ch)
    {
        return ch >= 'A' && ch <= 'Z';
    }
    bool is_lower_case_letter(const char ch)
    {
        return ch >= 'a' && ch <= 'z';
    }
    bool is_alpha(const unsigned char c)
    {
        return is_upper_case_letter(c) || is_lower_case_letter(c);
    }
    bool is_comment_line(const OneDataLine& line)
    { return line.front() == '#' || line.front() == '!'; }


private:
    BufferLineReader    m_Reader;
	const u1_t*		encode_table;
};

#endif // FASTA_READER_H
