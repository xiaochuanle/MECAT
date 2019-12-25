#ifndef FASTA_READER_H
#define FASTA_READER_H

#include "buffer_line_iterator.h"
#include "sequence.h"

class FastaReader {
    public:
	typedef Sequence::str_t str_t;
	typedef BufferLineReader::OneDataLine OneDataLine;
    public:
	explicit FastaReader(const char* const fasta_file_name) : m_Reader(fasta_file_name), encode_table(get_dna_encode_table()) { }
	idx_t read_one_seq(Sequence& seq);
	std::streampos tellg() const {
		return m_Reader.tellg();
	}
	void seekg(std::streampos pos) {
		m_Reader.seekg(pos);
	}
    private:
	void x_parse_defline(const OneDataLine& line, str_t& header);
	void x_parse_data_line(const OneDataLine& line, str_t& seq);
	void x_check_data_line(const OneDataLine& line);
	static bool is_header_line(const OneDataLine& line) {
		return line.size() > 0 && (line.front() == '>' || line.front() == '@');
	}
	bool is_nucl(const unsigned char c) const {
		return encode_table[c] < 16;
	}
	bool is_ambig_nucl(const unsigned char c) const {
		const int r(encode_table[c]);
		return 3 < r && r < 16;
	}
	static bool is_upper_case_letter(const char c) {
		return 'A' <= c && c <= 'Z';
	}
	static bool is_lower_case_letter(const char c) {
		return 'a' <= c && c <= 'z';
	}
	static bool is_alpha(const char c) {
		return is_upper_case_letter(c) || is_lower_case_letter(c);
	}
	static bool is_comment_line(const OneDataLine& line) {
		return line.front() == '#' || line.front() == '!';
	}
    private:
	BufferLineReader m_Reader;
	const u1_t* const encode_table;
};

#endif // FASTA_READER_H
