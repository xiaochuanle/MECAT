#include "sequence.h"

std::ostream& operator<<(std::ostream& out, const Sequence& seq)
{
	out << ">";
	const Sequence::str_t& header = seq.header();
	idx_t n = header.size();
	for (idx_t i = 0; i < n; ++i) out << header[i];
	out << "\n";
	
	const Sequence::str_t& s = seq.sequence();
	n = s.size();
	for (idx_t i = 0; i < n; ++i) out << s[i];
	out << "\n";
	
	return out;
}
