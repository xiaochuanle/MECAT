#include "fasta_reader.h"

#include <algorithm>
#include <sstream>

idx_t FastaReader::read_one_seq(Sequence& seq)
{
    seq.clear();
    bool need_defline = true;
    while (++m_Reader)
    {
        const OneDataLine& odl = m_Reader.get_line();
        if (odl.size() == 0) continue;
        const int c = odl.front();
        if (c == '>' || c == '@')
        {
            if (need_defline)
            {
                x_parse_defline(odl, seq.header());
                need_defline = false;
                continue;
            }
            else
            {
                m_Reader.unget_line();
                break;
            }
        }
        else if (c == '+')
        {
			if (!++m_Reader)
				ERROR("FastaReader: quality score line is missing at around line %lld", (long long)m_Reader.line_number());
            break;
        }
        else if (is_comment_line(odl))
        {
            continue;
        }
        else if (need_defline)
        {
			ERROR("FastaReader: Input doesn't start with a defline or comment around line %lld", (long long)m_Reader.line_number());
        }

        x_parse_data_line(odl, seq.sequence());
    }

    if (seq.size() == 0 && seq.header().size() > 0)
    {
		ERROR("FastaReader: Near line %lld, sequence data is missing.", (long long)m_Reader.line_number());
    } 

    if (seq.header().size() == 0 && seq.sequence().size() == 0) return -1;
    return seq.size();
}


void FastaReader::x_parse_data_line(const OneDataLine& line, str_t& seq)
{
    x_check_data_line(line);
    const idx_t len = line.size();
    idx_t curr_pos = seq.size();
    seq.resize(seq.size() + len);
    idx_t pos = 0;
    for (pos = 0; pos < len; ++pos)
    {
        const int c = line[pos];
        if (c == ';') break;
        if (is_nucl(c) || c == '-')
        {
           seq[curr_pos++] = c;
        }
        else if (!isspace(c))
        {
			ERROR("FastaReader: There are invalid residue(s) around position %d of line %lld.", (int)(pos + 1), (long long)m_Reader.line_number());
        } 
    }

    seq.resize(curr_pos);
}


void FastaReader::x_parse_defline(const OneDataLine& line, str_t& header)
{
    header.clear();
    header.push_back(line.begin() + 1, line.size() - 1);
    if (header.size() == 0)
    {
        const idx_t line_number = m_Reader.line_number();
		ERROR("A sequence is given an empty header around line %lld.", (long long)line_number);
    }
}

void FastaReader::x_check_data_line(const OneDataLine& line)
{
    idx_t good = 0, bad = 0, len = line.size();
    idx_t ambig_nucl = 0;
    for (idx_t pos = 0; pos < len; ++pos)
    {
        const unsigned char c = line[pos];
        if (is_alpha(c) || c == '*')
        {
            ++good;
            if (is_nucl(c) && is_ambig_nucl(c)) ++ambig_nucl;
        }
        else if (c == '-')
        {
            ++good;
        }
        else if (isspace(c) || (c >= '0' && c <= '9'))
        {

        }
        else if (c == ';')
        {
            break;
        }
        else
        {
            ++bad;
        }
    }

    if (bad >= good / 3 && (len > 3 || good == 0 || bad > good))
    {
		ERROR("FastaReader: Near line %lld, there's a line that doesn't look like plausible data, but it's not marked as defline or commnet.",
			  (long long)m_Reader.line_number());
    }
}


