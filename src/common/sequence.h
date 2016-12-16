#ifndef SEQUENCE_H
#define SEQUENCE_H

#include "buffer_line_iterator.h"

#include <iostream>

class Sequence
{
public:
    typedef BufferLineReader LineReader;
    typedef BufferLineReader::OneDataLine str_t;

public:
    str_t& header() { return header_; }
    const str_t& header() const { return header_; }
    str_t& sequence() { return sequence_; }
    const str_t& sequence() const { return sequence_; }
    idx_t size() { return sequence_.size(); }
    idx_t size() const { return sequence_.size(); }
    void clear() { header_.clear(); sequence_.clear(); }

    idx_t read_one_seq(LineReader& line_reader);

private:
    str_t header_;
    str_t sequence_;
};

std::ostream& operator<<(std::ostream& out, const Sequence& seq);

#endif // SEQEUNCE_H
