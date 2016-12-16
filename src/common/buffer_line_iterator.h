#ifndef BUFFER_LINE_ITERATOR_H
#define BUFFER_LINE_ITERATOR_H

#include <cstring>
#include <stdint.h>

#include <iostream> // std::ios std::istream
#include <fstream>  // std::filebuf
#include <string>
#include <vector>

#include "defs.h"
#include "pod_darr.h"

class BufferLineReader
{
public:
    typedef PODArray<char>  OneDataLine;

public:
    BufferLineReader(const char* file_name);
    ~BufferLineReader();
    OneDataLine& get_line()
    { return line_; }
    bool eof() const;
    bool operator++();
    void unget_line()
    {
        --line_number_;
        unget_line_ = true;
    }
    idx_t line_number() { return line_number_; }
    idx_t line_number() const { return line_number_; }

private:
    void x_load_long();
    bool x_read_buffer();
    
private:
    std::filebuf    fb_;
    std::istream*   ins_;
    static const idx_t kBufferSize = 1024 * 1024 * 8;
    char*           buf_;
    idx_t         cur_;
    idx_t         buf_sz_;
    bool            done_; 
    OneDataLine     line_;
    bool            unget_line_;
    idx_t         line_number_;
};

#endif // BUFFER_LINE_ITERATOR_H
