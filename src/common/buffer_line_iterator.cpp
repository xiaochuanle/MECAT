#include "buffer_line_iterator.h"

BufferLineReader::BufferLineReader(const char* file_name)
{
    if (!fb_.open(file_name, std::ios::in)) ERROR("cannot open file \'%s\' for reading", file_name);

    ins_ = new std::istream(&fb_);
    buf_ = new char[kBufferSize];
    done_ = false;
    unget_line_ = false;
    line_number_ = 0;
    x_read_buffer();
}


bool BufferLineReader::eof() const
{
     return done_ && (cur_ >= buf_sz_);
}


bool BufferLineReader::operator++()
{
    ++line_number_;

    if (unget_line_) 
    {
        unget_line_ = false;
        return true;
    }    

#define buffer_read_ret (!(done_ && line_.size() == 0))

    line_.clear();
    const idx_t start = cur_;
    const idx_t end = buf_sz_;
    for (idx_t p = start; p < end; ++p)
    {
        const int c = buf_[p];
        if (c == '\n')
        {
            line_.push_back(buf_ + start, p - start);
            cur_ = ++p;
            if (p == end)
                x_read_buffer();
            return buffer_read_ret;
        }
        else if (c == '\r')
        {
            line_.push_back(buf_ + start, p - start);
            cur_ = ++p;
            if (p == end)
            {
                if (x_read_buffer())
                {
                    p = cur_;
                    if (buf_[p] == '\n')            
                        cur_ = p + 1;
                }
                return buffer_read_ret;
            }
            if (buf_[p] != '\n') return buffer_read_ret;
            cur_ = ++p;
            if (p == end)
                x_read_buffer();
            return buffer_read_ret;
        }
    }

    x_load_long();
    return buffer_read_ret;
}


void BufferLineReader::x_load_long()
{
    idx_t start = cur_;
    idx_t end = buf_sz_;
    line_.push_back(buf_ + start, end - start);
    while (x_read_buffer())
    {
        start = cur_;
        end = buf_sz_;
        for (idx_t p = start; p < end; ++p)
        {
            const int c = buf_[p];
            if (c == '\r' || c == '\n')
            {
                line_.push_back(buf_ + start, p - start);
                if (++p == end)
                {
                    if (x_read_buffer())
                    {
                        p = cur_;
                        end = buf_sz_;
                        if (p < end && c == '\r' && buf_[p] == '\n') { ++p; cur_ = p; }
                    }
                }
                else
                {
                    if (c == '\r' && buf_[p] == '\n')
                    {
                        if (++p == end) 
                        {
                            x_read_buffer();
                            p = cur_;
                        }
                    }
                    cur_ = p;
                }
                return;
            }
        }
        line_.push_back(buf_ + start, end - start);
    }
}


bool BufferLineReader::x_read_buffer()
{
    std::streambuf* sb = ins_->rdbuf();
    bool ok = sb && ins_->good();
    std::streamsize r = ok ? sb->sgetn(buf_, kBufferSize) : 0;
    buf_sz_ = (idx_t)r;
    cur_ = 0;
    if (buf_sz_ == 0)
    {
        done_ = true;
        return false;
    }
    else if (buf_sz_ < kBufferSize)
        done_ = true;
    return true;
}

BufferLineReader::~BufferLineReader()
{
    delete ins_;
    delete[] buf_;
}

