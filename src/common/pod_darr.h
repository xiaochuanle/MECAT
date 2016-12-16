#ifndef DARR_H
#define DARR_H

#include <cassert>
#include <cstring>

#include <iterator>
#include <iostream>
#include <sstream>

#include "defs.h"

template <typename T>
class PODArray
{
public:
    typedef typename std::iterator<std::random_access_iterator_tag, T>            short_for_iterator;
    typedef typename std::iterator_traits<short_for_iterator>::iterator_category  iterator_category;
    typedef typename std::iterator_traits<short_for_iterator>::value_type         value_type;
    typedef typename std::iterator_traits<short_for_iterator>::pointer            pointer;
    typedef typename std::iterator_traits<short_for_iterator>::reference          reference;
    typedef typename std::iterator_traits<short_for_iterator>::difference_type    difference_type;

public:
    PODArray(const idx_t alloc_size = 1024)
    {
        used_size_ = 0;
        alloc_size_ = std::max((idx_t)1024, alloc_size);
        safe_malloc(data_, T, alloc_size_);
    }
    void destroy()
    {
        if (alloc_size_)
        {
            safe_free(data_);
        }
        used_size_ = alloc_size_ = 0;
    }
    ~PODArray()
    {
        destroy();
    }

public:
    pointer begin() { return data_; }
    const pointer begin () const { return data_; }
    pointer end() { return data_ + used_size_; }
    const pointer end() const { return data_ + used_size_; }
    reference front() { return *data_; }
    const reference front() const { return *data_; }
    reference back() { return data_[used_size_ - 1]; }
    const reference back() const { return data_[used_size_ - 1]; }
    reference operator[](const idx_t idx) { return data_[idx]; }
    const reference operator[](const idx_t idx) const { return data_[idx]; }
    idx_t size() { return used_size_; }
    idx_t size() const { return used_size_; }
    pointer data() { return data_; }
    const pointer data() const { return data_; }
    void clear() { used_size_ = 0; }

    void push_back(const T* p, const idx_t size)
    {
        if (used_size_ + size > alloc_size_)
        {
            idx_t new_alloc_size = alloc_size_ ? alloc_size_ : 1024;
            while (used_size_ + size > new_alloc_size) new_alloc_size *= 2;
            T* new_data = NULL;
            safe_malloc(new_data, T, new_alloc_size);
            memcpy(new_data, data_, sizeof(T) * used_size_);
            safe_free(data_);
            data_ = new_data;
            alloc_size_ = new_alloc_size;
        }
        memcpy(data_ + used_size_, p, sizeof(T) * size);
        used_size_ += size;
    }
	
	void push_back(const T& t)
	{
		push_back(&t, 1);
	}

    void pop_back()
    {
        --used_size_;
    }

    void resize(const idx_t new_size)
    {
        if (new_size > alloc_size_)
        {
            idx_t new_alloc_size = alloc_size_ ? alloc_size_ : 1024;
            while (new_size > new_alloc_size) new_alloc_size *= 2;
            T* new_data = NULL;
            safe_malloc(new_data, T, new_alloc_size);
            memcpy(new_data, data_, sizeof(T) * used_size_);
            safe_free(data_);
            data_ = new_data;
            alloc_size_ = new_alloc_size;
        }

        used_size_ = new_size;
    }
	
	void expand()
	{
		resize(used_size_ + 1);
	}

    void reserve(const idx_t size)
    {
        if (alloc_size_ < size)
        {
            idx_t new_alloc_size = alloc_size_ ? alloc_size_ : 1024;
            while (new_alloc_size < size) new_alloc_size *= 2;
            T* new_data = NULL;
            safe_malloc(new_data, T, new_alloc_size);
            memcpy(new_data, data_, sizeof(T) * used_size_);
            safe_free(data_);
            data_ = new_data;
            alloc_size_ = new_alloc_size;
        }
    }

private:
    idx_t     used_size_;
    idx_t     alloc_size_;
    T*          data_;
};

#endif // DARR_H
