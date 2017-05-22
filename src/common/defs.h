#ifndef DEFS_H
#define DEFS_H

#include <cassert>
#include <cstdlib>
#include <stdint.h>

#include <iostream>
#include <sstream>

typedef int8_t      i1_t;
typedef uint8_t     u1_t;
typedef int16_t     i2_t;
typedef uint8_t     u2_t;
typedef int32_t     i4_t;
typedef uint32_t    u4_t;
typedef int64_t     i8_t;
typedef uint64_t    u8_t;
typedef i8_t        idx_t;

typedef u1_t		uint1;
typedef idx_t		index_t;

#define INVALID_IDX (-1)

#define LOG(out, ...) { \
    fprintf(out, "[%s, %u] ", __func__, __LINE__); \
    fprintf(out, __VA_ARGS__); \
    fprintf(out, "\n"); \
}

#define ERROR(...) { \
    LOG(stderr, __VA_ARGS__); \
    abort(); \
}

#define d_assert(r) assert(r)
#define RDBL 1
#if RDBL > 0
#define r_assert(r) \
do { \
    const bool __ra__rv__ = (r); \
    if (!__ra__rv__) \
        ERROR("assertion \'%s\' failed", #r); \
} while(0)
#else
#define r_assert(r) assert(r)
#endif

#define open_fstream(file, path, mode) \
do { \
    (file).open(path, mode); \
    if (!(file)) ERROR("failed to open file \'%s\' with mode \'%s\'", path, #mode); \
} while(0)

#define do_fstream(type, name, path, mode) \
do { \
	type name; \
	open_fstream(name, path, mode); \
} while(0)

#define close_fstream(file) (file).close()

#define sb_write(sb, buf, size) \
do { \
    const char* __sw__buf__ = (const char*)(buf); \
    std::streamsize __sw__ws__ = sb->sputn(__sw__buf__, size); \
    if (__sw__ws__ != (size)) ERROR("file write error"); \
} while(0)

#define sb_read(sb, buf,size) \
do { \
    char* __sr__buf__ = (char*)(buf); \
    std::streamsize __sr__rs__ = sb->sgetn(__sr__buf__, size); \
    if (__sr__rs__ != (size)) ERROR("file read error"); \
} while(0)

#ifndef ABS
#define ABS(a) ((a) >= 0 ? (a) : (-(a)))
#endif
#ifndef MAX
#define MAX(a, b) ((a) >= (b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a, b) ((a) <= (b) ? (a) : (b))
#endif

#define safe_malloc(arr, type, count) \
do { \
    size_t __sm__sz__ = sizeof(type) * (count); \
    (arr) = (type *)malloc(__sm__sz__); \
    if (!(arr)) ERROR("malloc fail"); \
} while(0)

#define safe_calloc(arr, type, count) \
do { \
    size_t __sc__sz__ = sizeof(type) * (count); \
    (arr) = (type *)calloc(1, __sc__sz__); \
    if (!(arr)) ERROR("calloc fail"); \
} while(0)

#define safe_realloc(arr, type, count) \
do { \
	size_t __sr__size__ = sizeof(type) * count; \
	arr = (type *)realloc(arr, __sr__size__); \
	if (!arr) \
	{ \
		ERROR("failed to realloc memory."); \
		exit(1); \
	} \
} while(0)

#include <new>
#define snew(arr, type, count) \
	do { \
    size_t snew_cnt = static_cast<size_t>(count); \
		try { \
        (arr) = new type[snew_cnt]; \
		} catch (std::bad_alloc& snew_e) { \
        ERROR("%s", snew_e.what()); \
    } \
} while(0)

#define sznew(arr, type, count) \
		do { \
		size_t zsnew_cnt = static_cast<size_t>(count); \
		snew(arr, type, count); \
		size_t zsnew_s = zsnew_cnt * sizeof(type); \
		memset(static_cast<void*>(arr), 0, zsnew_s); \
	} while(0)

#define dsnew(arr, type, count) type * arr; snew(arr, type, count)

#define dsznew(arr, type, count) type * arr; sznew(arr, type, count)

#define sfree(arr) delete[] arr

#define SAFE_WRITE(p, type, count, out) \
	do { \
	size_t __sw__es__ = sizeof(type); \
	size_t __sw__wn__ = fwrite(p, __sw__es__, (size_t)count, out); \
    if (__sw__wn__ != (size_t)count) \
	{ \
		ERROR("write error!"); \
		exit(1); \
	} \
} while (0)

#define SAFE_READ(p, type, count, in) \
		do { \
		size_t __sr__es__ = sizeof(type); \
		size_t __sr__rn__ = fread(p, __sr__es__, (size_t)count, in); \
		if (__sr__rn__ != (size_t)count) \
		{ \
			ERROR("read error!"); \
			exit(1); \
		} \
	} while (0)


#define safe_free(arr) free(arr)

#include <sys/time.h>

struct Timer
{
    struct timeval start;
    struct timeval end;

    void go() { gettimeofday(&start, NULL); }
    void stop() { gettimeofday(&end, NULL); }
    double elapsed() { return end.tv_sec - start.tv_sec + 1.0 * (end.tv_usec - start.tv_usec) / 1000000; }
};

struct DynamicTimer
{
	DynamicTimer(const char* func) : m_func(func) 
	{ 
		if (m_func) fprintf(stderr, "[%s] begins.\n", m_func);
		timer.go(); 
	}
	~DynamicTimer() 
	{
		timer.stop();
		fprintf(stderr, "[%s] takes %.2f secs.\n", m_func, timer.elapsed());
	}
	
private:
	const char* m_func;
	Timer timer;
};

const u1_t* get_dna_encode_table();
const char* get_dna_decode_table();
const u1_t* get_dna_complement_table();
#define GAP '-'
#define GAP_CODE 4
#define GAP_CHAR '-'
#define FWD 0
#define REV 1
#define REVERSE_STRAND(s) (1-(s))
#define MAX_SEQ_SIZE 500000
#define MAX_INVALID_END_SIZE 200
#define MIN_EXTEND_SIZE 500
#define MIN_OVERLAP_SIZE 1000
#define TECH_PACBIO 0
#define TECH_NANOPORE 1

#endif //  DEFS_H
