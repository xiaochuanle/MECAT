#ifndef DW_H
#define DW_H

#include <math.h>	// ceil()
#include <stdint.h>	// int64_t, uint8_t
#include <string>	// string
#include <vector>	// vector<>

class Alignment {
    public:
	// current_size tracks actual buffer use
	std::vector<uint8_t> q_aln_str, t_aln_str;
	int aln_q_e, aln_t_e, current_size;
    public:
	explicit Alignment(const size_t max_size) : q_aln_str(max_size), t_aln_str(max_size) { }
	~Alignment() { }
	size_t size() const {
		return q_aln_str.size();
	}
	void clear() {
		current_size = 0;
	}
};

class OutputStore {
    public:
	// these track actual buffer use - buffer use starts at
	// buffer_start, with left going down, and right going up
	std::vector<uint8_t> q_buffer, t_buffer;
	int buffer_start, left_size, right_size;
	int query_start, query_end;
	int target_start, target_end;
    public:
	explicit OutputStore(const size_t max_size) : q_buffer(max_size), t_buffer(max_size) { }
	~OutputStore() { }
	size_t size() const {
		return q_buffer.size();
	}
	void clear(const int i) {
		buffer_start = i;
		left_size = right_size = 0;
	}
};

struct DPathData {
	// using shorts to reduce size of largest buffer
	short int x1, y1, x2, y2, pre_k;
	void set(const int i, const int j, const int k, const int l, const int m) {
		x1 = i;
		y1 = j;
		x2 = k;
		y2 = l;
		pre_k = m;
	}
};

struct DPathIndex {
	int d_offset, min_k;
	void set(const int i, const int j) {
		d_offset = i;
		min_k = j;
	}
};

struct PathPoint {
	int x, y;
	void set(const int i, const int j) {
		x = i;
		y = j;
	}
};

// if the end of query/target is within this or less, extend the whole distance
#define SEGMENT_BORDER 100

// use static sizing here rather than dynamic, as dynamic is 50% slower in practice;
// that said, reserving the space instead of initializing it may be viable

class DiffRunningData {
    public:
	static const int segment_size = 500;		// 500 is "small", 1000 is "large"
	std::vector<int> DynQ, DynT;
	std::vector<DPathIndex> d_path_index;
	std::vector<PathPoint> aln_path;
	std::vector<DPathData> d_path;
	Alignment align;
	OutputStore result;
    public:
	explicit DiffRunningData(const double error_rate, const int64_t max_read_size) :
		// if the definitions of k_offset, band_tolerance, max_band_size
		// in dw.cpp are changed, you'll need to update these to reflect them
		DynQ(ceil((segment_size + SEGMENT_BORDER) * 4 * error_rate) * 2),
		DynT(ceil((segment_size + SEGMENT_BORDER) * 4 * error_rate) * 2),
		// effectively a right triangle on a rectangle,
		// as it's bounded geometric growth (4 * error_rate limited to .3)
		d_path_index(ceil((segment_size + SEGMENT_BORDER) * 4 * error_rate)),
		aln_path(ceil((segment_size + SEGMENT_BORDER) * 4 * error_rate) * 4),
		d_path(4 * error_rate < .3 ? ceil((segment_size + SEGMENT_BORDER) * .3 + 1) * ceil((segment_size + SEGMENT_BORDER) * .3 + 2) / 2 + ceil((segment_size + SEGMENT_BORDER) * .3 + 1) * ceil((segment_size + SEGMENT_BORDER) * (4 * error_rate - .3)) : ceil((segment_size + SEGMENT_BORDER) * 4 * error_rate) * ceil((segment_size + SEGMENT_BORDER) * 4 * error_rate + 1) / 2),
		align((segment_size + SEGMENT_BORDER) * 2),
		result(max_read_size * 2) { }
	~DiffRunningData() { }
};

class M5Record {
    public:
	int64_t qoff, qend, soff, send;
	std::string qaln, saln;
	explicit M5Record() { }
	~M5Record() { }
};

int GetAlignment(const std::string& query, int query_start, const std::string& target, int target_start, DiffRunningData& drd, M5Record& m5, double error_rate, int min_aln_size);

#endif  // DW_H
