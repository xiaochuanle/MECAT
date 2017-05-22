#ifndef GAP_ALIGN_H
#define GAP_ALIGN_H

class GapAligner
{
public:
    GapAligner() {

    }
    virtual ~GapAligner() {

    }

    virtual bool go(const char* query, const int query_start, const int query_size,
                    const char* target, const int target_start, const int target_size,
                    const int min_aln_size) = 0;

    virtual int query_start() const = 0;
    virtual int query_end() const = 0;
    virtual int target_start() const = 0;
    virtual int target_end() const = 0;
    virtual double calc_ident() const = 0;
	virtual char* query_mapped_string() = 0;
	virtual char* target_mapped_string() = 0;
};

template <typename T>
inline T extract_char(const char* A, int i, bool forward)
{
    if (forward) {
        return static_cast<T>(A[i]);
    } else {
        return static_cast<T>(A[-i]);
    }
}

bool
retrieve_next_aln_block(const char* query,
						int qidx, 
						const int qsize, 
						const char* target,
						int tidx, 
						const int tsize, 
						const int desired_block_size, 
						const bool forward,
						const char*& Q,
						const char*& T,
						int& qblk, 
						int& tblk);

bool
trim_mismatch_end(const char* qaln, 
				  const char* taln, 
				  const int aln_size, 
				  const int mat_cnt, 
				  int& qcnt, 
				  int& tcnt, 
				  int& aln_cnt);

#endif // GAP_ALIGN_H
