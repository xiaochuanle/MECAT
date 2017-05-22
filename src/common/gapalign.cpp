#include "gapalign.h"

#include <algorithm>

#include "defs.h"

using namespace std;

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
						int& tblk)
{
	bool last_block;
	int qleft = qsize - qidx;
	int tleft = tsize - tidx;
	if (qleft < desired_block_size + 100 || tleft < desired_block_size + 100) {
		qblk = min(qleft, static_cast<int>(tleft + tleft * 0.2));
		tblk = min(tleft, static_cast<int>(qleft + qleft * 0.2));
		last_block = true;
	} else {
		qblk = desired_block_size;
		tblk = desired_block_size;
		last_block = false;
	}
	
	if (forward) {
		Q = query + qidx;
		T = target + tidx;
	} else {
		Q = query - qidx;
		T = target - tidx;
	}
	
	return last_block;
}

bool
trim_mismatch_end(const char* qaln, 
				  const char* taln, 
				  const int aln_size, 
				  const int mat_cnt, 
				  int& qcnt, 
				  int& tcnt, 
				  int& aln_cnt)
{
	int m = 0, k;
	for (k = aln_size - 1, qcnt = 0, tcnt = 0, aln_cnt = 0; k >= 0 && m < mat_cnt; --k) {
		++aln_cnt;
		if (qaln[k] != GAP_CODE) ++qcnt;
		if (taln[k] != GAP_CODE) ++tcnt;
		if (qaln[k] == taln[k]) {
			++m;
		} else {
			m = 0;
		}
	}
	return m == mat_cnt && k > 0;
}
