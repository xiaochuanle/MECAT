#include "xdrop_gapalign.h"
//#include "smart_assert.h"

#include <algorithm>

using namespace std;

#define MIN_SCORE (-100000000)

int
xdrop_align(const char* A, 
            const int M,
            const char* B,
            const int N,
            int matrix[][4],
            int gap_open,
            int gap_extend,
            int x_dropoff,
            u8* state_array,
            BlastGapDP* score_array,
            u8** edit_script,
            int* edit_start_offset,
            GapPrelimEditBlock* edit_block,
            const bool forward,
            int& ae,
            int& be)
{
    ae = 0;
    be = 0;
    edit_block->clear();
    if (M <= 0 || N <= 0) return 0;

    int i;
    int a_index;
    int b_index, b_size, first_b_index, last_b_index;
    int gap_open_extend = gap_open + gap_extend;
    int best_score;
    int* matrix_row = NULL;
    int score;
    int score_gap_row;
    int score_gap_col;
    int next_score;
    u8* edit_script_row = NULL;
    int orig_b_index;
    u8 script, next_script, script_row, script_col;
    int states_used = 0;

    if (x_dropoff < gap_open_extend) x_dropoff = gap_open_extend;
    edit_script[0] = state_array;
    edit_start_offset[0] = 0;
    edit_script_row = state_array;

    score = -gap_open_extend;
    score_array[0].best = 0;
    score_array[0].best_gap = -gap_open_extend;
    for (i = 1; i <= N; ++i) {
        if (score < -x_dropoff) break;
        score_array[i].best = score;
        score_array[i].best_gap = score - gap_open_extend;
        score -= gap_extend;
        edit_script_row[i] = SCRIPT_GAP_IN_A;
    }
    states_used = min(N, i + 1);

    b_size = i;
    best_score = 0;
    first_b_index = 0;

    for (a_index = 1; a_index <= M; ++a_index) {
        const u8 AC = extract_char<u8>(A, a_index - 1, forward);
        edit_script[a_index] = state_array + states_used + 1;
        edit_start_offset[a_index] = first_b_index;
        
        /* The inner loop assumes the current traceback
         * row begins at offset zero of B */
        edit_script_row = edit_script[a_index] - first_b_index;
        orig_b_index = first_b_index;
        matrix_row = matrix[ AC ];

        score = MIN_SCORE;
        score_gap_row = MIN_SCORE;
        last_b_index = first_b_index;

        for (b_index = first_b_index; b_index < b_size; ++b_index) {
            const u8 BC = extract_char<u8>(B, b_index, forward);
            score_gap_col = score_array[b_index].best_gap;
            next_score = score_array[b_index].best + matrix_row[ BC ];

            /* script, script_row and script_col contain
             * the actions specified by the dynmaic programming.
             * when the inner loop has finished, 'script' contains
             * all of the actions to perform, and is
             * written to edit_script[a_index][b_index]. */

            script = SCRIPT_SUB;
            script_row = SCRIPT_EXTEND_GAP_B;
            script_col = SCRIPT_EXTEND_GAP_A;

            if (score < score_gap_col) {
                script = SCRIPT_GAP_IN_B;
                score = score_gap_col;
            }

            if (score < score_gap_row) {
                script = SCRIPT_GAP_IN_A;
                score = score_gap_row;
            }

            if (best_score - score > x_dropoff) {
                if (first_b_index == b_index) ++first_b_index;
                else score_array[b_index].best = MIN_SCORE;
            } else {
                last_b_index = b_index;
                if (score > best_score) {
                    best_score = score;
                    ae = a_index;
                    be = b_index;
                }

                score_gap_col -= gap_extend;
                if (score_gap_col < (score - gap_open_extend)) {
                    score_array[b_index].best_gap = score - gap_open_extend;
                } else {
                    score_array[b_index].best_gap = score_gap_col;
                    script += script_col;
                }

                score_gap_row -= gap_extend;
                if (score_gap_row < (score - gap_open_extend)) {
                    score_gap_row = score - gap_open_extend;
                } else {
                    script += script_row;
                }

                score_array[b_index].best = score;
            }

            score = next_score;
            edit_script_row[b_index] = script;
        }

        if (first_b_index == b_size) break;

        if (last_b_index < b_size - 1) {
            b_size = last_b_index + 1;
        } else {
            while (score_gap_row >= (best_score - x_dropoff) && b_size < N) {
                score_array[b_size].best = score_gap_row;
                score_array[b_size].best_gap = score_gap_row - gap_open_extend;
                score_gap_row -= gap_extend;
                edit_script_row[b_size] = SCRIPT_GAP_IN_A;
                ++b_size;
            }
        }

        /* update memory allocator to reflect the exact
         * number of traceback cells this row needed. */
        states_used += std::max(b_index, b_size) - orig_b_index + 1;

        if (b_size < N) {
            score_array[b_size].best = MIN_SCORE;
            score_array[b_size].best_gap = MIN_SCORE;
            ++b_size;
        }
    }

    /* pick the optimal path through the now complete edit_script[][].
     * This is equivalent to flattening the 2-D array into a 1-D list of actions. */

    a_index = ae;
    b_index = be;
    script = SCRIPT_SUB;

    while (a_index > 0 || b_index > 0) {
        /* Retrieve the next action to perform.
         * Rows of the traceback array do not necessarily start
         * at the offset zero of B, so a correction is needed
         * to point to the correct position. */

        next_script = edit_script[a_index][b_index - edit_start_offset[a_index]];

        switch (script) {
            case SCRIPT_GAP_IN_A:
                script = next_script & SCRIPT_OP_MASK;
                if (next_script & SCRIPT_EXTEND_GAP_A)
                    script = SCRIPT_GAP_IN_A;
                break;

            case SCRIPT_GAP_IN_B:
                script = next_script & SCRIPT_OP_MASK;
                if (next_script & SCRIPT_EXTEND_GAP_B)
                    script = SCRIPT_GAP_IN_B;
                break;

            default:
                script = next_script & SCRIPT_OP_MASK;
                break;
        }

        if (script == SCRIPT_GAP_IN_A) {
            --b_index;
        } else if (script == SCRIPT_GAP_IN_B) {
            --a_index;
        } else {
            --a_index;
            --b_index;
        }

        edit_block->add((EGapAlignOpType)script, 1);
    }

    return best_score;
}

void
script_to_aligned_string(const char* query, 
						 const char* target,
						 int qe,
						 int te,
						 GapPrelimEditBlock* edit_block,
						 const bool forward,
						 string& qaln,
						 string& taln)
{
	qaln.clear();
	taln.clear();
	const char* q = query;
	const char* t = target;
	int inc = forward ? 1 : -1;
	for (int i = edit_block->num_ops - 1; i >= 0; --i) {
		EGapAlignOpType op_type = edit_block->edit_ops[i].op_type;
		int n = edit_block->edit_ops[i].num;
		switch (op_type) {
			case eGapAlignSub:
				for (int j = 0; j < n; ++j) {
					qaln += (*q);
					taln += (*t);
					q += inc;
					t += inc;
				}
				break;
			case SCRIPT_GAP_IN_A:
				for (int j = 0; j < n; ++j) {
					qaln += GAP_CODE;
					taln += (*t);
					t += inc;
				}
				break;
			case SCRIPT_GAP_IN_B:
				for (int j = 0; j < n; ++j) {
					qaln += (*q);
					taln += GAP_CODE;
					q += inc;
				}
				break;
			default:
				ERROR("i = %d, unknown op_type: %d", i, (int)op_type);
				break;
		}
	}
}

void 
align_ex(const char* query, 
		 const int query_size,
		 const char* target, 
		 const int target_size,
		 char* qbf, 
		 char* tbf,
		 XdropAlignParameters* xap,
		 int matrix[][4], 
		 u8* state_array, 
		 BlastGapDP* score_array,
		 u8** edit_script, 
		 int* edit_start_offset, 
		 GapPrelimEditBlock* edit_block,
		 string& qaln,
		 string& taln,
		 bool forward)
{
	const int kTailMatchBP = 4;
	const int kBlkSize = xap->block_size;
	int qidx = 0, tidx = 0;
	int qblk, tblk;
	int aln_qe, aln_te;
	qaln.clear();
	taln.clear();
	string tmp_qaln, tmp_taln;
	const char* Q;
	const char* T;
	while (1) {
		bool last_block = retrieve_next_aln_block(query,
												  qidx,
												  query_size,
												  target,
												  tidx,
												  target_size,
												  kBlkSize,
												  forward,
												  Q,
												  T,
												  qblk,
												  tblk);
		
		int score = xdrop_align(Q, 
								qblk, 
								T, 
								tblk, 
								matrix, 
								xap->gap_open, 
								xap->gap_extend, 
								xap->x_dropoff, 
								state_array, 
								score_array, 
								edit_script, 
								edit_start_offset, 
								edit_block, 
								forward,
							    aln_qe,
							    aln_te); 
		
		//cout << "qidx = " << qidx << ", qe = " << aln_qe << ", tidx = " << tidx << ", te = " << aln_te 
		//	 << ", score = " << score << endl;

		script_to_aligned_string(Q,
								 T, 
								 aln_qe, 
								 aln_te, 
								 edit_block, 
								 forward, 
								 tmp_qaln, 
								 tmp_taln);
		
		bool full_map = false;
		if (qblk - aln_qe <= 20 || tblk - aln_te <= 20) full_map = true;
		if ((!full_map) || last_block) {
			qaln += tmp_qaln;
			taln += tmp_taln;
			break;
		}
		int qcnt = 0, tcnt = 0, acnt = 0;
		const bool trim = trim_mismatch_end(tmp_qaln.data(), 
											tmp_taln.data(), 
											tmp_qaln.size(), 
											kTailMatchBP, 
											qcnt, 
											tcnt, 
											acnt);
		if (!trim) break;
		tmp_qaln.resize(tmp_qaln.size() - acnt);
		tmp_taln.resize(tmp_taln.size() - acnt);
		qaln += tmp_qaln;
		taln += tmp_taln;
		qidx += (aln_qe - qcnt);
		tidx += (aln_te - tcnt);
	}
}
		 
bool
XdropAligner::go(const char* query, const int query_start, const int query_size,
				 const char* target, const int target_start, const int target_size,
				 const int min_aln_size)
{	
	align_ex(query + query_start - 1, 
			 query_start,
			 target + target_start - 1,
			 target_start,
			 qbuf,
			 tbuf,
			 &param,
			 score_matrix,
			 state_array,
			 score_array,
			 edit_script,
			 edit_start_offset,
			 &edit_block,
			 left_qaln,
			 left_taln,
			 false);
	
	align_ex(query + query_start, 
			 query_size - query_start,
			 target + target_start,
			 target_size - target_start,
			 qbuf,
			 tbuf,
			 &param,
			 score_matrix,
			 state_array,
			 score_array,
			 edit_script,
			 edit_start_offset,
			 &edit_block,
			 right_qaln,
			 right_taln,
			 true);
	
	int aln_idx = 0;
	int i, j, k, n;
	const char* dt = "ACGT-";
	n = left_qaln.size() - 1;
	for (i = 0, j = 0, k = n - 1; k >= 0; --k, ++aln_idx) {
		unsigned char c = left_qaln[k];
		r_assert(c >= 0 && c <= 4);
		if (c != GAP_CODE) ++i;
		qaln[aln_idx] = dt[c];
		
		c = left_taln[k];
		r_assert(c >= 0 && c <= 4);
		if (c != GAP_CODE) ++j;
		taln[aln_idx] = dt[c];
	}
	qoff = query_start - i;
	toff = target_start - j;
	r_assert(qoff >= 0);
	r_assert(toff >= 0);
	
	n = right_qaln.size();
	for (i = 0, j = 0, k = 0; k < n; ++k, ++aln_idx) {
		unsigned char c = right_qaln[k];
		r_assert(c >= 0 && c <= 4);
		if (c != GAP_CODE) ++i;
		qaln[aln_idx] = dt[c];
		
		c = right_taln[k];
		r_assert(c >= 0 && c <= 4);
		if (c != GAP_CODE) ++j;
		taln[aln_idx] = dt[c];
	}
	qaln[aln_idx] = '\0';
	taln[aln_idx] = '\0';
	aln_size = aln_idx;
	qend = query_start + i;
	tend = target_start + j;
	r_assert(qend <= query_size);
	r_assert(tend <= target_size);

	return qend - qoff >= min_aln_size;
}
