#ifndef XDROP_GAPALIGN_H
#define XDROP_GAPALIGN_H

#include "defs.h"
#include "gapalign.h"
//#include "smart_assert.h"

#include <string>

typedef uint8_t u8;

/* operation types within the edit script */
enum EGapAlignOpType {
    eGapAlignDel = 0, /* deletion: a gap in query */
    eGapAlignDel2 = 1, /* Frame shift deletion of two nucleotides */
    eGapAlignDel1 = 2, /* Frame shift deletion of one nucleotide */
    eGapAlignSub = 3, /* substitution */
    eGapAlignIns1 = 4, /* Frame shift insertion of one nucleotide */
    eGapAlignIns2 = 5, /* Frame shift insertion of two nucleotides */
    eGapAlignIns = 6, /* insertion: a gap in subject */
    eGapAlignDecline = 7, /* Non-aligned region */
    eGapAlignInvalid = 8 /* invalid operation */
};

// values for the editing script operation in trace back
enum {
    SCRIPT_SUB          = eGapAlignSub, /* substitution */
    SCRIPT_GAP_IN_A     = eGapAlignDel, /* deletion */
    SCRIPT_GAP_IN_B     = eGapAlignIns, /* insertion */
    SCRIPT_OP_MASK      = 0x07, /* mask for edit script operation */

    SCRIPT_EXTEND_GAP_A = 0x10, /* continue a gap in A */
    SCRIPT_EXTEND_GAP_B = 0x40 /* continue a gap in B */
};

/* auxiliary structure for dynamic programming gapped extension */
struct BlastGapDP {
    int best; /* score of best path that ends in a match at this position */
    int best_gap; /* score of best path that ends in a gap at this position */
};

struct GapPrelimEditScript
{
    EGapAlignOpType op_type;
    int num;
};

struct GapPrelimEditBlock
{
    GapPrelimEditScript* edit_ops;
    int num_ops_allocated;
    int num_ops;
    EGapAlignOpType last_op;

    GapPrelimEditBlock() {
        num_ops_allocated = MAX_SEQ_SIZE;
        snew(edit_ops, GapPrelimEditScript, num_ops_allocated);
        num_ops = 0;
        last_op = eGapAlignInvalid;
    }

    ~GapPrelimEditBlock() {
        sfree(edit_ops);
    }

    void clear() {
        last_op = eGapAlignInvalid;
        num_ops = 0;
    }

    void add(EGapAlignOpType op_type, int n) {
        if (n == 0) return;
        if (last_op == op_type) {
            edit_ops[num_ops - 1].num += n;
        }
        else {
            last_op = op_type;
            edit_ops[num_ops].op_type = op_type;
            edit_ops[num_ops].num = n;
            ++num_ops;
        }
    }
};

struct XdropAlignParameters
{
	int reward;
	int penalty;
	int gap_open;
	int gap_extend;
	int x_dropoff;
	
	int block_size;
	int state_array_size;
	int score_array_size;
	int edit_script_row_size;
	
	void init(const int large_block) {
		reward = 1;
		penalty = -1;
		gap_open = 0;
		gap_extend = 1;
		x_dropoff = 30;
		
		if (large_block) {
			block_size = 1000;
		} else {
			block_size = 500;
		}
		
		state_array_size = 5000000;
		score_array_size = 4096;
		edit_script_row_size = 4096;
	}
};

class XdropAligner : public GapAligner
{
public:
	XdropAligner(const int large_block = 0) {
			param.init(large_block);
			snew(state_array, u8, param.state_array_size);
			snew(score_array, BlastGapDP, param.score_array_size);
			snew(edit_script, u8*, param.edit_script_row_size);
			snew(edit_start_offset, int, param.edit_script_row_size);
			build_score_matrix();
			snew(qbuf, char, param.block_size * 2);
			snew(tbuf, char, param.block_size * 2);
			snew(qaln, char, MAX_SEQ_SIZE);
			snew(taln, char, MAX_SEQ_SIZE);
		}
		
	virtual ~XdropAligner() {
		sfree(state_array);
		sfree(score_array);
		sfree(edit_script);
		sfree(edit_start_offset);
		sfree(qbuf);
		sfree(tbuf);
		sfree(qaln);
		sfree(taln);
	}
	
	void build_score_matrix() {
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				if (i == j) score_matrix[i][j] = param.reward;
				else score_matrix[i][j] = param.penalty;
			}
		}
	}
	
	virtual double calc_ident() const {
		if (aln_size == 0) return 0;
		
		int n = aln_size, ident = 0;
		for (int i = 0; i < n; ++i) {
			if (qaln[i] == taln[i]) ++ident;
		}
		return 100.0 * ident / n;
	}
	
	virtual bool go(const char* query, const int query_start, const int query_size,
			const char* target, const int target_start, const int target_size,
			const int min_aln_size);
	
	virtual int query_start() const {
		return qoff;
	}
	
	int query_end() const {
		return qend;
	}
	
	virtual int target_start() const {
		return toff;
	}
	
	virtual int target_end() const {
		return tend;
	}
	
	virtual char* query_mapped_string() {
		return qaln;
	}
	
	virtual char* target_mapped_string() {
		return taln;
	}
		
public:
	XdropAlignParameters	param;
	u8*						state_array;
	BlastGapDP*				score_array;
	u8**					edit_script;
	int*					edit_start_offset;
	GapPrelimEditBlock		edit_block;
	int						score_matrix[4][4];
	std::string				left_qaln;
	std::string				left_taln;
	std::string				right_qaln;
	std::string				right_taln;
	char*					qaln;
	char*					taln;
	int						aln_size;
	int						qoff;
	int						qend;
	int						toff;
	int						tend;
	char*					qbuf;
	char*					tbuf;
};

#endif // XDROP_GAPALIGN_H
