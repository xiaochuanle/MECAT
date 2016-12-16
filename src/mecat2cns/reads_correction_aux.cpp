#include "reads_correction_aux.h"

void normalize_gaps(const char* qstr, const char* tstr, const index_t aln_size, std::string& qnorm, std::string& tnorm, const bool push)
{
    qnorm.clear();
    tnorm.clear();
    const char kGap = '-';

#ifndef NDEBUG
    int qcnt = 0, tcnt = 0;
    for (index_t i = 0; i < aln_size; ++i)
    {
        const char qc = qstr[i];
        const char tc = tstr[i];
        if (qc != kGap) ++qcnt;
        if (tc != kGap) ++tcnt;
    }
#endif

    // convert mismatches to indels
    for (index_t i = 0; i < aln_size; ++i)
    {
        const char qc = qstr[i];
        const char tc = tstr[i];
        if (qc != tc && qc != kGap && tc != kGap)
        { qnorm += kGap; qnorm += qc; tnorm += tc; tnorm += kGap; }
        else
        { qnorm += qc; tnorm += tc; }
    }

    // push gaps to the right, but not pass the end
    if (push)
    {
        index_t qlen = qnorm.size();
        index_t tlen = tnorm.size();
        for (index_t i = 0; i < qlen - 1; ++i)
        {
            // push target gaps
            if (tnorm[i] == kGap)
            {
                index_t j = i;
                while (1)
                {
                    const char c = tnorm[++j];
                    if (c != kGap || j > qlen - 1)
                    {
                        if (c == qnorm[i]) { tnorm[i] = c; tnorm[j] = kGap; }
                        break;
                    }
                }
            }
            // push query gaps
            if (qnorm[i] == kGap)
            {
                index_t j = i;
                while (1)
                {
                    const char c = qnorm[++j];
                    if (c != kGap || j > tlen - 1)
                    {
                        if (c == tnorm[i]) { qnorm[i] = c; qnorm[j] = kGap; }
                        break;
                    }
                }
            }
        }
    }
    r_assert(qnorm.size() == tnorm.size());

#ifndef NDEBUG
    int qcnt2 = 0, tcnt2 = 0;
    for (std::string::const_iterator citer = qnorm.begin(); citer != qnorm.end(); ++citer)
        if ((*citer) != kGap) ++qcnt2;
    for (std::string::const_iterator citer = tnorm.begin(); citer != tnorm.end(); ++citer)
        if ((*citer) != kGap) ++tcnt2;
    d_assert(qcnt == qcnt2);
    d_assert(tcnt == tcnt2);
#endif
}

struct CmpExtensionCandidateBySid
{
	bool operator()(const ExtensionCandidate& a, const ExtensionCandidate& b)
	{
		return a.sid < b.sid;
	}
};

void
build_cns_thrd_data_can(ExtensionCandidate* ec_list, 
						const int nec,
						const idx_t min_rid,
						const idx_t max_rid,
						ReadsCorrectionOptions* prco,
						PackedDB* reads,
						std::ostream* out,
					    ConsensusThreadData** ppctd)
{
	const index_t num_reads = max_rid - min_rid + 1;
	const int num_threads = prco->num_threads;
    const index_t num_reads_per_thread = (num_reads + num_threads - 1) / num_threads;
	std::sort(ec_list, ec_list + nec, CmpExtensionCandidateBySid());
	idx_t max_id = min_rid;
	idx_t i = 0, j;
	int tid = 0;
	while (i < nec)
	{
		max_id += num_reads_per_thread;
		j = i + 1;
		while (j < nec && ec_list[j].sid < max_id) ++j;
		ppctd[tid] = new ConsensusThreadData(prco, tid, reads, ec_list + i, j - i, out);
		++tid;
		i = j;
	}
}
