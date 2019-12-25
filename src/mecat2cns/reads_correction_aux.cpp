#include "reads_correction_aux.h"

#include <cassert>	// assert()
#include <string>	// string

void normalize_gaps(const std::string& qstr, const std::string& tstr, std::string& qnorm, std::string& tnorm, const int push)
{
    qnorm.clear();
    tnorm.clear();

#ifndef NDEBUG
    int qcnt = 0, tcnt = 0;
    for (size_t i = 0; i < qstr.size(); ++i)
    {
        const char qc = qstr[i];
        const char tc = tstr[i];
        if (qc != GAP_CHAR) ++qcnt;
        if (tc != GAP_CHAR) ++tcnt;
    }
#endif

    // convert mismatches to indels
    for (size_t i = 0; i < qstr.size(); ++i)
    {
        const char qc = qstr[i];
        const char tc = tstr[i];
        if (qc != tc && qc != GAP_CHAR && tc != GAP_CHAR)
        { qnorm += GAP_CHAR; qnorm += qc; tnorm += tc; tnorm += GAP_CHAR; }
        else
        { qnorm += qc; tnorm += tc; }
    }

    // push gaps to the right, but not pass the end
    if (push)
    {
        int64_t qlen = qnorm.size();
        int64_t tlen = tnorm.size();
        for (int64_t i = 0; i < qlen - 1; ++i)
        {
            // push target gaps
            if (tnorm[i] == GAP_CHAR)
            {
                int64_t j = i;
                while (1)
                {
                    const char c = tnorm[++j];
                    if (c != GAP_CHAR || j > qlen - 1)
                    {
                        if (c == qnorm[i]) { tnorm[i] = c; tnorm[j] = GAP_CHAR; }
                        break;
                    }
                }
            }
            // push query gaps
            if (qnorm[i] == GAP_CHAR)
            {
                int64_t j = i;
                while (1)
                {
                    const char c = qnorm[++j];
                    if (c != GAP_CHAR || j > tlen - 1)
                    {
                        if (c == tnorm[i]) { qnorm[i] = c; qnorm[j] = GAP_CHAR; }
                        break;
                    }
                }
            }
        }
    }
    assert(qnorm.size() == tnorm.size());

#ifndef NDEBUG
    int qcnt2 = 0, tcnt2 = 0;
    for (std::string::const_iterator citer = qnorm.begin(); citer != qnorm.end(); ++citer)
        if ((*citer) != GAP_CHAR) ++qcnt2;
    for (std::string::const_iterator citer = tnorm.begin(); citer != tnorm.end(); ++citer)
        if ((*citer) != GAP_CHAR) ++tcnt2;
    assert(qcnt == qcnt2);
    assert(tcnt == tcnt2);
#endif
}

void allocate_ecs(ConsensusThreadData& data, ExtensionCandidate* const ec_list, const int64_t nec) {
	const int n(data.rco.num_threads);
	// split by number of ec's, rather than reads, since reads ids
	// are not contiguous and we could get empty lists
	for (int64_t i(0), k(0); k != n; ++k) {
		const int64_t start(i);
		// drop fractions here, as we'll likely add a few more ec's below
		i += (nec - i) / (n - k);
		if (i != nec) {			// include all ec's for the last read
			const int final_sid(ec_list[i].sid);
			for (++i; i != nec && ec_list[i].sid == final_sid; ++i) { }
		}
		data.data[k].num_candidates = i - start;
		data.data[k].candidates = ec_list + start;
	}
}

void allocate_ecs(ConsensusThreadData& data, ExtensionCandidateCompressed* const ec_list, const int64_t nec) {
	const int64_t n(data.rco.num_threads);
	// split by number of ec's, rather than reads, since reads ids
	// are not contiguous and we could get empty lists
	for (int64_t i(0), k(0); k != n; ++k) {
		const int64_t start(i);
		// drop fractions here, as we'll likely add a few more ec's below
		i += (nec - i) / (n - k);
		if (i != nec) {			// include all ec's for the last read
			const uint32_t final_sid(ec_list[i].sid);
			for (++i; i != nec && ec_list[i].sid == final_sid; ++i) { }
		}
		data.data[k].num_candidates = i - start;
		data.data[k].candidates = ec_list + start;
	}
}
