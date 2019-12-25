#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <fstream>

#include "defs.h"

struct ExtensionCandidate
{
	int qdir, qid, qext, qsize, qoff, qend;
	int sdir, sid, sext, ssize, soff, send;
	int score;
};

typedef ExtensionCandidate Overlap;

// candidates only use a few of these values, so make a smaller structure for them;
// sid and qid are used a lot, and we sort on score, so stash qdir inside msb of qext

struct ExtensionCandidateCompressed {
	uint32_t sid, qid, sext, qext_, score;
	int qdir() const {
		return qext_ & MSB_ ? REV : FWD;
	}
	int qext() const {
		return qext_ & ~MSB_;
	}
	void set_qext(const uint32_t new_qext, const int new_qdir) {
		qext_ = new_qdir == FWD ? new_qext : new_qext | MSB_;
	}
	void set(const ExtensionCandidate& a) {
		sid = a.sid;
		qid = a.qid;
		sext = a.sext;
		// sdir is forced to FWD, so swap qdir is sdir is REV
		qext_ = a.qdir == a.sdir ? a.qext : a.qext | MSB_;
		score = a.score;
	}
	void set_swap(const ExtensionCandidate& a) {		// swap s and q
		sid = a.qid;
		qid = a.sid;
		sext = a.qext;
		qext_ = a.sdir == a.qdir ? a.sext : a.sext | MSB_;
		score = a.score;
	}
	// (these two are set in alignment.cpp, as std::numeric_limits<>
	// can't be used at compile time)
	// used to check conversions from type int; need to use int64_t in case
	// int is only int32_t in size (which wouldn't hold uint32_t max)
	static const int64_t max_value;
	// account for using MSB for qdir
	static const int64_t max_qext;
    private:
	static const uint32_t MSB_ = 1 << (sizeof(uint32_t) * 8 - 1);
};

std::istream&
operator>>(std::istream& in, ExtensionCandidate& ec);

std::ostream&
operator<<(std::ostream& out, const ExtensionCandidate& ec);

struct M4Record
{
	idx_t qid;
	idx_t sid;
	double  ident;
	int     vscore;
	int     qdir;
	idx_t qoff;
	idx_t qend;
	idx_t qsize;
	int     sdir;
	idx_t soff;
	idx_t send;
	idx_t ssize;
	idx_t qext;
	idx_t sext;
};

#define m4qid(m)		((m).qid)
#define m4sid(m)		((m).sid)
#define m4ident(m)		((m).ident)
#define m4vscore(m)		((m).vscore)
#define m4qdir(m)		((m).qdir)
#define m4qoff(m)		((m).qoff)
#define m4qend(m)		((m).qend)
#define m4qsize(m)		((m).qsize)
#define m4sdir(m)		((m).sdir)
#define m4soff(m)		((m).soff)
#define m4send(m)		((m).send)
#define m4ssize(m)		((m).ssize)
#define m4qext(m)		((m).qext)
#define m4sext(m)		((m).sext)

std::istream& operator>>(std::istream& in, M4Record& m4);
std::ostream& operator<<(std::ostream& out, const M4Record& m4);

template <class ListT>
void load_m4records_from_m4_file(const char* m4_file_name, ListT& m4v)
{
	std::ifstream in;
	open_fstream(in, m4_file_name, std::ios::in);
	M4Record m4;
	m4qext(m4) = INVALID_IDX;
	m4sext(m4) = INVALID_IDX;
	while (in >> m4)
	{
		m4v.push_back(m4);
	}
	close_fstream(in);
}

inline void reverse_m4record(const M4Record& m4, M4Record& nm4)
{
	m4qid(nm4) = m4sid(m4);
	m4sid(nm4) = m4qid(m4);
	m4ident(nm4) = m4ident(m4);
	m4vscore(nm4) = m4vscore(m4);
	m4qdir(nm4) = m4sdir(m4);
	m4qoff(nm4) = m4soff(m4);
	m4qend(nm4) = m4send(m4);
	m4qsize(nm4) = m4ssize(m4);
	m4sdir(nm4) = m4qdir(m4);
	m4soff(nm4) = m4qoff(m4);
	m4send(nm4) = m4qend(m4);
	m4ssize(nm4) = m4qsize(m4);
	m4qext(nm4) = m4sext(m4);
	m4sext(nm4) = m4qext(m4);
}

inline void normalize_m4record(const M4Record& src, const bool subject_is_target, M4Record& dst)
{
	if (subject_is_target)
		dst = src;
	else
		reverse_m4record(src, dst);
	
	if (m4sdir(dst) == REV)
	{
		m4sdir(dst) = REVERSE_STRAND(m4sdir(dst));
		m4qdir(dst) = REVERSE_STRAND(m4qdir(dst));
	}
}

inline idx_t M4RecordOverlapSize(const M4Record& m4)
{
	const idx_t qs = m4qend(m4) - m4qoff(m4);
	const idx_t ss = m4send(m4) - m4soff(m4);
	return MIN(qs, ss);
}

inline bool overlap_aend_is_5prime(const M4Record& m4, const idx_t invalid_end_size)
{
	const idx_t ahg5 = m4qoff(m4);
	const idx_t ahg3 = m4qsize(m4) - m4qend(m4);
	const idx_t bhg5 = m4soff(m4);
	const idx_t bhg3 = m4ssize(m4) - m4send(m4);
	return (bhg5 > 0) && (ahg3 > 0) && (ahg5 <= invalid_end_size) && (bhg3 <= invalid_end_size);
}

inline bool overlap_aend_is_3prime(const M4Record& m4, const idx_t invalid_end_size)
{
	const idx_t ahg5 = m4qoff(m4);
	const idx_t ahg3 = m4qsize(m4) - m4qend(m4);
	const idx_t bhg5 = m4soff(m4);
	const idx_t bhg3 = m4ssize(m4) - m4send(m4);
	return (ahg5 > 0) && (bhg3 > 0) && (bhg5 <= invalid_end_size) && (ahg3 <= invalid_end_size);
}

inline bool overlap_5prime_is_partial(const M4Record& m4, const idx_t invalid_end_size)
{
	const idx_t ahg5 = m4qoff(m4);
	//const idx_t ahg3 = m4qsize(m4) - m4qend(m4);
	const idx_t bhg5 = m4soff(m4);
	//const idx_t bhg3 = m4ssize(m4) - m4send(m4);
	return (ahg5 > invalid_end_size) && (bhg5 > invalid_end_size);
}

inline bool overlap_3prime_is_partial(const M4Record& m4, const idx_t invalid_end_size)
{
	//const idx_t ahg5 = m4qoff(m4);
	const idx_t ahg3 = m4qsize(m4) - m4qend(m4);
	//const idx_t bhg5 = m4soff(m4);
	const idx_t bhg3 = m4ssize(m4) - m4send(m4);
	return (ahg3 > invalid_end_size) && (bhg3 > invalid_end_size);
}

inline bool overlap_is_partial(const M4Record& m4, const idx_t invalid_end_size)
{
	return overlap_5prime_is_partial(m4, invalid_end_size) || overlap_3prime_is_partial(m4, invalid_end_size);
}

inline bool overlap_a_is_contained(const M4Record& m4, const idx_t invalid_end_size)
{
	const idx_t ahg5 = m4qoff(m4);
	const idx_t ahg3 = m4qsize(m4) - m4qend(m4);
	//const idx_t bhg5 = m4soff(m4);
	//const idx_t bhg3 = m4ssize(m4) - m4send(m4);
	return (ahg5 <= invalid_end_size) && (ahg3 <= invalid_end_size);
}

inline bool overlap_a_is_container(const M4Record& m4, const idx_t invalid_end_size)
{
	//const idx_t ahg5 = m4qoff(m4);
	//const idx_t ahg3 = m4qsize(m4) - m4qend(m4);
	const idx_t bhg5 = m4soff(m4);
	const idx_t bhg3 = m4ssize(m4) - m4send(m4);
	return (bhg5 <= invalid_end_size) && (bhg3 <= invalid_end_size);
}

inline void
m4_to_candidate(const M4Record& m4, ExtensionCandidate& ec)
{
	ec.qdir = m4qdir(m4);
	ec.qid = m4qid(m4);
	ec.qext = m4qext(m4);
	ec.qsize = m4qsize(m4);
	ec.qoff = m4qoff(m4);
	ec.qend = m4qend(m4);
	ec.sdir = m4sdir(m4);
	ec.sid = m4sid(m4);
	ec.sext = m4sext(m4);
	ec.ssize = m4ssize(m4);
	ec.soff = m4soff(m4);
	ec.send = m4send(m4);
	ec.score = m4vscore(m4);
}

#endif // ALIGNMENT_H
