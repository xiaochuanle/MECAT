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

struct M5Record
{
	idx_t qid;		// 1) qname
	idx_t qsize;  	// 2) qlength
	idx_t qstart; 	// 3) qstart
	idx_t qend;		// 4) qend
	int     qdir;		// 5) qstrand
	idx_t sid;		// 6) sname
	idx_t ssize;		// 7) slength
	idx_t sstart;		// 8) sstart
	idx_t send;		// 9) send
	int     sdir;		// 10) sstrand
	int     score;		// 11) score
	int     mat;		// 12) match
	int		mis; 	    // 13) mismatch
	int     ins;		// 14) insertion
	int		dels;		// 15) deletion
	int		mapq;		// 16) mapQ
	char*   pm_q;		// 17) aligned query
	char*   pm_p;		// 18) aligned pattern
	char*   pm_s;		// 19) aligned subject
	double	ident;		// 20) identity percentage
	idx_t	qext;
	idx_t 	sext;
};

#define m5qid(m) 		((m).qid)
#define m5qsize(m)		((m).qsize)
#define m5qoff(m)		((m).qstart)
#define m5qend(m)		((m).qend)
#define m5qdir(m)		((m).qdir)
#define m5sid(m)		((m).sid)
#define m5ssize(m)		((m).ssize)
#define m5soff(m)		((m).sstart)
#define m5send(m)		((m).send)
#define m5sdir(m)		((m).sdir)
#define m5score(m)		((m).score)
#define m5mat(m)		((m).mat)
#define m5mis(m)		((m).mis)
#define m5ins(m)		((m).ins)
#define m5dels(m)		((m).dels)
#define m5mapq(m)		((m).mapq)
#define m5qaln(m)		((m).pm_q)
#define m5pat(m)		((m).pm_p)
#define m5saln(m)		((m).pm_s)
#define m5ident(m)		((m).ident)
#define m5qext(m)		((m).qext)
#define m5sext(m)		((m).sext)

void PrintM5Record(std::ostream& out, const M5Record& m5, const int printAln);
void InitM5Record(M5Record& m5);
void DestroyM5Record(M5Record& m5);

inline M5Record* NewM5Record(idx_t maxAlnSize)
{
	M5Record* m5 = new M5Record;
	m5qaln(*m5) = new char[maxAlnSize];
	m5saln(*m5) = new char[maxAlnSize];
	m5pat(*m5) = new char[maxAlnSize];
	
	return m5;
}

inline M5Record* DeleteM5Record(M5Record* m5)
{
	if (!m5) return NULL;
	delete[] m5qaln(*m5);
	delete[] m5saln(*m5);
	delete[] m5pat(*m5);
	delete m5;
	return NULL;
}

inline int M5RecordOvlpSize(const M5Record& m)
{
	int oq = m5qend(m) - m5qoff(m);
	int os = m5send(m) - m5soff(m);
	return std::max(oq, os);
}

//struct Overlap
//{
//    index_t qid, qoff, qend, qsize, qext;
//    int qdir;
//    index_t sid, soff, send, ssize, sext;
//    int sdir;
//};

typedef ExtensionCandidate Overlap;

struct CompareOverlapBySid
{
    bool operator()(const Overlap& a, const Overlap& b)
    {
        return a.sid < b.sid;
    }
};

struct CnsResult
{
	idx_t id;
	idx_t range[2];
	std::string seq;
};

#endif // ALIGNMENT_H
