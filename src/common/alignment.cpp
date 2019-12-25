#include "alignment.h"

#include <sstream>
#include <string>

using namespace std;

std::istream&
operator>>(std::istream& in, ExtensionCandidate& ec)
{
	std::string line;
	if (!getline(in, line)) return in;
	istringstream ins(line);
	ins >> ec.qid >> ec.sid >> ec.qdir >> ec.sdir >> ec.qext >> ec.sext >> ec.score >> ec.qsize >> ec.ssize;
	return in;
}

std::ostream&
operator<<(std::ostream& out, const ExtensionCandidate& ec)
{
	const char delim = '\t';
	out << ec.qid << delim
		<< ec.sid << delim
		<< ec.qdir << delim
		<< ec.sdir << delim
		<< ec.qext << delim
		<< ec.sext << delim
		<< ec.score << delim
		<< ec.qsize << delim
		<< ec.ssize << std::endl;
	return out;
}

std::istream& operator>>(std::istream& in, M4Record& m4)
{
	std::string line;
	if (!getline(in, line)) return in;
	std::istringstream ins(line);
	ins >> m4qid(m4)
		>> m4sid(m4)
		>> m4ident(m4)
		>> m4vscore(m4)
		>> m4qdir(m4)
		>> m4qoff(m4)
		>> m4qend(m4)
		>> m4qsize(m4)
		>> m4sdir(m4)
		>> m4soff(m4)
		>> m4send(m4)
		>> m4ssize(m4);
	
	// qext and sext are optional
	if (!(ins >> m4qext(m4))) return in;
	ins >> m4sext(m4);
	return in;
}

std::ostream& operator<<(std::ostream& out, const M4Record& m4)
{
	const char sep = '\t';
	
	out << m4qid(m4)    << sep
	    << m4sid(m4)    << sep
	    << m4ident(m4) << sep
	    << m4vscore(m4)   << sep
	    << m4qdir(m4)   << sep
	    << m4qoff(m4)   << sep
	    << m4qend(m4)   << sep
	    << m4qsize(m4)  << sep
	    << m4sdir(m4)   << sep
	    << m4soff(m4)   << sep
	    << m4send(m4)   << sep
	    << m4ssize(m4)	<< sep
		<< m4qext(m4)	<< sep
		<< m4sext(m4)	<< "\n";
	
	return out;
}

const int64_t ExtensionCandidateCompressed::max_value = std::numeric_limits<uint32_t>::max();
const int64_t ExtensionCandidateCompressed::max_qext = std::numeric_limits<uint32_t>::max() >> 1;
