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

void PrintM5Record(std::ostream& out, const M5Record& m5, const int printAln)
{
	out << "(" << m5qid(m5) << ", " << m5qsize(m5) << ", " << m5qoff(m5) << ", " << m5qend(m5) << ", " << m5qdir(m5) << ")"
		<< " x "
		<< "(" << m5sid(m5) << ", " << m5ssize(m5) << ", " << m5soff(m5) << ", " << m5send(m5) << ", " << m5sdir(m5) << ")"
		<< ", score = " << m5score(m5) 
		<< ", mapq = " << m5mapq(m5)
		<< ", (" << m5mat(m5) << ", " << m5mis(m5) << ", " << m5ins(m5) << ", " << m5dels(m5) << ")\n";
	
	if (printAln)
	{
		out << m5qaln(m5) << "\n";
		out << m5pat(m5) << "\n";
		out << m5saln(m5) << "\n";
	}
}

void InitM5Record(M5Record& m5)
{
	m5qid(m5) = m5qsize(m5) = m5qoff(m5) = m5qend(m5) = 0;
	m5sid(m5) = m5ssize(m5) = m5soff(m5) = m5send(m5) = 0;
	m5qdir(m5) = m5sdir(m5) = 2;
	m5score(m5) = m5mat(m5) = m5mis(m5) = m5ins(m5) = m5dels(m5) = -1;
	m5mapq(m5) = 0;
	m5qaln(m5) = m5pat(m5) = m5saln(m5) = NULL;
}

void DestroyM5Record(M5Record& m5)
{
	if (m5qaln(m5)) delete[] m5qaln(m5);
	InitM5Record(m5);
}

/*
std::istream& operator>>(std::istream& in, ReferenceMapping& rm)
{
	if (!in) return in;
	if (!(in >> rmqid(rm))) return in;
	in >> rmsid(rm)
	   >> rmqdir(rm)
	   >> rmqoff(rm)
	   >> rmqend(rm)
	   >> rmqext(rm)
	   >> rmqsize(rm)
	   >> rmsdir(rm)
	   >> rmsoff(rm)
	   >> rmsend(rm)
	   >> rmsext(rm)
	   >> rmssize(rm);
	return in;
}

std::ostream& operator<<(std::ostream& out, const ReferenceMapping& rm)
{
	constexpr char delim = '\t';
   out << rmqid(rm) << delim
	   << rmsid(rm) << delim
	   << rmqdir(rm) << delim
	   << rmqoff(rm) << delim
	   << rmqend(rm) << delim
	   << rmqext(rm) << delim
	   << rmqsize(rm) << delim
	   << rmsdir(rm) << delim
	   << rmsoff(rm) << delim
	   << rmsend(rm) << delim
	   << rmsext(rm) << delim
	   << rmssize(rm)
	   << "\n";
   return out;
}
*/
