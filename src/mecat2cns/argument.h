#ifndef ARGUMENT_H
#define ARGUMENT_H

#include <string>

#include "../common/defs.h"

class Argument
{
public:
    Argument(const std::string* an, const std::string* ad) : arg_name(an), arg_desc(ad) {}
    virtual ~Argument() {}
    virtual int ProcessArgument(int argc, char** argv) = 0;
    
protected:
    const std::string* arg_name;
    const std::string* arg_desc;
};

class IntegerArgument : public Argument
{
public:
    IntegerArgument(const std::string* an, const std::string* ad, const index_t v) : Argument(an, ad), val(v) {}
    virtual ~IntegerArgument() {}
    virtual int ProcessArgument(int argc, char** argv);
    index_t value() { return val; }

private:
    index_t val;
};

class DoubleArgument : public Argument
{
public:
    DoubleArgument(const std::string* an, const std::string* ad, double v) : Argument(an, ad), val(v) {}
    virtual ~DoubleArgument() {}
    virtual int ProcessArgument(int argc, char** argv);
    double value() { return val; }

private:
    double val;
};

class BooleanArgument : public Argument
{
public:
    BooleanArgument(const std::string* an, const std::string* ad, const bool v, const bool na)
        : Argument(an, ad), val(v), need_arg(na) {}
    virtual ~BooleanArgument() {}
    virtual int ProcessArgument(int argc, char** argv);
    bool value() { return val; }

private:
    bool val;
    bool need_arg;
};

class StringArgument : public Argument
{
public:
    StringArgument(const std::string* an, const std::string* ad, const char* v) : Argument(an, ad), val(v) {}
    ~StringArgument() {}
    virtual int ProcessArgument(int argc, char** argv);
    const char* value() { return val; }

private:
    const char* val;
};

#endif // ARGUMENT_H
