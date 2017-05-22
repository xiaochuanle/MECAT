#ifndef INDEL_DETECTION_SMART_ASSERT_H
#define INDEL_DETECTION_SMART_ASSERT_H

#include <string>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>
#include <map>
#include <cassert>
#include <cstdlib>

namespace smart_assert_private
{
    // allows finding if a value is of type 'const char*' and is null
    // if so, we cannot print it to an ostream directly !!!
    template <class T>
    struct is_null_finder
    {
        bool is(const T&) const
        {
            return false;
        }
    };

    template <>
    struct is_null_finder<char*>
    {
        bool is(char* const& val)
        {
            return val == 0;
        }
    };

    template <>
    struct is_null_finder<const char*>
    {
        bool is(const char* const& val)
        {
            return val == 0;
        }
    };
}

class assert_context
{
public:
    typedef std::string string;
    typedef std::pair<string, string> val_and_str;
    typedef std::vector<val_and_str> vals_array;
public:
    void set_file_line(const char* file, int line)
    {
        m_file = file;
        m_line = line;
    }
    const string& get_context_file() const { return m_file; }
    int get_context_line() const { return m_line; }

    void set_expr(const string& expr) { m_expr = expr; }
    const string& get_expr() const { return m_expr; }

    const vals_array& get_vals_array() const { return m_vals; }
    void add_val(const string& val, const string& str)
    {
        m_vals.push_back(val_and_str(val, str));
    }

private:
    string m_file;
    int m_line;
    string m_expr;
    vals_array m_vals;
};

struct Assert
{
public:
    Assert& SMART_ASSERT_A;
    Assert& SMART_ASSERT_B;

    Assert(const char* expr)
            : SMART_ASSERT_A(*this),
              SMART_ASSERT_B(*this),
              m_need_handling(true)
    {
        m_context.set_expr(expr);
    }
    Assert(const Assert& other)
            : SMART_ASSERT_A(*this),
              SMART_ASSERT_B(*this),
              m_context(other.m_context),
              m_need_handling(true)
    {
        other.m_need_handling = false;
    }
    ~Assert()
    {
        if(m_need_handling) handle_assert();
    }

    template <class type>
    Assert& print_current_val(const type& val, const char* msg)
    {
        std::ostringstream out;
        smart_assert_private::is_null_finder<type> f;
        bool bIsNull = f.is(val);
        if (!bIsNull) out << val;
        else out << "null";
        m_context.add_val(out.str(), msg);
        return *this;
    }

    Assert& print_context(const char* file, int line)
    {
        m_context.set_file_line(file, line);
        return *this;
    }

private:
    void handle_assert()
    {
        std::ostream& out = std::cerr;
        out << "Assertion failed at \'" << m_context.get_context_file() << ":" << m_context.get_context_line() << "\'\n";
        out << "Expression: \'" << m_context.get_expr() << "\'\n";
        const assert_context::vals_array& aVals = m_context.get_vals_array();
        if (!aVals.empty())
        {
            out << "Values: " << "\n";
            assert_context::vals_array::const_iterator first = aVals.begin(), last = aVals.end();
            while (first != last)
            {
                out << "\t" << first->second << " = \'" << first->first << "\'\n";
                ++first;
            }
        }
        abort();
    }

private:
    assert_context m_context;
    mutable bool m_need_handling;
};

inline Assert
make_assert(const char* expr) { return Assert(expr); }

#define SMART_ASSERT_OP(x, next) \
    SMART_ASSERT_A.print_current_val((x), #x).SMART_ASSERT_ ## next \

#define SMART_ASSERT_A(x) SMART_ASSERT_OP(x, B)
#define SMART_ASSERT_B(x) SMART_ASSERT_OP(x, A)

#define D_ASSERT(expr) \
    if ((expr)); \
    else make_assert(#expr).print_context(__FILE__, __LINE__).SMART_ASSERT_A \

#define R_ASSERT(expr) \
    if ((expr)); \
    else make_assert(#expr).print_context(__FILE__, __LINE__).SMART_ASSERT_A \

#define d_assert D_ASSERT
#define r_assert R_ASSERT

#endif //INDEL_DETECTION_SMART_ASSERT_H
