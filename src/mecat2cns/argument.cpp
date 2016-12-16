#include "argument.h"

#include <sstream>

#define error_and_exit(msg) { std::cerr << msg << "\n"; abort(); }

#define no_argument \
    do { \
        std::ostringstream err_msg; \
        err_msg << "no argument is supplied to option \'" << (*arg_name) << "\'"; \
        error_and_exit(err_msg.str()); \
    } while(0);

#define argument_process_error \
    do { \
        std::ostringstream err_msg; \
        err_msg << "argument \'" << argv[0] << "\' for option \'" << (*arg_name) << "\' processed failed."; \
        error_and_exit(err_msg.str()); \
    } while(0);

int IntegerArgument::ProcessArgument(int argc, char** argv)
{
    if (!argc) no_argument;
    std::istringstream ins(argv[0]);
    ins >> val;
    if (!ins) argument_process_error;
    return 1;
}

int DoubleArgument::ProcessArgument(int argc, char** argv)
{
    if (!argc) no_argument;
    std::istringstream ins(argv[0]);
    ins >> val;
    if (!ins) argument_process_error;
    return 1;
}

int BooleanArgument::ProcessArgument(int argc, char** argv)
{
    if (!need_arg)
    {
        val = true;
        return 0;
    }
    else
    {
        if (!argc) no_argument;

        if (argv[0][0] == '0') 
        {
            val = false;
            return 1;
        }
        else if (argv[0][0] == '1')
        {
            val = true;
            return 1;
        }
        else
        {
            std::ostringstream err_msg;
            err_msg << "argument to option \'" << (*arg_name) << "\' must be \'0\' or \'1\'";
            error_and_exit(err_msg.str());
        }
    }
}

int StringArgument::ProcessArgument(int argc, char** argv)
{
    if (!argc) no_argument;

    val = argv[0];
    return 1;
}
