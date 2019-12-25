#ifndef PW_IMPL_H
#define PW_IMPL_H

#include "../common/alignment.h"
#include "../common/lookup_table.h"

#include "pw_options.h"	// options_t
#include <string>	// string

void process_one_volume(const options_t* options, int svid, int evid, const std::string &volume_results_name, volume_names_t* vn);

#endif // PW_IMPL_H
