#ifndef MEAP_CORRECTION_H
#define MEAP_CORRECTION_H

#include "reads_correction_aux.h"

namespace ns_meap_cns {

void
meap_consensus_one_read_m4(ConsensusThreadData* ctd, const index_t read_id, const index_t sid, const index_t eid);
	
void
meap_consensus_one_read_can(ConsensusThreadData* ctd, const index_t read_id, const index_t sid, const index_t eid);

} // namespace ns_meap_cns

#endif // MEAP_CORRECTION_H
