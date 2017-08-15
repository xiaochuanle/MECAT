#ifndef MEAP_CORRECTION_H
#define MEAP_CORRECTION_H

#include "reads_correction_aux.h"

namespace ns_meap_cns {

void
consensus_one_read_m4_pacbio(ConsensusThreadData* ctd, const idx_t read_id, const idx_t sid, const idx_t eid);

void
consensus_one_read_m4_nanopore(ConsensusThreadData* ctd, const idx_t read_id, const idx_t sid, const idx_t eid);

void
consensus_one_read_can_pacbio(ConsensusThreadData* ctd, const idx_t read_id, const idx_t sid, const idx_t eid);

void
consensus_one_read_can_nanopore(ConsensusThreadData* ctd, const idx_t read_id, const idx_t sid, const idx_t eid);

} // namespace ns_meap_cns

#endif // MEAP_CORRECTION_H
