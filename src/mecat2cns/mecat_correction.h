#ifndef MEAP_CORRECTION_H
#define MEAP_CORRECTION_H

#include <stdint.h>	// int64_t

#include "reads_correction_aux.h"	// ConsensusPerThreadData, ConsensusThreadData

void consensus_one_read_m4_pacbio(ConsensusThreadData& ctd, ConsensusPerThreadData &cptd, const int64_t read_id, const int64_t sid, const int64_t eid);

void consensus_one_read_m4_nanopore(ConsensusThreadData& ctd, ConsensusPerThreadData &cptd, const int64_t read_id, const int64_t sid, const int64_t eid);

void consensus_one_read_can_pacbio(ConsensusThreadData& ctd, ConsensusPerThreadData &cptd, const int64_t read_id, const int64_t sid, const int64_t eid);

void consensus_one_read_can_nanopore(ConsensusThreadData& ctd, ConsensusPerThreadData &cptd, const int64_t read_id, const int64_t sid, const int64_t eid);

#endif // MEAP_CORRECTION_H
