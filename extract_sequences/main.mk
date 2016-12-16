ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET       := libes.a

SOURCES      := \
		src/common/AS_GKP_buildPartition.C \
		src/common/AS_GKP_checkFrag.C \
 		src/common/AS_GKP_checkLibrary.C \
		src/common/AS_GKP_checkLink.C \
		src/common/AS_GKP_checkPlace.C \
		src/common/AS_GKP_dump.C \
		src/common/AS_GKP_edit.C \
		src/common/AS_GKP_errors.C \
		src/common/AS_GKP_illumina.C \
		src/common/AS_global.C \
		src/common/AS_MSG_pmesg1.C \
		src/common/AS_MSG_pmesg2.C \
		src/common/AS_MSG_pmesg.C \
		src/common/AS_PER_encodeSequenceQuality.C \
		src/common/AS_PER_genericStore.C \
		src/common/AS_PER_gkLibrary.C \
		src/common/AS_PER_gkStore.C \
		src/common/AS_PER_gkStore_clearRange.C \
		src/common/AS_PER_gkStore_fragments.C \
		src/common/AS_PER_gkStore_IID.C \
		src/common/AS_PER_gkStore_partition.C \
		src/common/AS_PER_gkStore_PLC.C \
		src/common/AS_PER_gkStore_stats.C \
		src/common/AS_PER_gkStore_UID.C \
		src/common/AS_PER_gkStream.C \
		src/common/AS_UTL_alloc.C \
		src/common/AS_UTL_fasta.C \
		src/common/AS_UTL_fileIO.C \
		src/common/AS_UTL_Hash.C \
		src/common/AS_UTL_heap.C \
		src/common/AS_UTL_reverseComplement.C \
		src/common/AS_UTL_stackTrace.C \
		src/common/AS_UTL_UID.C \
		src/common/AS_UTL_Var.C \


SRC_INCDIRS  := src 

SUBMAKEFILES := src/extract_sequences/extract_sequences.mk \
				src/fasta2fastq/fasta2fastq.mk \
				src/fastqToCA/fastqToCA.mk\
				src/gatekeeper/gatekeeper.mk
