
#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET       := libcanu.a

SOURCES      := AS_global.C \
                \
                AS_UTL/AS_UTL_decodeRange.C \
                AS_UTL/AS_UTL_fasta.C \
                AS_UTL/AS_UTL_fileIO.C \
                AS_UTL/AS_UTL_reverseComplement.C \
                AS_UTL/AS_UTL_stackTrace.C \
                \
                AS_UTL/AS_UTL_alloc.C \
                \
                AS_UTL/bitEncodings.C \
                AS_UTL/bitPackedFile.C \
                AS_UTL/bitPackedArray.C \
                AS_UTL/dnaAlphabets.C \
                AS_UTL/md5.C \
                AS_UTL/mt19937ar.C \
                AS_UTL/readBuffer.C \
                AS_UTL/speedCounter.C \
                AS_UTL/stddev.C \
                AS_UTL/sweatShop.C \
                AS_UTL/timeAndSize.C \
                AS_UTL/kMer.C \
                \
                stores/gkStore.C \
                stores/gkStoreEncode.C \
                \
                stores/ovOverlap.C \
                stores/ovStore.C \
                stores/ovStoreFile.C \
                \
                stores/tgStore.C \
                stores/tgTig.C \
                stores/tgTigSizeAnalysis.C \
                stores/tgTigMultiAlignDisplay.C \
                \
                meryl/libmeryl.C \
                \
                overlapInCore/overlapReadCache.C \
                \
                overlapErrorAdjustment/analyzeAlignment.C \
                \
                overlapInCore/liboverlap/Binomial_Bound.C \
                overlapInCore/liboverlap/Display_Alignment.C \
                overlapInCore/liboverlap/prefixEditDistance.C \
                overlapInCore/liboverlap/prefixEditDistance-allocateMoreSpace.C \
                overlapInCore/liboverlap/prefixEditDistance-extend.C \
                overlapInCore/liboverlap/prefixEditDistance-forward.C \
                overlapInCore/liboverlap/prefixEditDistance-reverse.C \
                \
                utgcns/libNDalign/NDalign.C \
                \
                utgcns/libNDalign/Binomial_Bound.C \
                utgcns/libNDalign/NDalgorithm.C \
                utgcns/libNDalign/NDalgorithm-allocateMoreSpace.C \
                utgcns/libNDalign/NDalgorithm-extend.C \
                utgcns/libNDalign/NDalgorithm-forward.C \
                utgcns/libNDalign/NDalgorithm-reverse.C \
                \
                utgcns/libcns/abAbacus-addRead.C \
                utgcns/libcns/abAbacus-appendBases.C \
                utgcns/libcns/abAbacus-applyAlignment.C \
                utgcns/libcns/abAbacus-baseCall.C \
                utgcns/libcns/abAbacus-mergeRefine.C \
                utgcns/libcns/abAbacus-refine.C \
                utgcns/libcns/abAbacus-refreshMultiAlign.C \
                utgcns/libcns/abAbacus.C \
                utgcns/libcns/abColumn.C \
                utgcns/libcns/abMultiAlign.C \
                utgcns/libcns/unitigConsensus.C \
                utgcns/libpbutgcns/Alignment.C	\
                utgcns/libpbutgcns/AlnGraphBoost.C  \
                utgcns/libpbutgcns/SimpleAligner.C \
                utgcns/libNDFalcon/dw.C

SRC_INCDIRS  := . \
                AS_UTL \
                stores \
                alignment \
                utgcns/libNDalign \
                utgcns/libcns \
                utgcns/libpbutgcns \
                utgcns/libNDFalcon \
                utgcns/libboost \
                meryl/libleaff \
                overlapInCore \
                overlapInCore/liboverlap

SUBMAKEFILES := stores/gatekeeperCreate.mk \
                stores/gatekeeperDumpFASTQ.mk \
                stores/gatekeeperDumpMetaData.mk \
                stores/gatekeeperPartition.mk \
		\
                stores/ovStoreBuild.mk \
                stores/ovStoreBucketizer.mk \
                stores/ovStoreDump.mk \
                stores/ovStoreStats.mk \
		\
                stores/tgStoreDump.mk \
                stores/tgStoreLoad.mk \
                stores/tgStoreCoverageStat.mk \
		\
                overlapBasedTrimming/trimReads.mk \
                overlapBasedTrimming/splitReads.mk \
                \
                overlapErrorAdjustment/findErrors.mk \
                overlapErrorAdjustment/correctOverlaps.mk \
                \
                bogart/bogart.mk \
                bogart/buildGraph.mk \
                \
                utgcns/utgcns.mk \
                \
		mecat2asmpw/fastaconvert.mk \
		mecat2asmpw/mecat2asmpw.mk \
                mecat2asmpw/mecat2asmpw50.mk \
                mecat2asmpw/mecat2trimpw.mk \
                mecat2asmpw/mecat2trimpw50.mk \
		mecat2asmpw/mecat2asmpwConvert.mk \
		partition_reads/partition_reads.mk
