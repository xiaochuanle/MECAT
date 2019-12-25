ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET       := libmecat.a

SOURCES      := common/alignment.cpp \
		common/buffer_line_iterator.cpp \
		common/defs.cpp \
		common/diff_gapalign.cpp \
		common/fasta_reader.cpp \
		common/gapalign.cpp \
		common/lookup_table.cpp \
		common/sequence.cpp \
		common/split_database.cpp \
		common/xdrop_gapalign.cpp

SRC_INCDIRS  := common \

SUBMAKEFILES := mecat2pw/pw.mk \
		mecat2ref/mecat2ref.mk \
		mecat2cns/mecat2cns.mk \
		filter_reads/filter_reads.mk

# Note: -O2 performed about 1% worse than -O3

TGT_CXXFLAGS := -D _FILE_OFFSET_BITS=64 -std=c++11 -O3 -mfpmath=sse -march=native -flto -fno-fat-lto-objects -fno-builtin -mmmx -msse -msse2 -mssse3 -msse4.1 -msse4.2 -mavx -maes -mpopcnt -mfxsr -mxsave -mxsaveopt
TGT_LDFLAGS  := ${TGT_CXXFLAGS}
ARFLAGS      += --plugin /mnt/local/gnu/libexec/gcc/x86_64-pc-linux-gnu/9.1.0/liblto_plugin.so.0
