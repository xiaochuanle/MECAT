ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := mecat2pw
SOURCES  := pw.cpp pw_impl.cpp pw_options.cpp

SRC_INCDIRS  := ../common .

TGT_CXXFLAGS := -D _FILE_OFFSET_BITS=64 -std=c++11 -O3 -march=native -flto -fno-fat-lto-objects -fno-builtin
TGT_LDFLAGS  := -L${TARGET_DIR} ${TGT_CXXFLAGS}
TGT_LDLIBS   := -lmecat
TGT_PREREQS  := libmecat.a

SUBMAKEFILES :=
