CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing
CC = gcc

all: dextract dexta undexta dexqv undexqv

dextract: DEXTRACTOR
	${CC} $(CFLAGS) -I$(PATH_HDF5)/include -L${PATH_HDF5}/lib -o ${BUILD_DIR}/dextract DEXTRACTOR/dextract.c DEXTRACTOR/DB.c DEXTRACTOR/QV.c -lhdf5

dexta: DEXTRACTOR
	${CC} ${CFLAGS} -o ${BUILD_DIR}/dexta DEXTRACTOR/dexta.c DEXTRACTOR/DB.c DEXTRACTOR/QV.c

undexta: DEXTRACTOR
	${CC} ${CFLAGS} -o ${BUILD_DIR}/undexta DEXTRACTOR/undexta.c DEXTRACTOR/DB.c DEXTRACTOR/QV.c

dexqv: DEXTRACTOR
	${CC} ${CFLAGS} -o ${BUILD_DIR}/dexqv DEXTRACTOR/dexqv.c DEXTRACTOR/DB.c DEXTRACTOR/QV.c

undexqv: DEXTRACTOR
	${CC} ${CFLAGS} -o ${BUILD_DIR}/undexqv DEXTRACTOR/undexqv.c DEXTRACTOR/DB.c DEXTRACTOR/QV.c

DEXTRACTOR:
	git clone https://github.com/PacificBiosciences/DEXTRACTOR.git
