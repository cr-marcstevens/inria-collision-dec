AM_CFLAGS = $(M4RI_CFLAGS)
LIBS = $(M4RI_LIBS)

bin_PROGRAMS = dumer4 dumer6 dumer8 mmt4 mmt8 generic_dumer

BUILT_SOURCES = dumer4.c dumer6.c dumer8.c mmt4.c mmt8.c isd_cmdline.h cpucycles/cpucycles.c cpucycles/cpucycles.h
EXTRA_DIST = isd.ggo gen_dumer.py gen_mmt.py cpucycles/ workfactor myworkfactor.sh easy_challenge.txt

isd_cmdline.h: isd_cmdline.c
isd_cmdline.c: isd.ggo
	gengetopt -i isd.ggo -F isd_cmdline

dumer4.c: gen_dumer.py
	${PYTHON} gen_dumer.py 4 > dumer4.c
dumer6.c: gen_dumer.py
	${PYTHON} gen_dumer.py 6 > dumer6.c
dumer8.c: gen_dumer.py
	${PYTHON} gen_dumer.py 8 > dumer8.c

mmt4.c: gen_mmt.py
	${PYTHON} gen_mmt.py 4 > mmt4.c
mmt8.c: gen_mmt.py
	${PYTHON} gen_mmt.py 8 > mmt8.c

cpucycles/cpucycles.c cpucycles/cpucycles.h: 
	cd cpucycles/; sh do

COMMON_SRC_FILES = \
	main.c \
	isd.c isd.h \
	libisd.c libisd.h \
	ciht.c ciht.h \
	support.c support.h \
	final_test.c final_test.h \
	process_solutions.c process_solutions.h \
	sparse_words_list.c sparse_words_list.h \
	io.c io.h \
	prng.c prng.h \
	measure.c measure.h \
	custom_brilliantrussian.c custom_brilliantrussian.h \
	isd_cmdline.c isd_cmdline.h \
	cpucycles/cpucycles.c cpucycles/cpucycles.h \
	syscall_macros.h \
	sub_isd.h

generic_dumer_SOURCES = \
	${COMMON_SRC_FILES} \
	generic_dumer.c

dumer4_SOURCES = \
	${COMMON_SRC_FILES} \
	dumer4.c

dumer6_SOURCES = \
	${COMMON_SRC_FILES} \
	dumer6.c

dumer8_SOURCES = \
	${COMMON_SRC_FILES} \
	dumer8.c

mmt4_SOURCES = \
	${COMMON_SRC_FILES} \
	mmt4.c

mmt8_SOURCES = \
	${COMMON_SRC_FILES} \
	mmt8.c
