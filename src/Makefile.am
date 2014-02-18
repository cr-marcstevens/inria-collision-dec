AM_CPPFLAGS = -I$(M4RI_INCLUDE)

bin_PROGRAMS = dumer4 dumer6 dumer8 mmt4 mmt8 generic_dumer

BUILT_SOURCES = dumer4.c dumer6.c dumer8.c mmt4.c mmt8.c isd_cmdline.h
EXTRA_DIST = isd.ggo gen_dumer.py gen_mmt.py

generic_dumer_SOURCES = \
	generic_dumer.c \
	main.c \
	isd.c \
	libisd.c \
	final_test.c \
	process_solutions.c \
	sparse_words_list.c \
	io.c \
	prng.c \
	isd_cmdline.c \
	measure.c \
	custom_brilliantrussian.c \
	counterht.c \
	custom_brilliantrussian.h \
	isd_cmdline.h \
	libisd.h \
	prng.h \
	sparse_words_list.h \
	syscall_macros.h \
	counterht.h \
	final_test.h \
	io.h \
	isd.h \
	measure.h \
	process_solutions.h \
	sub_isd.h \
	cpucycles/cpucycles.h



isd_cmdline.h: isd_cmdline.c
isd_cmdline.c: isd.ggo
	gengetopt -i isd.ggo -F isd_cmdline

dumer4.c: gen_dumer.py
	python gen_dumer.py 4 > dumer4.c
dumer6.c: gen_dumer.py
	python gen_dumer.py 6 > dumer6.c
dumer8.c: gen_dumer.py
	python gen_dumer.py 8 > dumer8.c

mmt4.c: gen_mmt.py
	python gen_mmt.py 4 > mmt4.c
mmt8.c: gen_mmt.py
	python gen_mmt.py 8 > mmt8.c

dumer4_SOURCES = \
	dumer4.c \
	main.c \
	isd.c \
	libisd.c \
	final_test.c \
	process_solutions.c \
	sparse_words_list.c \
	io.c \
	prng.c \
	isd_cmdline.c \
	measure.c \
	custom_brilliantrussian.c \
	counterht.c \
	custom_brilliantrussian.h \
	isd_cmdline.h \
	libisd.h \
	prng.h \
	sparse_words_list.h \
	syscall_macros.h \
	counterht.h \
	final_test.h \
	io.h \
	isd.h \
	measure.h \
	process_solutions.h \
	sub_isd.h \
	cpucycles/cpucycles.h

dumer6_SOURCES = \
	dumer6.c \
	main.c \
	isd.c \
	libisd.c \
	final_test.c \
	process_solutions.c \
	sparse_words_list.c \
	io.c \
	prng.c \
	isd_cmdline.c \
	measure.c \
	custom_brilliantrussian.c \
	counterht.c \
	custom_brilliantrussian.h \
	isd_cmdline.h \
	libisd.h \
	prng.h \
	sparse_words_list.h \
	syscall_macros.h \
	counterht.h \
	final_test.h \
	io.h \
	isd.h \
	measure.h \
	process_solutions.h \
	sub_isd.h \
	cpucycles/cpucycles.h

dumer8_SOURCES = \
	dumer8.c \
	main.c \
	isd.c \
	libisd.c \
	final_test.c \
	process_solutions.c \
	sparse_words_list.c \
	io.c \
	prng.c \
	isd_cmdline.c \
	measure.c \
	custom_brilliantrussian.c \
	counterht.c \
	custom_brilliantrussian.h \
	isd_cmdline.h \
	libisd.h \
	prng.h \
	sparse_words_list.h \
	syscall_macros.h \
	counterht.h \
	final_test.h \
	io.h \
	isd.h \
	measure.h \
	process_solutions.h \
	sub_isd.h \
	cpucycles/cpucycles.h

mmt4_SOURCES = \
	mmt4.c \
	main.c \
	isd.c \
	libisd.c \
	final_test.c \
	process_solutions.c \
	sparse_words_list.c \
	io.c \
	prng.c \
	isd_cmdline.c \
	measure.c \
	custom_brilliantrussian.c \
	iaht.c \
	waht.c \
	nocolht.c \
	custom_brilliantrussian.h \
	isd_cmdline.h \
	libisd.h \
	prng.h \
	sparse_words_list.h \
	syscall_macros.h \
	iaht.h \
	waht.h \
	nocolht.h \
	final_test.h \
	io.h \
	isd.h \
	measure.h \
	process_solutions.h \
	sub_isd.h \
	cpucycles/cpucycles.h


mmt8_SOURCES = \
	mmt8.c \
	main.c \
	isd.c \
	libisd.c \
	final_test.c \
	process_solutions.c \
	sparse_words_list.c \
	io.c \
	prng.c \
	isd_cmdline.c \
	measure.c \
	custom_brilliantrussian.c \
	iaht.c \
	waht.c \
	nocolht.c \
	custom_brilliantrussian.h \
	isd_cmdline.h \
	libisd.h \
	prng.h \
	sparse_words_list.h \
	syscall_macros.h \
	iaht.h \
	waht.h \
	nocolht.h \
	final_test.h \
	io.h \
	isd.h \
	measure.h \
	process_solutions.h \
	sub_isd.h \
	cpucycles/cpucycles.h

dumer4_LDADD = $(M4RI_LIB)/libm4ri.a cpucycles/libcpucycles.a -lm
dumer6_LDADD = $(M4RI_LIB)/libm4ri.a cpucycles/libcpucycles.a -lm
dumer8_LDADD = $(M4RI_LIB)/libm4ri.a cpucycles/libcpucycles.a -lm
mmt4_LDADD = $(M4RI_LIB)/libm4ri.a cpucycles/libcpucycles.a -lm
mmt8_LDADD = $(M4RI_LIB)/libm4ri.a cpucycles/libcpucycles.a -lm
generic_dumer_LDADD = $(M4RI_LIB)/libm4ri.a cpucycles/libcpucycles.a -lm
