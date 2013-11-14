AM_CPPFLAGS = -I$(M4RI_INCLUDE)

bin_PROGRAMS = dumer4

dumer4_SOURCES = \
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
	dumer4.c \
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
	cpucycles.h

dumer4_LDADD = $(M4RI_LIB)/libm4ri.a cpucycles/libcpucycles.a -lm
