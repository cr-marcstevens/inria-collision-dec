CC = /usr/lib/gcc-snapshot/bin/gcc
CFLAGS = -Wall -Wextra -std=gnu99 -O3 -march=native -flto -g -mpopcnt -I$(M4RI)
LDFLAGS =
LDLIBS = -lpng -lm
M4RI = $(HOME)/Travaux/dev/m4ri

all: dumer4 dumer6 mmt4
#all: dumer4
 
debug: CFLAGS += -DDEBUG -g -O0 -mno-popcnt
debug: M4RI = $(HOME)/Travaux/dev/m4ri_dbg
debug: all

COMMON = main.c isd.c libisd.c final_test.c process_solutions.c sparse_words_list.c cpucycles/cpucycles.o io.c prng.c isd_cmdline.c measure.c $(M4RI)/.libs/libm4ri.a

dumer4: $(COMMON) dumer4.c htable.c
dumer6: $(COMMON) dumer6.c htable.c
mmt4: $(COMMON) mmt4.c iaht.c waht.c nocolht.c

profgen: CFLAGS += -fprofile-generate
profgen: all
profuse: CFLAGS += -fprofile-use
profuse: all

%_cmdline.c %_cmdline.h: %.ggo
	gengetopt -i $^ -F $*_cmdline

clean: 
	-rm dumer4 dumer6 mmt4
