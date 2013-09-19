CC = gcc
CFLAGS = -Wall -Wextra -std=gnu99 -I$(M4RI)
LDFLAGS =
LDLIBS = -lpng -lm $(M4RI)/.libs/libm4ri.a cpucycles/cpucycles.o
M4RI = m4ri
M4RI_DBG = m4ri_dbg

-include Makefile.local

.SUFFIXES:

default: CFLAGS += -O3 -march=native -g -mpopcnt
default: all
 
debug: CFLAGS += -DDEBUG -g -O0 -mno-popcnt
debug: M4RI = $(M4RI_DBG)
debug: all

stat: CFLAGS += -DSTAT
stat: all

profile: CFLAGS += -g -mno-popcnt
profile: LDLIBS += -lprofiler
profile: M4RI = $(M4RI_DBG)
profile: all

TARGET = dumer4 dumer6 dumer8 mmt4 mmt8 mmt12 BJMM compute_threshold dumer4_col dumer6_col dumer8_col
all: $(TARGET)

SOURCES = main.c isd.c libisd.c final_test.c process_solutions.c sparse_words_list.c io.c prng.c isd_cmdline.c measure.c custom_brilliantrussian.c

dumer%.c: gen_dumer.py
	python gen_dumer.py $* > $@

dumer%: $(SOURCES) dumer%.c counterht.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

dumer%_col: CFLAGS += -DMANAGE_COL
dumer%_col: $(SOURCES) dumer%.c counterht_col.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

mmt%.c: gen_mmt.py
	python gen_mmt.py $* > $@

mmt%: $(SOURCES) mmt%.c iaht.c waht.c nocolht.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

BJMM: $(SOURCES) BJMM.c BJMMtools.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

cpucycles/cpucycles.o: cpucycles/*.c cpucycles/*.h
	cd cpucycles; sh do; cd ..

%_cmdline.c %_cmdline.h: %.ggo
	gengetopt -i $^ -F $*_cmdline

compute_threshold: compute_threshold.c
	$(CC) $(CFLAGS) -o $@ $^

clean: 
	-rm -f $(TARGET)
