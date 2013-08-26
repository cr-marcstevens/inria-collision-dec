include Makefile.local
CC ?= gcc
CFLAGS += -Wall -Wextra -std=gnu99 -I$(M4RI)
LDFLAGS +=
LDLIBS += -lpng -lm
M4RI ?= $(HOME)/Travaux/dev/m4ri
M4RIG ?= $(HOME)/Travaux/dev/m4ri_dbg

TARGET = dumer4 dumer6 dumer8 mmt4 mmt8 mmt12 BJMM compute_threshold

default:

all: $(TARGET)

default: CFLAGS += -O3 -march=native -g -mpopcnt
default: all
 
debug: CFLAGS += -DDEBUG -g -O0 -mno-popcnt
debug: M4RI = $(M4RIG)
debug: all

stat: CFLAGS += -DSTAT
stat: all

profile: CFLAGS += -g -mno-popcnt
profile: LDLIBS += -lprofiler
profile: M4RI = $(M4RIG)
profile: all

COMMON = main.c isd.c libisd.c final_test.c process_solutions.c sparse_words_list.c cpucycles/cpucycles.o io.c prng.c isd_cmdline.c measure.c $(M4RI)/.libs/libm4ri.a custom_brilliantrussian.c

dumer%.c: gen_dumer.py
	python gen_dumer.py $* > $@

dumer4: $(COMMON) dumer4.c counterht_nocol.c
dumer6: $(COMMON) dumer6.c counterht_nocol.c
dumer8: $(COMMON) dumer8.c counterht_nocol.c


mmt%.c: gen_mmt.py
	python gen_mmt.py $* > $@

mmt4: $(COMMON) mmt4.c iaht.c waht.c nocolht.c
mmt8: $(COMMON) mmt8.c iaht.c waht.c nocolht.c
mmt12: $(COMMON) mmt12.c iaht.c waht.c nocolht.c

BJMM: $(COMMON) BJMM.c BJMMtools.c

compute_threshold: compute_threshold.c

profgen: CFLAGS += -fprofile-generate
profgen: all
profuse: CFLAGS += -fprofile-use
profuse: all

test: $(TARGET)
	for i in $(TARGET); do echo $$i; ./$$i -i small_challenge.txt -l 15 -I 1; done

%_cmdline.c %_cmdline.h: %.ggo
	gengetopt -i $^ -F $*_cmdline

clean: 
	-rm -f $(TARGET)
