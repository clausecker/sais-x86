CFLAGS = -O3 -Wall -Wno-parentheses -g
ASFLAGS= -g

.PHONY: all
all: suftest suftest.orig

suftest.orig: suftest.o sais.orig.o
	$(CC) $(LDFLAGS) -o suftest.orig suftest.o sais.orig.o $(LDLIBS)

suftest: suftest.o sais.o
	$(CC) $(LDFLAGS) -o suftest suftest.o sais.o $(LDLIBS)

clean:
	rm -f suftest suftest.orig sais.o sais.orig.o suftest.o

sais.o: sais.S sais-chr.S
sais.orig.o: sais.orig.c sais-chr.orig.c sais.h
suftest.o: sais.h
