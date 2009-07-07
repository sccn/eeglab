CC=gcc
CFLAGS=-I. -Wall
LFLAGS=-lncurses
DEPS = ACQ_MessagePacket.h

all: AcqBuffer FieldtripBuffer

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

AcqBuffer: AcqBuffer.o
	gcc -o $@ $^ $(CFLAGS) $(LFLAGS)

FieldtripBuffer: FieldtripBuffer.o 
	gcc -o $@ $^ $(CFLAGS) $(LFLAGS)

