rootPath = ../
include ../include.mk
CFLAGS =-g -Wall -Werror
CC = gcc

all:maftocigar
maftocigar: alignmentFormats.o maftocigar.o
	${CC} ${CFLAGS} -o maftocigar alignmentFormats.o maftocigar.o ${sonLibPath}/sonLib.a
	cp maftocigar ${binPath}/maftocigar

maftocigar.o: maftocigar.c
	${CC} -c ${CFLAGS} -I ${sonLibPath} maftocigar.c
alignmentFormats.o: alignmentFormats.c alignmentFormats.h
	${CC} -c ${CFLAGS} -I ${sonLibPath} alignmentFormats.c
clean:
	cd ${binPath} && rm maftocigar
	rm maftocigar
	rm *.o
	
