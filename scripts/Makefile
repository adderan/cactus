rootPath = ../
include ../include.mk
sources = maftocigar.c
CFLAGS =-g -Wall -Werror
CC = gcc


all: ${sources}
	${CC} -o maftocigar ${CFLAGS} -I ${sonLibPath} ${sources} ${sonLibPath}/sonLib.a
	cp maftocigar ${binPath}
	
clean:
	cd ${binPath} && rm maftocigar
	rm maftocigar
	
