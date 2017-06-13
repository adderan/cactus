rootPath = ../
include ../include.mk

cflags += ${tokyoCabinetIncl}

all : ${binPath}/cactus_adjustSamplingRates

${binPath}/cactus_adjustSamplingRates : *.c ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_adjustSamplingRates cactus_adjustSamplingRates.c ${libPath}/cactusLib.a ${basicLibs}
