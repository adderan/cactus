rootPath = ../
include ../include.mk

#host=kolossus-10
#batchSystem = parasol 
#parasolCommand=/cluster/home/markd/pub/parasol -host=swarm-10
#tempDir=/scratch/benedict/
#maxThreads=30

host=localhost
batchSystem = singleMachine
parasolCommand=parasol
tempDir=./
maxCores=4

port=1978
databaseOptions=-ls -tout 200000 -th 64 -bgsi 100000000 #opts=ls#bnum=30m#msiz=50g#ktopts=p
minRecordSize=0
maxRecordSize=100
keysPerJob=10001
totalJobs=30

toil=${tempDir}/toil
log = ${tempDir}/log.txt
databaseDir=${tempDir}/testDb

all : ${binPath}/dbTestScript 

${binPath}/dbTestScript  : ${basicLibsDependencies} *.c
	${cxx} ${cflags} -I inc -I${libPath} -Wno-error -o ${binPath}/dbTestScript *.c ${basicLibs}

clean :
	rm -rf ${binPath}/dbTestScript

test :
	ktserver -log ${log} -host ${host} -port ${port} ${databaseOptions} &
	rm -rf ${databaseDir} ${jobTree}
	python dbTestScript.py --maxCores=${maxCores} --parasolCommand="${parasolCommand}" --batchSystem ${batchSystem} --logLevel DEBUG --jobStore ${toil} --stats --keysPerJob ${keysPerJob} --totalJobs ${totalJobs} --host ${host} --port ${port} --databaseDir ${databaseDir} --databaseOption "${databaseOptions}" --minRecordSize ${minRecordSize} --maxRecordSize ${maxRecordSize}
	toil stats --jobStore ${toil} --outputFile ./jobTreeStatsTest.xml
	ktremotemgr report -host ${host} -port ${port}
	rm -rf ${databaseDir} ${toil} ${log}
	ps ax | grep 'ktserver' | cut -f1 -d' ' | xargs kill
