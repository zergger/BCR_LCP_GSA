CC = g++


FASTQ = 0  #Use fastQ file and handle quality score sequences
SAP = 0    #compute the reduced SAP array
RLO = 0    #compute the RLO-BWT

LCP = 0    #compute the LCP array
DA = 0     #compute the DA array

STORE_INDICES_DOLLARS = 0 #Store the indexes of each distint end-marker in eBWT string 
#For each end-marker, the binary file contains
#- indexes of the sequences
#- position in the eBWT string
#- symbol with which the associated suffix begins, i.e. the first symbol in the string.

PRINT = 0 #print in the terminal the BWT along with other data structures (if computed), such as LCP, DA, GSA, SAP array, and so on
#if STORE_INDICES_DOLLARS = 1 and PRINT = 1
# print the index of each end-markers symbol

IO_BUFFER_BYTES ?= 1048576
TRANSPOSE_BUFFER_BYTES ?= 134217728
STORE_LENGTHS ?= 0
COMPRESSED_CYC ?= 1
COMPRESSED_PARTIALS ?= 1
FINAL_BWT_FORMATS ?= 1

DEFINES = -DFASTQ=$(FASTQ) -DSAP=$(SAP) -DRLO=$(RLO) -DLCP=$(LCP) -DDA=$(DA) -DSTORE_ENDMARKER_POS=$(STORE_INDICES_DOLLARS) -DprintFinalOutput=$(PRINT) -DSIZEBUFFER=$(IO_BUFFER_BYTES) -DTRANSPOSE_BUFFER_BYTES=$(TRANSPOSE_BUFFER_BYTES) -DSTORE_LENGTH_IN_FILE=$(STORE_LENGTHS)
EBWT2INDEL_DEFINES = $(DEFINES) -DBCR_COMPRESSED_CYC=$(COMPRESSED_CYC) -DBCR_COMPRESSED_PARTIALS=$(COMPRESSED_PARTIALS) -DBCR_FINAL_BWT_FORMATS=$(FINAL_BWT_FORMATS)

CPPFLAGS = -Wall -ansi -pedantic -g -O3 -std=c++11 $(DEFINES)


BCR_BWTCollection_obs = BCR_BWTCollection.o BWTCollection.o BCRexternalBWT.o Tools.o Sorting.o TransposeFasta.o Timer.o -lz
BCR_sources = BCR_BWTCollection.cpp BWTCollection.cpp BCRexternalBWT.cpp Tools.cpp Sorting.cpp TransposeFasta.cpp Timer.cpp
BCR_ebwt2indel = build/ebwt2indel/BCR_LCP_GSA

BCR_BWTCollection: $(BCR_BWTCollection_obs)
	$(CC) -o BCR_LCP_GSA $(BCR_BWTCollection_obs)

.PHONY: ebwt2indel
ebwt2indel:
	mkdir -p build/ebwt2indel
	$(CC) -Wall -ansi -pedantic -g -O3 -std=c++11 $(EBWT2INDEL_DEFINES) -o $(BCR_ebwt2indel) $(BCR_sources) -lz

clean:
	rm -f core *.o *~ BCR_LCP_GSA

depend:
	$(CC) -MM *.cpp *.c > dependencies.mk

include dependencies.mk
