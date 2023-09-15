#/*****************************************************************************
#   Copyright (C) 1993-2008 by Phil Green.                          
#   All rights reserved.                           
#                                                                           
#   This software is part of a beta-test version of the swat/cross_match/phrap 
#   package.  It should not be redistributed or
#   used for any commercial purpose, including commercially funded
#   sequencing, without written permission from the author and the
#   University of Washington.
#   
#   This software is provided ``AS IS'' and any express or implied
#   warranties, including, but not limited to, the implied warranties of
#   merchantability and fitness for a particular purpose, are disclaimed.
#   In no event shall the authors or the University of Washington be
#   liable for any direct, indirect, incidental, special, exemplary, or
#   consequential damages (including, but not limited to, procurement of
#   substitute goods or services; loss of use, data, or profits; or
#   business interruption) however caused and on any theory of liability,
#   whether in contract, strict liability, or tort (including negligence
#   or otherwise) arising in any way out of the use of this software, even
#   if advised of the possibility of such damage.
#                                                                         
#*****************************************************************************/

# Makefile for swat, phrap, cross_match, phrapview; also swprobs (currently inactivated)
#CC= icc
#CFLAGS= -O2 -wd266,880
CC= cc
CFLAGS= -O2
LFLAGS= -lm

SWATOBJS= swat.o weibull.o

ALPHOBJS= alpha_swat.o weibull.o fast_smith_wat.o

COUNOBJS= swat_counts.o weibull.o smith_wat_counts.o

PHOBJS= contigs.o tig_node.o nodes.o

CROBJS= 

CLOBJS= cluster.o

LOCOOBJS= loco.o 

SWOBJS= alignments.o db.o smith_wat.o full_smith_wat.o quick_smith_wat.o utilities.o nw.o full_nw.o profile.o parameters.o 

SCANOBJS= readin.o words.o segments.o recursive_swat.o pairs.o cand_pairs.o diffs.o names.o anomalies.o qual.o tags.o analyze_splices.o

ERFOBJS= erf_test.o

PROBOBJS= probs.o profile.o utilities.o

JKBSWOBJS= exp_file.o seqInfo.o expFileIO.o array.o compress.o scf_bits.o

CALFOBJS= calf_merge.o  calf_utils.o

all:  swat phrap cross_match phrapview cluster loco calf_merge

#manyreads:  
#	touch swat.h;
#	make CFLAGS="$(CFLAGS) -DMANYREADS" phrap cross_match gcphrap; was prev commented out!!
#	make CFLAGS="$(CFLAGS) -DMANYREADS" phrap cross_match;
#	mv phrap phrap.manyreads;
#	mv cross_match cross_match.manyreads;
#	touch swat.h;
#	make CFLAGS="$(CFLAGS) -DLONGREADS" phrap;
#	mv phrap phrap.longreads;
#	touch swat.h;
#	make phrap cross_match;

#all:
#	make CFLAGS="-O2 -DMANYREADS" phrap cross_match
#	mv phrap phrap.manyreads
#	mv cross_match cross_match.manyreads
#	make clean
#	make CFLAGS="-O2 -DLONGREADS" phrap
#	mv phrap phrap.longreads
#	make clean
#	make CFLAGS="-O2"  swat phrap cross_match phrapview cluster loco 


clean:
	rm $(SWATOBJS) $(SWOBJS) $(PHOBJS) $(SCANOBJS) $(CROBJS) cross_match.o gccross_match.o phrap.o gcphrap.o call_subs.o gccall_subs.o cluster.o loco.o

tar:
	tar cvf distrib.tar makefile *.h *.c phrap.doc general.doc swat.doc PAM250 BLOSUM62 BLOSUM50 mat50 mat70 mb_matrix penalty2 vector.seq phrapview
#	compress distrib.tar N.B. INSTRUCTIONS NEED TO BE CHANGED!!
	gzip distrib.tar
	uuencode distrib.tar.gz distrib.tar.gz > distrib.uu

tar.other:
	tar cvf other.tar NOTES *.perl 
	tar cvf swat.misc.tar erf_test.c z_dist.c smith_wat.c.old smith_wat.c.sparse 
	compress swat.misc.tar other.tar

calf_merge: makefile $(CALFOBJS) 
	$(CC) $(CFLAGS) -o $@ $(CALFOBJS) $(LFLAGS)

phrapview: makefile
	chmod a+x phrapview

swat: makefile $(SWATOBJS) $(SWOBJS)
	$(CC) $(CFLAGS) -o $@ $(SWATOBJS) $(SWOBJS) $(LFLAGS)
	chmod o-r swat

alpha_swat: makefile $(ALPHOBJS) $(SWOBJS) 
	$(CC) $(CFLAGS) -o $@ $(ALPHOBJS) $(SWOBJS) $(LFLAGS) 
	chmod o-r alpha_swat

swat_counts: makefile $(COUNOBJS) $(SWOBJS) 
	$(CC) $(CFLAGS) -o $@ $(COUNOBJS) $(SWOBJS) $(LFLAGS) 
	chmod o-r swat_counts

swprobs: makefile $(PROBOBJS)
	$(CC) $(CFLAGS) -o $@ $(PROBOBJS) $(LFLAGS)

phrap: makefile phrap.o call_subs.o $(PHOBJS) $(SCANOBJS) $(SWOBJS)
	$(CC) $(CFLAGS) -o $@ phrap.o call_subs.o $(PHOBJS) $(SCANOBJS) $(SWOBJS) $(LFLAGS)
	chmod o-r phrap

gcphrap: makefile gcphrap.o gccall_subs.o $(PHOBJS) $(SCANOBJS) $(SWOBJS) $(JKBSWOBJS)
	$(CC) $(CFLAGS) -o $@ gcphrap.o gccall_subs.o $(PHOBJS) $(SCANOBJS) $(SWOBJS) $(JKBSWOBJS) $(LFLAGS)
	chmod o-r gcphrap

cross_match: makefile cross_match.o call_subs.o $(CROBJS) $(SCANOBJS) $(SWOBJS)
	$(CC) $(CFLAGS) -o $@ cross_match.o call_subs.o $(CROBJS) $(SCANOBJS) $(SWOBJS) $(LFLAGS)
	chmod o-r cross_match

gccross_match: makefile gccross_match.o gccall_subs.o $(CROBJS) $(SCANOBJS) $(SWOBJS) $(JKBSWOBJS)
	$(CC) $(CFLAGS) -o $@ gccross_match.o gccall_subs.o $(CROBJS) $(SCANOBJS) $(SWOBJS) $(JKBSWOBJS) $(LFLAGS)
	chmod o-r gccross_match

cluster: makefile $(CLOBJS) call_subs.o $(SCANOBJS) $(SWOBJS)
	$(CC) $(CFLAGS) -o $@ $(CLOBJS) call_subs.o $(SCANOBJS) $(SWOBJS) $(LFLAGS)
	chmod o-r cluster

loco: makefile $(LOCOOBJS) $(SWOBJS)
	$(CC) $(CFLAGS) -o $@ $(LOCOOBJS) $(SWOBJS) $(LFLAGS)
	chmod o-r loco

erf_test: $(ERFOBJS)
	$(CC) $(CFLAGS) -o $@ $(ERFOBJS) $(LFLAGS)

calf_merge.o: calf.h

calf_utils.o: calf.h

loco.o: swat.h

cluster.o: swat.h

alignments.o: swat.h

readin.o: swat.h 

db.o: swat.h

phrap.o: swat.h 

gcphrap.o: phrap.c swat.h 
	$(CC) $(CFLAGS) -c -o $@  phrap.c -DGCPHRAP

call_subs.o: swat.h

gccall_subs.o: call_subs.c swat.h
	$(CC) $(CFLAGS) -c -o $@  call_subs.c -DGCPHRAP

cross_match.o: swat.h 

gccross_match.o: cross_match.c swat.h 
	$(CC) $(CFLAGS) -c -o $@  cross_match.c -DGCPHRAP

smith_wat.o: swat.h

full_smith_wat.o: smith_wat.c swat.h
	$(CC) $(CFLAGS) -c -o $@  smith_wat.c -DFINDALIGN

quick_smith_wat.o: smith_wat.c swat.h
	$(CC) $(CFLAGS) -c -o $@  smith_wat.c -DQUICKALIGN

smith_wat_counts.o: smith_wat.c swat.h
	$(CC) $(CFLAGS) -c -o $@  smith_wat.c -DCOUNTS

swat.o: swat.h

weibull.o: swat.h

nw.o: swat.h

probs.o: swat.h

profile.o: swat.h

full_nw.o: nw.c swat.h
	$(CC) $(CFLAGS) -c -o $@  nw.c -DFINDALIGN

fast_smith_wat.o: fast_smith_wat.s 
	$(CC) $(CFLAGS) -c -o $@ fast_smith_wat.s 

alpha_swat.o: swat.c swat.h
	$(CC) $(CFLAGS) -c -o $@ swat.c -DALPHA

swat_counts.o: swat.c swat.h
	$(CC) $(CFLAGS) -c -o $@ swat.c -DCOUNTS

contigs.o: swat.h 

segments.o: swat.h 

parameters.o: swat.h 

words.o: swat.h 

recursive_swat.o: swat.h 

names.o: swat.h 

nodes.o: swat.h 

diffs.o: swat.h 

pairs.o: swat.h 

analyze_splices.o: swat.h 

cand_pairs.o: swat.h 

qual.o: swat.h 

tags.o: swat.h 

utilities.o: 

anomalies.o: swat.h 

fast_smith_wat.s: fast_smith_wat.c
	cc -S -machine_code -source_listing -O3 fast_smith_wat.c
#N.B. fast_smith_wat.s MUST BE EDITTED BEFORE USE!!!

tig_node.o: swat.h
