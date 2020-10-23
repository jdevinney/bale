#/******************************************************************
#
#
#  Copyright(C) 2020, Institute for Defense Analyses
#  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
# 
#
#  All rights reserved.
#  
#   This file is a part of Bale.  For license information see the
#   LICENSE file in the top level directory of the distribution.
#  
# 
# *****************************************************************/ 

#CFLAGS= -O3
#CFLAGS= -std=c99
LFLAGS= -lm
CC= gcc -Wall ${CFLAGS}

SRC= spmat_utils.h spmat_utils.c 
  SRC+=histo.c 
  SRC+=ig.c 
  SRC+=randperm.c 
  SRC+=toposort.c 
  SRC+=triangle.c 
  SRC+=unionfind.c 
  SRC+=transpose_matrix.c
  SRC+=permute_matrix.c
  SRC+=sssp.c
  SRC+=sssp_dijsktra.c
  SRC+=sssp_bellmanford.c
  SRC+=demo_spmat.c
  SRC+=opts_demo.c

EXE= histo 
  EXE+=ig 
  EXE+=randperm 
  EXE+=topo 
  EXE+=triangles 
  EXE+=unionfind
  EXE+=transpose_matrix
  EXE+=permute_matrix
  EXE+=sssp
  EXE+=demo_spmat
  EXE+=opts_demo

DOC= Doxyfile

OBJ= spmat_utils.o std_options.o

GRAPHS= demo.mm 10_ring.mm 2comps.mm

READMES= README.md histo.md ig.md randperm.md 
  READMES+= permute_matrix.md transpose_matrix.md
  READMES+= std_options.md spmat_utils.md
  READMES+= toposort.md triangle.md sssp.md unionfind.md
  READMES+= opts_demo.md

all: std_options.o spmat_utils.o ${EXE}

std_options.o: std_options.c std_options.h
	${CC} -g -c std_options.c

spmat_utils.o: spmat_utils.c spmat_utils.h
	${CC} -c spmat_utils.c -lm

demo_spmat: demo_spmat.c spmat_utils.o 
	${CC} demo_spmat.c spmat_utils.o ${LFLAGS} -o demo_spmat

opts_demo: opts_demo.c spmat_utils.o std_options.o
	${CC} opts_demo.c spmat_utils.o std_options.o ${LFLAGS} -o opts_demo

histo: histo.c std_options.o spmat_utils.o 
	${CC} -g histo.c std_options.o spmat_utils.o ${LFLAGS} ${LDFLAGS} -o histo

ig: ig.c std_options.o spmat_utils.o 
	${CC} ig.c std_options.o spmat_utils.o ${LFLAGS} ${LDFLAGS} -o ig

randperm: randperm.c std_options.o spmat_utils.o
	${CC} randperm.c std_options.o spmat_utils.o ${LFLAGS} ${LDFLAGS} -o randperm

topo: toposort.c std_options.o spmat_utils.o
	${CC} toposort.c std_options.o spmat_utils.o ${LFLAGS} ${LDFLAGS} -o topo

triangles: triangle.c std_options.o spmat_utils.o
	${CC} -g triangle.c std_options.o spmat_utils.o ${LFLAGS} ${LDFLAGS} -o triangles

unionfind: unionfind.c std_options.o spmat_utils.o
	${CC} unionfind.c std_options.o spmat_utils.o ${LFLAGS} ${LDFLAGS} -o unionfind

permute_matrix: permute_matrix.c std_options.o spmat_utils.o
	${CC} permute_matrix.c std_options.o spmat_utils.o ${LFLAGS} ${LDFLAGS} -o permute_matrix

transpose_matrix: transpose_matrix.c std_options.o spmat_utils.o
	${CC} -g transpose_matrix.c spmat_utils.o std_options.o ${LFLAGS} ${LDFLAGS} -o transpose_matrix

sssp: sssp.c sssp_dijsktra.c sssp_delta.c sssp_bellmanford.c std_options.o spmat_utils.o
	${CC} sssp.c sssp_dijsktra.c sssp_delta.c sssp_bellmanford.c std_options.o spmat_utils.o ${LFLAGS} ${LDFLAGS} -o sssp

doc: ${READMES}
	doxygen

test: std_options.o spmat_utils.o ${EXE}
	pytest -s

clean:
	rm -f *.out ${EXE} ${OBJ}
	rm -rf latex html
