all: histo ig topo spmat triangle

histo: histo.chpl
	chpl -o histo histo.chpl

ig: ig.chpl
	chpl -o ig ig.chpl

topo: topo.chpl
	chpl -o topo topo.chpl

spmat: spmat.chpl
	chpl -o spmat spmat.chpl

triangle: triangle.chpl
	chpl -o triangle triangle.chpl

clean:	
	rm -f	histo      ig      topo      spmat      triangle      \
		histo_real ig_real topo_real spmat_real triangle_real
