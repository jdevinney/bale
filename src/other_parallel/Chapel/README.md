# bale (Chapel)

So far we have 5 different Chapel files in bale_chapel. None of these codes is using software aggregation libraries (like exstack or conveyors) because we haven't written them yet in Chapel. Hopefully that will be a future project! Our Chapel implementations mostly follow the AGP model as an experiment to see how clean they look in Chapel.

* histo.chpl :  This file implements histogram in a few different chapel styles.
* ig.chpl: This file implements index_gather in a few different chapel styles
* spmat.chpl: This file implements most of the sparse matrix functionality of the spmat library in bale_classic. We take advantage of Chapel iterators to make operations on sparse matrices much cleaner and more intuitive looking than the equivalent C code.
* topo.chpl: 2 different implementations of toposort in Chapel
* triangle.chpl: an implementation of triangle counting in Chapel