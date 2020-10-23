# Triangle counting

## Definition

This uses matrix algebra approach to counting triangles in a graph.
See the book, "Graph Algorithms in the Language of Linear Algebra",
edited by Gilbert, and Kepner for more details on our approach to this problem.

## Discussion

We have implemented two algorithms to count triangles in a graph. The
first computes (L & L * U) and the second computes (L & U * L), where
'&' means element-wise AND and '*' is ordinary matrix
multiplication. 

#### Why is it in bale?

Our implementation of triangle counting has several interesting properties. For one thing, the algorithm allows an owner of a row to "push" it's nonzeros to other rows that require that data and have them perform the computation, or to "pull" (or get with PGAS reads) the data itself and do the computation locally. This freedom is not found in other bale apps. The "push" model maps onto aggregation quite nicely while the "pull" model is much more natural using the AGP model. Further, the two different models, "push" and "pull", can have drastic consequences on communication volume depending on the algorithm used (either (L & L * U) or (L & U * L)).  Finally, as others have shown in publications, this application is interesting because there are ways to reduce communication by permuting the rows and columns of the input matrix if one wanted to speed up the calculation.

#### From the Book?

The pull version of the AGP code is quite simple, but suffers from the same inscrutability that our other AGP matrix codes do. The push versions of the aggregated codes are also pretty simple, but could be improved with more modern language tools (like iterators for the sparse matrix data structure).