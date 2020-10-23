## unionfind

### Definition

This is a directory for a future application called unionfind. Here
is a description of that application. It takes in a very sparse graph
(represented as a symmetrix adjacency matrix) and finds the components
of the graph. What counts as "finding the components" is up for
discussion, but it must be "easy" to query any vertex in the graph and
find out which component it lives in. After this initial computation,
waves of new edges are added to the graph and the components must be
updated. Whether these edges are added synchronously or asynchronously
is also up for discussion. We imagine it would be easier to consider
them coming in synchronous waves -- a batch of edges arrives, the
components are updated before the next batch arrives.

There is serial version of union-find in the "cousins" directory 
`other_serial/C`.



