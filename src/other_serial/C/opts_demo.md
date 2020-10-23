## opt_demo

### Definition
Consistency among the command line options for the
bale apps was hard to maintain.  We are now using
argp to solve this problem.  This file contains four
complete programs (or the boilerplate for programs)
that are protect by #ifdef.  You have to un-comment
the appropriate #define and remake `opt_demo` to run
each of them.  The intent is to show that you can
modify default values for given options and add new
options without getting into the details of argp.

### References
