# std_options

This library handles option parsing for the bale
[apps](../apps/README.md). This library is built on top of the arg_parse library in glibc. There are two main classes of options:
standard options, and standard graph options. Standard options are
included in all bale apps and give the user control of things like,
buffer count (for aggregation libraries), RNG seeds, and an implementation
mask (this controls which "implementations" are run for each app). The
standard graph options allow the user to control the input graph for
apps that require a matrix or graph. Run any bale app with --help to
see more information.

