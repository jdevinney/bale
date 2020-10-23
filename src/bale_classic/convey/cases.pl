#!/usr/bin/perl --

=pod

=for html <style type="text/css">h1{font-size:150%;color:#990000;}
h2{font-size:120%;color:#990000;}</style>

=for html <div style="text-align: center; font-size: 120%; color: #006600; margin-bottom: 1ex;">
<strong>UNCLASSIFIED//FOR OFFICIAL USE ONLY</strong></div>

=head1 NAME

cases - a source code preprocessor

=head1 SYNOPSIS

  perl cases.pl [options] [infile] > outfile

=head1 DESCRIPTION

The B<cases> script is a simple C/C++ code preprocessor that can
replicate a block of code while specializing one or more variables in
each copy, so that the compiler can optimize each special case
separately.  It recognizes constructions that look similar to directives
for the standard C preprocessor.  It reads from the given file, if
present, and otherwise from standard input, and it writes to standard
output.

Only complete tokens (preceded and followed by non-word characters) are
specialized.  B<Cases> is oblivious of source code syntax, however: it
makes substitutions even within literal strings and comments.  If you
have a string constant in which you do not want substitutions to occur,
you should move it outside the relevant B<cases> directives, declaring a
new variable to hold its value if necessary.

=head2 Directive Syntax

Each directive occupies either a single line of the input file, or
multiple lines joined by continuation marks.  A directive (except #encase)
is continued onto the next line if it ends with a backslash, optionally
followed by whitespace.  Comments are allowed at the end of directive
lines and are stripped off before looking for continuation marks.  By
default, comments are introduced by the double slash C<//>, but this
string can by changed via the B<--comment> option.

Whitespace is permitted, and ignored, before and after each punctuation
character in a directive (equal sign, parentheses, commas).  Whitespace
may also appear between '#' and the directive name (C<cases>,
C<endcases>, etc.).  Lettercase is not significant in directive names.

All directives except #defcase and #encase must occur in matched pairs,
and the pairs must be properly nested.  The following directives are
recognized.

=over

=item B<#calc> ... B<#endcalc>

A construction of the form

    #calc <name> = <value>
    (lines of code)
    #endcalc

where <value> is an integer-valued expression, is processed by evaluating
the expression and replacing the given <name> by the decimal representation
of the resulting integer.  It is equivalent to

    #loop <name> = <value>, <value>
    (lines of code)
    #endloop

A construction of the form

    #calc (<name1>,...,<namek>) = (<value1>,...,<valuek>)
    (lines of code)
    #endcalc

evaluates each of the given <value>s and substitutes the results for the
given <name>s.  It is equivalent to k nested #calc constructions.

Expressions may contain integer constants, parentheses, and all the
C operators that operate on integers (see EXPRESSIONS below).  The
effect of any other expression is undefined.  An expression I<cannot>
contain identifiers whose values are given by #define directives.

=item B<#cases> ... B<#endcases>

A construction of the form

    #cases <name> = <value1>, <value2>, ..., <valuen>
    (lines of code)
    #endcases

is expanded to code like this:

    if (<name> == (<value1>)) {
    (lines of code, with the token <name> replaced by <value1>)
    } else
    if (<name> == (<value2>)) {
    (lines of code, with the token <name> replaced by <value2>)
    } else
    ...
    if (<name> == (<valuen>)) {
    (lines of code, with the token <name> replaced by <valuen>)
    }

The last <value> can be simply <name>, in which case the final
vacuous test if (<name> == <name>) is omitted and the final block of
generated code is a copy of the input code.  This enables the code
to handle the general case as well as n-1 special cases.

The #cases construction can also specialize several variables at once,
using the following syntax:

    #cases (<name1>,...,<namek>) = (<v11>,...,<v1k>), ..., (<vn1>,...,<vnk>)
    ...
    #endcases

Each case can be specific, partially general, or wholly general.  For
example, the construction

    #cases (foo,bar) = (1,1), (foo,1), (foo,bar)
    (lines of code)
    #endcases

expands to

    if (foo == (1) && bar == (1)) {
    (lines of code, with both foo and bar replaced by 1)
    } else
    if (bar == (1)) {
    (lines of code, with bar replaced by 1)
    } else
    {
    (lines of code, unchanged)
    }

The fully general case, if present, must be the last case in the list.

=item B<#defcases> ... B<#endcases>

This directive replicates the intervening lines without adding any
surrounding material like an I<if> statement.  Its syntax is identical to
that of #cases: it can substitute values for either a single name or
multiple names.  The primary purpose of this directive is to create
multiple declarations or definitions.

=item B<#defcase> 

Acts like #defcases, except that only one "case" must be present, and no
matching #endcases directive is required or allowed.  The scope of this
directive extends from its location to the end of the input.  Directives
must still be properly nested, so #defcase cannot appear within any
matched pair of other directives.  (The script converts #defcase to
#defcases and adds a matching #endcases at the end of the input.)

=item B<#edit> ... B<#endedit>

A construction of the form

    #edit <name> = <value>
    (lines of code)
    #endedit

where <value> is a string-valued expression, is processed by evaluating
the expression and replacing the given <name> by the resulting string.
A construction of the form

    #edit (<name1>,...,<namek>) = (<value1>,...,<valuek>)
    (lines of code)
    #endedit

evaluates each of the given <value>s and substitutes the results for the
given <name>s.  It is equivalent to k nested #edit constructions.

Expressions may contain string constants enclosed in single quotes,
the concatenation (dot) operator, parentheses, decimal numbers,
and the Perl operators chr, lc, lcfirst, substr, uc, ucfirst, x.  In the
simple form with only one <name>, commas are also allowed.  An expression
I<cannot> contain identifiers whose values are given by #define directives.
To use an identifier substituted by an enclosing #defcases or similar
directive, enclose it in single quotes.

=item B<#encase "filename">

Replace this directive by the contents of the given file, and process
#encase directives in that file recursively.  This directive is
analogous to #include in the C preprocessor.  If the filename does not
begin with a slash, B<cases> looks for the file in the current directory
followed by any directories provided via B<-E> or B<--encase> options,
in the order that they appear on the command line.

=item B<#hexcalc> ... B<#endcalc>

Equivalent to #calc ... #endcalc except that the expressions, after
evaluation, are converted to uppercase hexadecimal (rather than decimal)
digit strings, with no leading C<0x>.

=item B<#hexloop> ... B<#endloop>

See #loop ... #endloop below.  The #hexloop directive differs from #loop
only in that the integer values are converted to uppercase hexadecimal
(rather than decimal) digit strings, and these strings are of fixed
length: they contain just enough digits to represent the upper limit of
the iteration.

=item B<#ifcases> ... B<#endcases>

A synonym for #cases ... #endcases.

=item B<#loop> ... B<#endloop>

A construction of the form

    #loop <name> = <start>, <end>, <step>
    (lines of code)
    #endloop

where <start>, <end>, and <step> are integer-valued expressions, expands
to floor((<end>-<start>)/<step>)+1 consecutive copies of the lines of
code, with <name> replaced by the decimal representations of <start>,
<start>+<step>, <start>+2*<step>, and so on, in succession.  The
replication count can be zero: If <step> is positive and <end>-<start>
is negative, or vice versa, then the surrounded lines of code are
suppressed entirely.  The <step> expression and the preceding comma may
be omitted, in which case <step> defaults to 1.  Both <start> and <end>
must evaluate to nonnegative integers, and <step>, if present, must
evaluate to a nonzero integer.

=item B<#swcases> ... B<#endcases>

A construction of the form

    #swcases <name> = <value1>, <value2>, ..., <valuen>
    (lines of code)
    #endcases

is expanded to code like this:

    case <value1>: {
    (lines of code, with the token <name> replaced by <value1>)
    } break;
    case <value2>: {
    (lines of code, with the token <name> replaced by <value2>)
    } break;
    ...
    case <value2>: {
    (lines of code, with the token <name> replaced by <valuen>)
    } break;

If one of the <value>s is <name> itself, the corresponding block
of output code has the form

    default: {
    (lines of code)
    } break;

No containing I<switch> statement is generated.

The #swcases ... #endcases construction cannot specialize more
than one variable at once.

=back

=head2 General Features

The #encase directives are processed recursively before any others.
Therefore no substitutions can be made in the names of files to be
incorporated.

Whenever B<cases> substitutes a value for a name, it removes any
occurrence of ## that immediately precedes or follows the name without
intervening whitespace.  This feature, which lets B<cases> construct new
tokens, is commonly used within #defcases .. #endcases to build function
names.  Note that if two names are being substituted with no intervening
characters, they must be separated by ####, because each substitution
removes ##.

The #*cases ... #endcases constructs can be nested.  The order in which
these constructs are expanded is undefined; if you are so devious as to
depend on it, you deserve what you get.  Likewise, when a #*cases or
#*calc directive involves multiple identifiers, the order in which the
identifiers are substituted is undefined.

The #*loop ... #endloop, #*calc ... #endcalc, and #edit ... #endedit
constructs can also be nested and interspersed with #*cases and #endcases
directives.  All #*cases directives are processed before any #*loop,
#*calc, or #edit directive is processed, and then the #*loop, #*calc, and
#edit directives are processed from outermost to innermost.  Consequently,
loop limits can involve names that are specialized by containing
#*loop ... #endloop and #*cases ... #endcases blocks.

=head2 Expression Syntax

At the time it is evaluated, an integer-valued expression may contain
only integer constants (decimal and/or hexadecimal), parentheses, and
the following operators:

=over

=item *

the arithmetic operators +, -, *, /, %;

=item *

the bitwise operators <<, >>, &, |, ^, ~;

=item *

the relational operators <, <=, ==, !=, >, >=, and the
boolean operators !, &&, and ||, all of which return 0 or 1; and

=item *

the ternary conditional operator ?:.

=back

The precedence of operators is the same as in C.


=head1 OPTIONS

=over

=item B<-a>, B<--arithmetic>

After performing all other substitutions, replace each occurrence of
##[[<expr>]] with the evaluation of <expr>, written in uppercase
hexadecimal, and then replace each occurrence of ##[<expr>] with the
evaluation of <expr>, written in decimal.

=item B<-c>[=]I<string>, B<--comment>[=]I<string>

Sets the I<string> used to identify the beginning of a comment
in directives recognized by B<cases>.  The same string is used
by B<--number> to introduce line-number comments.  Its default
value is C<//>.

=item B<-d> <name>=<value>, B<--defcase> <name>=<value>

Prepends the directive #defcase <name>=<value> to the input,
unless a #defcase with the same <name> is already present.
In the latter case, the <value> on the command line overrides
the <value> in the input.

=item B<-E>[=]I<directory>, B<--encase>[=]I<directory>

Add I<directory> to the search path for files named in #encase
directives.

=item B<-f>, B<--file>

Replaces each occurrence of the token C<__FILE__> by the name of the input
file in which it appears.  At the top level this name is given by the
command-line argument, which must be present.  Within an "encased" file,
the name is the one given in the #encase directive.

=item B<-h>, B<--header>

After processing #encase directives, simply convert each #defcase
directive to a corresponding #define directive for the C preprocessor,
output those #define directives, and stop.

=item B<-l>, B<--line>

Replaces each occurrence of the token C<__LINE__> by the number of the
source line on which it appears.  (The first line is line number 1.)
For C code, the B<-p> option is probably more useful.

=item B<-n>, B<--number>

Before transforming the input, appends to each input line (except for
lines recognized as directives by B<cases>) a comment that records the
original line number.  The line number is surrounded by equal signs.
For C code, the B<-p> option is probably more useful.

=item B<-p>, B<--preserve>

Inserts #line directives into the output so that the C compiler and
follow-on tools will know the line numbers in the original source file.
The source file must be provided as an argument rather than through
standard input.

=item B<-v>, B<--version>

Print the version number and exit.

=back


=head1 EXAMPLES

=over

=item 1.

Suppose we have an n-long array a[] of 64w-bit items, where w is unknown
at compile time but probably very small, and we want to shift each item
left by k bits, where 0 < k < 64.  The following code should be
efficient:

    #cases w = 1, 2, 3, w
      for (i = 0; i < n; i++) {
        for (j = 0; j < w-1; j++)
          a[w*i+j] = (a[w*i+j] << k) | (a[w*i+j+1] >> 64-k);
        a[w*i+w-1] <<= k;
      }
    #endcases

For w = 1, 2, and 3, the compiler should fully unroll the loop over j
and then simplify the body of the loop over i.  One could write this
code in a more verbose style with explicit temporary variables, but it
shouldn't make any difference for the cases of fixed w.  (It might
improve the general case, as the compiler probably doesn't know that w
is positive.)

=item 2.

The following code defines functions for performing insertion sort
on half-word (i.e., 32-bit), one-word (64-bit), two-word, and three-word
items, respectively.  The sort key is contained in the first word of
each item and is specified by a mask.

    typedef struct { uint32_t key; } halfword;
    typedef struct { uint64_t key; } oneword;
    typedef struct { uint64_t key; uint64_t data; } twoword;
    typedef struct { uint64_t key; uint64_t data[2]; } threeword;

    #defcases Type = halfword, oneword, twoword, threeword
    void isort_##Type(Type *arr, int64_t n, uint64_t mask)
    {
      int64_t i, j;
      Type ins;
      for (i = 1; i < n; i++) {
        ins = arr[i];
        for (j = i; j > 0 && (arr[j-1].key & mask) >
                    (ins.key & mask)); j--)
          arr[j] = arr[j-1];
        arr[j] = ins;
      }
    }
    #endcases

=back


=begin comment

=head1 DIAGNOSTICS

[All error messages should be explained here.]

=end comment


=head1 RESTRICTIONS

=over

=item 1.

Each <name> must be an alphanumeric string (underscores are allowed).

=item 2.

The <name> tokens are replaced wherever they are preceded and followed
by non-word characters -- even inside C literal strings.

=item 3.

No <value> can contain a comma, except in a simple #edit directive.
In #*cases directives, no <value> can contain a parenthesis.

=back


=head1 BUGS

There is no support for Fortran or any other compiled language.

Error messages are relatively uninformative.  In particular, they fail
to mention the original line number of the offending directive.

The #line directives introduced by the B<--preserve> option can
break up a sequence of source lines that are meant to be joined by
continuation marks.  B<Cases> makes some attempt to avoid this problem:
it omits the #line directive that would otherwise be inserted
after a #*cases or #*loop directive if the last source line preceding
the #*cases or #*loop directive ends in a backslash.  This is not
a complete solution, however.


=head1 NOTES

Most of what B<cases> does to C code can be done in several other ways,
though less conveniently.  One can replace the body of a #*cases ...
#endcases or #*loop ... #endloop block by a macro and then invoke the
macro several times.  Or one can put the body into a separate file and
#include it several times, preceded and followed by appropriate #define
and #undef statements.  Or in C++, one can use templates.

These other solutions usually require some refactoring of the code,
which is inconvenient for experimentation and makes the code less
maintainable.  In contrast, with B<cases> one can often drop #ifcases
and #endcases directives around a block of code and see immediate
performance gains.


=head1 SEE ALSO

macdo(1CCR)


=head1 CLASSIFICATION

The source code is B<UNCLASSIFIED>.


=head1 COPYRIGHT

Copyright (c) 2020, Institute for Defense Analyses
4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500

All rights reserved.

This file is part of the conveyor package. For license information,
see the LICENSE file in the top level directory of the distribution.

=for html <div style="text-align: center; font-size: 120%; color: #006600; margin-top: 1ex;">
<strong>UNCLASSIFIED//FOR OFFICIAL USE ONLY</strong></div>

=cut


use integer;
use strict;
use Getopt::Long;

my $version = '2.2.2';
my $versiondate = '2018-12-03';

# Define a regular expression for matching directives
my $directives = '(calc|hexcalc|endcalc|loop|hexloop|endloop|edit|endedit|' .
    'defcase|encase|cases|defcases|ifcases|swcases|endcases)';

# Define templates here.  A template comprises a function that wraps
# a block after substitution, and three strings used to join the
# blocks together: an initial string, an interstitial string, and
# a final string.
my %templates = ( xif =>  [ \&ifCase, '', " else\n", "\n" ],
		  xdef => [ sub { $_[2] } , '', '', '' ],
                  xsw =>  [ \&swCase, '', " break;\n", " break;\n" ] );
$templates{'x'} = $templates{'xif'};

# Main routine starts here.  Process command-line options and arguments.
my ($doline, $doeval, $dofile, $donumber, $dopreserve, $comment) = (0, 0, 0, 0, 0, '//');
my ($doheader, $doversion) = (0, 0);
my @encasedirs = ('.');
my @extradefs = ();
if (! &GetOptions('arithmetic' => \$doeval, 'comment=s' => \$comment,
		  'defcase=s@' => \@extradefs, 'encase=s@' => \@encasedirs,
		  'file' => \$dofile, 'header' => \$doheader,
		  'line' => \$doline, 'number' => \$donumber,
		  'preserve' => \$dopreserve, 'version' => \$doversion)) {
    die "usage: perl cases.pl [options] [infile] > outfile\n" .
	"options: --comment=<string>, --encase=<directory>, --evaluate, --file,\n" .
        "          --header, --line, --number, --preserve, --version\n";
}

if ($doversion) {
    print "cases.pl version $version of $versiondate\n";
    exit;
}
die "cases: one file a time, please\n" if (1 < @ARGV);
die "cases: --file requires a filename\n" if (0 == @ARGV && $dofile);
die "cases: --preserve requires a filename\n" if (0 == @ARGV && $dopreserve);

# Read the input; handle line number and file name substitutions,
# and recursively process #encase directives
my $filename = $ARGV[0];
my @lines = <>;
&processFile(\@lines, $filename, 0);

# Use continuation marks to splice directives together
do {
    my $n = scalar @lines;
    my @newlines;
    for (my $i = 0; $i < $n; ) {
        my $j = 1;
	if ($lines[$i] =~ m/^\#\s*$directives\b/io) {
	    for (; $lines[$i] =~ m/\\\s*$/o && $i + $j < $n; $j++) {
	        substr($lines[$i],rindex($lines[$i],"\\")) = $lines[$i+$j];
	    }
        }
	push @newlines, $lines[$i];
	$i += $j;
    }
    @lines = @newlines;
};

# Convert --defcase options into #defcase directives
foreach my $extradef (reverse @extradefs) {
    my $ok = ($extradef =~ m/^(\w+)\s*=\s*(.*)$/o);
    unless ($ok) {
	print STDERR "malformed -d option: $extradef\n";
	next;
    }
    my $name = $1;
    my $value = $2;
    my $i = 0;
    for (; $i <= $#lines; $i++) {
	if ($lines[$i] =~ m/^\#\s*defcase\s+(\w+)\s*=/io) {
	    last if ($1 eq $name);
	}
    }
    my $line = "#defcase $extradef\n";
    $lines[$i] = $line if ($i <= $#lines);
    unshift @lines, $line if ($i > $#lines);
}

# In header mode, just convert #defcase to #define and exit
if ($doheader) {
    my $whence = $filename || 'standard input';
    print "/* Produced by perl cases.pl --header from $whence */\n";
    for (my $i = 0; $i <= $#lines; $i++) {
	if ($lines[$i] =~ m/^\#\s*defcase\s+(\w+)\s*=\s*(.*?)\s*$/io) {
	    print "\#define $1 $2\n";
	}
    }
    exit;
}

# Otherwise, convert #defcase to #defcases and add matching #endcases
do {
    my $count = 0;
    for (my $i = 0; $i <= $#lines; $i++) {
	$count += ($lines[$i] =~ s/^\#\s*defcase\b/\#defcases/io);
    }
    for (; $count > 0; $count--) {
	push @lines, "\#endcases\n";
    }
};

# Break the text into alternating chunks of program text (even indices)
# and #cases/#endcases directives (odd indices).  The "-1" ensures
# that the final chunk is program text, possibly null.
my $text = join '',@lines;
my @chunks = split /^(\#\s*(?:end||if|def|sw)cases\b.*?\n)/im, $text, -1;

# Repeatedly find the first innermost #cases...#endcases block and
# expand it.
for (my $i = 1; $i < @chunks - 3; $i += 2) {
    if ($chunks[$i] !~ m/^\#\s*end/io && $chunks[$i+2] =~ m/^\#\s*end/io) {
        my ($kind,$args) = ($chunks[$i] =~ m/^\#\s*([a-z]*)cases\b(.*)\n/io);
        my $newchunk = $chunks[$i-1]
            . &doCases($templates{'x'.$kind},$args,$chunks[$i+1])
            . $chunks[$i+3];
        splice(@chunks, $i-1, 5, $newchunk);
        $i -= ($i > 1 ? 4 : 2);
    }
}
die "unmatched directive: $chunks[1]\n" if (1 < @chunks);
$text = $chunks[0];

# Now process loop and calc directives in a similar way
@chunks = split /^(\#\s*(?:end||hex)(?:loop|calc|edit)\b.*?\n)/im, $text, -1;

# Repeatedly find an *outermost* #loop...#endloop or #calc...#endcalc
# or #edit...#endedit block and expand it
while (1 < @chunks) {
    if (4 >= @chunks || $chunks[1] =~ m/^\#\s*end/io) {
	die "unmatched directive: $chunks[1]\n";
    }
    my ($i,$depth) = (3,1);
    for (; $i < @chunks - 1; $i += 2) {
	$depth += ($chunks[$i] =~ m/^\#\s*end/io) ? -1 : 1;
	last if ($depth == 0);
    }
    die "unmatched directive: $chunks[1]\n" if ($depth > 0);
    my ($kind,$stem,$args) = ($chunks[1] =~ m/^\#\s*([a-z]*)(loop|calc|edit)\b(.*)\n/io);
    my %funcs = ('calc' => \&doCalc, 'edit' => \&doEdit, 'loop' => \&doLoop);
    my $func = $funcs{ lc($stem) };
    $text = $chunks[0] . &$func($kind, $args, (join '',@chunks[2..($i-1)])) . $chunks[$i+1];
    my @newlist = split /^(\#\s*(?:end||hex)(?:loop|calc|edit)\b.*?\n)/im, $text, -1;
    splice(@chunks, 0, $i+2, @newlist);
}
$text = $chunks[0];

if ($doeval) {
    $text =~ s/\#\#\[\[(.*?)\]\]/sprintf('%lX',&evalCExpr($1))/eg;
    $text =~ s/\#\#\[(.*?)\]/&evalCExpr($1)/eg;
}

print $text;
exit;


sub processFile {
    my ($liner, $filename, $level) = @_;

# Make sure the last line ends with a newline
    $liner->[-1] .= "\n" unless ($liner->[-1] =~ m/\n$/s);

# Strip comments from directives, including their continuation lines
    if ($comment) {
	my $n = $#$liner + 1;
	my $continuation;
	for (my $i = 0; $i < $n; $i++) {
	    if ($liner->[$i] =~ m/^\#\s*$directives\b/io) {
		do {
		    my $where = index $liner->[$i], $comment;
		    substr($liner->[$i],$where) = "\n" if ($where > 0);
		    $continuation = ($liner->[$i] =~ m/\\\s*$/o);
		    $i++ if ($continuation);
		} while ($continuation && $i < $n);
	    }
	}
    }

# Replace __LINE__ and/or __FILE__ if requested
    for (my $lineno = 1; $doline && $lineno <= @$liner; $lineno++) {
	$liner->[$lineno-1] =~ s/\b__LINE__\b/$lineno/g;
    }
    for (my $lineno = 1; $dofile && $lineno <= @$liner; $lineno++) {
	$liner->[$lineno-1] =~ s/\b__FILE__\b/\"$filename\"/g;
    }

# Add line numbers if requested, except to lines whose last non-blank
# character is a continuation mark (backslash)
    if ($donumber) {
	my $n = $#$liner + 1;
	for (my $i = 0; $i < $n; $i++) {
	    if ($liner->[$i] =~ m/^\#\s*$directives\b/io) {
		$i++ while ($liner->[$i] =~ m/\\\s*$/o && $i < $n);
	    }
	    elsif ($liner->[$i] !~ m/\\\s*$/o) {
		chop $liner->[$i];
		$liner->[$i] .= sprintf("\t\t%s =%d=\n", $comment, $i+1);
	    }
	}
    }

# Add #line directives if requested, skipping past continuation lines.
# Omit #line if the previous non-directive line ended with a continuation mark.
    if ($dopreserve) {
	my $n = $#$liner + 1;
	my @newlines;
	my $continu = 0;
	for (my $i = 0; $i < $n; $i++) {
	    push @newlines, $liner->[$i];
	    if ($liner->[$i] =~ m/^\#\s*$directives\b/io) {
		for (; $i+1 < $n && $liner->[$i] =~ m/\\\s*$/o; $i++) {
		    push @newlines, $liner->[$i+1];
		}
		my $encase = ($liner->[$i] =~ m/^\#\s*encase\b/io);
		push @newlines, sprintf("\#line %d\n", $i+2) unless ($encase or $continu);
		push @newlines, sprintf("\#line %d \"%s\"\n", $i+2, $filename) if
		    ($encase and not $continu);
	    }
	    else {
		$continu = ($liner->[$i] =~ m/\\\s*$/o);
	    }
	}
	@$liner = @newlines;
	unshift @$liner, "\#line 1 \"$filename\"\n";
    }

# Handle #encase directives recursively
    for (my $i = $#$liner; $i >= $0; $i--) {
	if ($liner->[$i] =~ m/\#\s*encase\s+"\s*(.*?)\s*"/io) {
	    my $subname = $1;
	    die "\#encase of $subname exceeds recursion limit of 10\n" if ($level >= 9);
	    my @encased;
	    my $subpath = &loadFile(\@encased, $subname);
	    &processFile(\@encased, $subpath, $level + 1);
	    splice @$liner, $i, 1, @encased;
	}
    }
}

# Return the path at which the file was found
sub loadFile {
    my ($liner, $filename) = @_;
    my @encpaths = ($filename);
    @encpaths = map "$_/$filename", @encasedirs unless ($filename =~ m|^/|o);
    foreach my $path (@encpaths) {
	if (-r $path) {
	    open ENCFILE, "<$path";
	    @$liner = <ENCFILE>;
	    close ENCFILE;
	    return $path;
	}
    }
    die "unable to open $filename for #encasement\n";
}

sub doCalc {
    my ($kind,$calc,$body) = @_;
    $calc =~ s/^\s(.*?)\s*$/$1/;
    if ($calc =~ m/^(\w+)\s*=\s*(.*?)\s*$/o) {
	my ($id, $value) = ($1, $2);
	$value = &evalCExpr($value);
	$value = sprintf('%lX',$value) if ($kind eq 'hex');
	return &singleSubst($body,$id,$value);
    }
    elsif ($calc =~ m/\(\s*([\w\s,]+?)\s*\)\s*=\s*\((.*?)\)\s*$/o) {
	my ($lhs, $rhs) = ($1, $2);
	my @ids = split /\s*,\s*/, $lhs;
	my @values = map { &evalCExpr($_) } split /\s*,\s*/, $rhs;
	@values = map { sprintf('%lX',$_) } @values if ($kind eq 'hex');
	return &multiSubst($body,\@ids,\@values);
    }
    else {
	print STDERR "unrecognized statement: #calc $calc\n";
	return '';
    }
}

sub doEdit {
    my ($kind,$edit,$body) = @_;
    $edit =~ s/^\s(.*?)\s*$/$1/;
    if ($edit =~ m/^(\w+)\s*=\s*(.*?)\s*$/o) {
	my ($id, $value) = ($1, $2);
	$value = &evalPerlExpr($value);
	return &singleSubst($body,$id,$value);
    }
    elsif ($edit =~ m/\(\s*([\w\s,]+?)\s*\)\s*=\s*\((.*?)\)\s*$/o) {
        my ($lhs, $rhs) = ($1, $2);
	my @ids = split /\s*,\s*/, $lhs;
	my @values = map { &evalPerlExpr($_) } split /\s*,\s*/, $rhs;
	return &multiSubst($body,\@ids,\@values);
    }
    else {
        print STDERR "unrecognized statement: #edit $edit\n";
        return '';
    }
}

sub doLoop {
    my ($kind,$loop,$body) = @_;
    $loop =~ s/^\s(.*?)\s*$/$1/;
    if ($loop =~ m/^(\w+)\s*=\s*(.*?)$/o) {
	my $id = $1;
	my @limits = split /\s*,\s*/, $2;
	if (2 != @limits && 3 != @limits) {
	    print STDERR "bad limits in loop: $loop\n";
	    return '';
	}
	@limits = map { &evalCExpr($_) } @limits;
	$limits[2] = 1 if (2 == @limits || 0 == $limits[2]);
	my $count = int(($limits[1] - $limits[0]) / $limits[2]);
	$count = int(($limits[0] - $limits[1]) / -$limits[2]) if ($limits[2] < 0);
	my @range = 0 .. $count;
	@range = map { $limits[0] + ($limits[2] * $_) } @range;
	if ($kind eq 'hex') {
	    my $format = '%0' . length(sprintf "%lX", $range[$#range]) . 'lX';
	    @range = map { sprintf $format, $_ } @range;
	}
	$body = join '', map { &singleSubst($body,$id,$_) } @range;
    }
    else {
	print STDERR "unrecognized directive: \#${kind}loop $loop\n";
	return '';
    }
    return $body;
}

sub evalCExpr {
    my ($expr) = @_;
    if ($expr !~ m/^([-+~!*\/\%<>=\&^|?:()\s]|\b0x[0-9a-f]+\b|\b[0-9]+\b)+$/io) {
	die "ill-formed integer expression: $expr\n";
    }
    my $value = eval $expr;
    die "bad expression \"$expr\": $@\n" if ($@);
    return $value;
}

# perform a limited set of string-manipulation operations
sub evalPerlExpr {
    my ($expr) = @_;
    my $simplify = $expr;
    $simplify =~ s/'([^'\\]|\\.)*?'//go;
    if ($simplify !~ m/^([.(,)\s]|\b[0-9]+\b|\b(chr|lc|lcfirst|substr|uc|ucfirst|x)\b)+$/io) {
        die "ill-formed string expression: $expr\n";
    }
    my $value = eval $expr;
    die "bad expression \"$expr\": $@\n" if ($@);
    return $value;
}

sub singleSubst {
    my ($text, $id, $value) = @_;
    $text =~ s/(\#\#)?\b$id\b(\#\#)?/$value/g;
    return $text;
}

sub multiSubst {
    my ($text, $ids, $values) = @_;
    if (scalar @$values != scalar @$ids) {
	print STDERR "(@$values) cannot be assigned to (@$ids)\n";
	return '';
    }
    for (my $i = 0; $i < @$ids; $i++) {
	$text =~ s/(\#\#)?\b$$ids[$i]\b(\#\#)?/$$values[$i]/g;
    }
    return $text;
}


sub doCases {
    my ($template,$cases,$body) = @_;
    $cases =~ s/^\s*(.*?)\s*$/$1/;
    my @blox;
    if ($cases =~ m/^(\w+)\s*=\s*(.*?)\s*$/o) {
	my $id = $1;
	my @values = split /\s*,\s*/, $2;
	@blox = map { &singleCases($$template[0],$body,$id,$_) } @values;
    }
    elsif ($cases =~ m/\(\s*([\w\s,]+?)\s*\)\s*=\s*(.*?)\s*$/o) {
	my ($lhs, $rhs) = ($1, $2);
	$rhs =~ s/\s*\(\s*(.*?)\s*\)\s*$/$1/;
	my @ids = split /\s*,\s*/, $lhs;
	my @tuples = split /\s*\)\s*,\s*\(\s*/, $rhs;
	@blox = map { &multiCases($$template[0],$body,\@ids,$_) } @tuples;
    }
    else {
	print STDERR "unrecognized statement: #cases $cases\n";
	return '';
    }
    my $result = join $$template[2], @blox;
    return $$template[1] . $result . $$template[3];
}

sub singleCases {
    my ($func, $text, $id, $value) = @_;
    $text =~ s/(\#\#)?\b$id\b(\#\#)?/$value/g;
    return &$func([ $id ], [ $value ], $text);
}

sub multiCases {
    my ($func, $text, $ids, $tuple) = @_;
    my @values = split /\s*,\s*/, $tuple, -1;
    if (scalar @values != scalar @$ids) {
	print STDERR "(@values) cannot be assigned to (@$ids)\n";
	return '';
    }
    for (my $i = 0; $i < @$ids; $i++) {
	$text =~ s/(\#\#)?\b$$ids[$i]\b(\#\#)?/$values[$i]/g;
    }
    return &$func($ids, \@values, $text);
}

# The arguments of a template function are:
# (0) a reference to a list of the identifiers
# (1) a reference to the corresponding list of values
# (2) the text to be wrapped, substitutions having been done

sub ifCase {
    my ($ids, $vals, $text) = @_;
    my @conds = map { $$ids[$_] . ' == (' . $$vals[$_] . ')' }
                grep { $$ids[$_] ne $$vals[$_] } (0 .. $#$vals);
    return (0 == @conds ? '' : 'if (' . (join ' && ', @conds) . ') ')
	. "{\n$text}";
}

sub swCase {
    my ($ids, $vals, $text) = @_;
    die "\#swcases cannot perform multiple substitution\n" if ($#$vals > 0);
    return ($$ids[0] eq $$vals[0] ? 'default' : 'case ' . $$vals[0])
	. ":{\n$text}";
}
