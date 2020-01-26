# maximal/waves/chowder/ README

This 'chowder' sub-directory is used to hold sequence files currently under
attention or slated for attention during development of the 'maximal' program.
Because the benchmark scripts (typically named 'surfboard-cleanup_set_*') use
sequences in the 'waves/tubespit' sub-directory, some chowder sequences will be
located in duplicate forms in both sub-directories with the same file name. In those
cases, the duplicate file in the chowder sub-directory indicates that the
sequence contained in the file is either no longer solved correctly or else is
not solved optimally even if solved correctly. If a sequece has never been
solved correctly it should be located in the 'chowder' sub-directory even if it
is referenced in a 
testing script. 

Often an interesting sequence problem is whittled down to a "knot", a snippet of 
a fuller sequence that isolates the interesting problem under investigation. In 
some of these cases, a duplicate sub-sequence to one maintained in the other 
'waves' sub-directories will be also stored in 'chowder' with some knot desigation
in the file name.

*Last updated*: 1/26/2010 AJE
