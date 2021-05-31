#!/usr/bin/python

import sys

if len(sys.argv) == 1:
    print '''
\n%s\n
A grep-like function. 
Takes a file that's a list to find in another file, specifying the column index to search in the targetfile. 
Ex: python grep_column.py listfile.txt targetfile.txt 0. Zero-indexed, whitespace delimited. 
%s 
''' % ( '#' * 20, '#' * 20)
    sys.exit(1)

else:
    try:
        listfile, targetfile, ncol = sys.argv[1:]
        ncol = int(ncol)
    except ValueError:
        print '\n%s\nError. Add listfile, targetfile, ncol, ie python grep_column.py listfile.txt targetfile.txt 0. %s ''' % ('#' * 20, '#' * 20)
        sys.exit(1)
    except TypeError:
        print '\n%s\nError. Add listfile, targetfile, ncol, ie python grep_column.py listfile.txt targetfile.txt 0. %s ''' % ('#' * 20, '#' * 20)
        sys.exit(1)

searchdict = dict([[line.strip(),1] for line in file(listfile).readlines()])
for line in file(targetfile):
    try:
        if line.strip().split()[ncol] in searchdict:
            print line.strip()
    except IndexError:
        # problem with line, continue
        continue
    
