#!/usr/bin/python

"""Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.
"""

"""
Compute statistics about the evolution of the Chaste code base.
Looks at number of source/test files, lines of code, and number
of test suites & test cases.
"""

import glob
import os
import time
import shutil
import sys

def get_files(dirs):
    """Return a list of all hpp/cpp files within the given folders."""
    file_pairs=[]
    for root_dir in dirs:
        for root, _, files in os.walk(root_dir):
            for filename in files:
                if (filename[-2:] == 'pp'):
                    file_pairs.append([root,filename])
    return file_pairs


def file_stats(file_pairs):
    """Compute total statistics for the given files.
    Each entry in file_pairs is a (path, leafname) pair.
    Returns (number of files, total lines of code, num suites, num tests).
    The lines of code includes comment lines.
    """
    loc = 0
    nfiles = 0
    nsuites = 0
    ntests = 0
    for path, filename in file_pairs:
        loc += int(os.popen('wc -l '+path+'/'+filename).read().split()[0])
        nfiles += 1
        if (filename[:4] == 'Test'):
            nsuites+=1
            ntests+= int(os.popen('egrep -c -i "void\s+Test" '+path+'/'+filename).read().split()[0])
    return (nfiles, loc, nsuites, ntests)

def print_stats():
    """Calulate and display stats about the checkout in the current directory."""
    rev1_epoch=1113991642
    svn_info = os.popen('svn info').read().split('\n')
    rev_line=svn_info[4]
    revision = rev_line.split()[1]
    date_line = svn_info[9].split()
    date_time=date_line[3]+' '+date_line[4]
    pattern = '%Y-%m-%d %H:%M:%S'
    epoch = int(time.mktime(time.strptime(date_time, pattern))) - rev1_epoch
    epoch_weeks = epoch / (3600*7*24.0)
  
    source_dirs=glob.glob('*/src')
    test_dirs=glob.glob('*/test')+glob.glob('*/tests')
  
    test_files=get_files(test_dirs)
    source_files=get_files(source_dirs)

    test_stats=file_stats(test_files)
    source_stats=file_stats(source_files)

    print revision,'\t',epoch_weeks,'\t',source_stats[0],'\t',source_stats[1],\
        '\t',test_stats[0],'\t',test_stats[1],'\t',test_stats[2],'\t',test_stats[3],\
        '\t',source_stats[1]+test_stats[1]

def print_header():
    """Print a TSV header line corresponding to the output of print_stats."""
    print '#rev\ttime (weeks)\tsrc_files\tsrc_loc\ttest_files\ttests_loc\ttest_suites\ttests\ttotal_loc'

def run():
    """Do the processing."""
    svn_revision = os.popen("svnversion").read().strip()
    if (svn_revision[-1] == 'M'):
        svn_revision = svn_revision[0:-1]
    last_revision = int(svn_revision)

    dir='../temp_lines_of_code'

    print '### Starting a fresh checkout in', dir
    print '###'
    if os.path.isdir(dir):
        print '### Erasing previous', dir
        shutil.rmtree(dir)
    print '###'
    os.system('svn co -r 1 https://chaste.cs.ox.ac.uk/svn/chaste/trunk '+dir+' > /dev/null')
    os.chdir(dir)

    print_header()
    sys.stdout.flush()

    step = 10
    start = 10
    for rev in range(start,last_revision,step):
        os.system('svn up --non-interactive -r '+str(rev)+' > /dev/null')
        print_stats()
        sys.stdout.flush()

if __name__ == '__main__':
    run()

# Cut'n'paste for gnuplot:
#
# set xlabel 'weeks'                                                                                                                    
# set term png                                                                                                                          
# set out 'loc.png'                                                                                                                     
# plot 'nohup.out' u 2:4  w l title 'lines of source', 'nohup.out' u 2:6 w l title 'lines of tests', 'nohup.out' u 2:9 w l title 'total'
# exit      