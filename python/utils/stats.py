#!/usr/bin/python

"""Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
import datetime
import subprocess

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

def print_stats(hash, date_line):
    """Calulate and display stats about the checkout in the current directory."""
    rev1_epoch = 1113991642
    revision = hash
    date_line = date_line[:-6] # Ignore daylight saving offset '.+0100' for portability
    pattern = '%a %b %d %H:%M:%S %Y'
    rev_time = time.strptime(date_line, pattern)
    rev_date = time.strftime('%Y-%m-%d', rev_time)
    epoch = time.mktime(rev_time) - rev1_epoch
    epoch_weeks = epoch / (3600*7*24.0)
  
    source_dirs = glob.glob('*/src')
    test_dirs = glob.glob('*/test')+glob.glob('*/tests')
  
    test_files = get_files(test_dirs)
    source_files = get_files(source_dirs)

    test_stats = file_stats(test_files)
    source_stats = file_stats(source_files)

    print revision,'\t',epoch_weeks,'\t',source_stats[0],'\t',source_stats[1],\
        '\t',test_stats[0],'\t',test_stats[1],'\t',test_stats[2],'\t',test_stats[3],\
        '\t',source_stats[1]+test_stats[1],'\t',rev_date

def print_header():
    """Print a TSV header line corresponding to the output of print_stats."""
    print '#rev\ttime (weeks)\tsrc_files\tsrc_loc\ttest_files\ttests_loc\ttest_suites\ttests\ttotal_loc\tdate'

def run(startDate):
    """Do the processing."""

    dir='../temp_lines_of_code'

    print '### Start date =', start_date
    print '###'
    erase_old = True
    #erase_old = False
    if erase_old:
        print '### Starting a fresh checkout in', dir
        print '###'
        if os.path.isdir(dir):
            print '### Erasing previous', dir
            shutil.rmtree(dir)
        print '### Checking out...'
        os.system('git clone -q -b develop https://chaste.cs.ox.ac.uk/git/chaste.git '+dir+' > /dev/null')
        print '### Checked out'
        print '###'
    print '### Switch to', dir
    print '###'
    os.chdir(dir)

    print_header()
    sys.stdout.flush()
    
    step = 3
    date = startDate + datetime.timedelta(20)
    old_hash = 'No hash'
    while (date <= datetime.date.today()):
        p = subprocess.Popen('git log develop --abbrev-commit --max-count=1 --before=' + str(date), shell=True, stdout=subprocess.PIPE, stderr=None)
        (out, _) = p.communicate()
        out_lines = out.split('\n')
        date_line = out_lines[2]
        if (date_line[:5] != 'Date:'):
            date_line = out_lines[3]
        date_line = date_line[8:] # Subtract "Date:..."
        hash =  out.split()[1]
        hash = hash[:-1] # Drop space
        if hash != old_hash:
            os.system('git checkout --quiet '+hash)
            print_stats(hash, date_line)
            sys.stdout.flush()
        date = date +datetime.timedelta(step)
        old_hash = hash
    
        
    

if __name__ == '__main__':
    start_date = datetime.date(2016,07,01)
    start_date = datetime.date(2017,02,01) # February
    # if len(sys.argv) > 1:
    #     start_date = int(sys.argv[1])
    run(start_date)

""" Cut'n'paste for gnuplot:

*****************
NB: 
SEE notforrelease/docs/loc.gnu 
*****************

set xlabel 'weeks'                                                                                                                    
set term png                                                                                                                          
set out 'loc.png'                                                                                                                     
plot 'loc.txt' u 2:4  w l title 'lines of source', 'loc.txt' u 2:6 w l title 'lines of tests', 'loc.txt' u 2:9 w l title 'total'
exit      
"""
