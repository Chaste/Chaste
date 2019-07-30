#!/usr/bin/env python

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

import logging
import os
import sys
import validator
import translate

try:
    end = sys.argv.index('--')
    args = sys.argv[1:end]
    options = sys.argv[end+1:]
except ValueError:
    args = sys.argv[1:]
    options = []

models = [arg for arg in args if arg[0] != '-']
options.extend([arg for arg in args if arg[0] == '-'])

if '-h' in options:
    print "do_trans.py applies PE and LT optimisations to CellML models."
    print "Usage: %s [options] model ..." % sys.argv[0]
    print "Options:"
    print "  -q   redirect PyCml output to files"
    print "  -H   generate a CellML file for comparison with the Haskell implementation of PE"
    print "  -e   generate error analysis code"
    print "  -j   generate jacobian for error analysis code"
    print "  -np  don't use PE when generating Matlab code"
    print "  other options (i.e. args starting with '-' or appearing after '--')"
    print "  passed through to PyCml"
    sys.exit(0)

quiet = gen_error = gen_jacobian = haskell = False
do_pe = ' -p'
if '-q' in options:
    quiet = True
    options.remove('-q')
if '-H' in options:
    haskell = True
    options.remove('-H')
if '-e' in options:
    gen_error = True
    options.remove('-e')
if '-j' in options:
    gen_jacobian = True
    options.remove('-j')
if '-np' in options:
    do_pe = ''
    options.remove('-np')

if not models:
    print "Need a file to transform!"
    sys.exit(1)

options = ' ' + ' '.join(options)
prefix = ''
if gen_error or gen_jacobian:
    transforms = {}
    if gen_error:
        transforms.update({'-t Matlab' + do_pe: '.m',
                           '-t Matlab -l' + do_pe: '_lt.m',
                           '-t Maple --compute-full-jacobian': '.mpl'})
    if gen_jacobian:
        transforms.update({'J': '_J.m'})
elif haskell:
    transforms = {'-t Haskell': '.hs',
                  '-p -t Haskell': '_pe.hs'}
    prefix = 'Cml_'
else:
    transforms = {'': '.hpp',
                  '-l': '_lut.hpp',
                  '-p': '_pe.hpp',
                  '-l -p': '_pe_lut.hpp'}

v = validator.CellMLValidator()
notifier = translate.NotifyHandler(level=logging.WARNING_TRANSLATE_ERROR)
logging.getLogger('validator').addHandler(notifier)

def transform(model_file):
    if not os.path.exists(model_file):
        print model_file, "does not exist!"
        return
    model_base = os.path.splitext(model_file)[0]
    if prefix:
        # We need to add a prefix to the file name
        d = os.path.dirname(model_base)
        b = os.path.basename(model_base)
        model_base = os.path.join(d, prefix+b)

    # Check model is valid for transformation
    notifier.reset()
    valid = v.validate(model_file, warn_on_units_errors=True,
                       show_errors=False, show_warnings=False)
    if not valid or notifier.messages:
        print model_file, "is not translatable!"
        return

    for opt, suffix in transforms.iteritems():
        if quiet:
            f = model_base + suffix[:suffix.rfind('.')]
            fs = [f + '-out.txt', f + '-err.txt']
            q = ' >' + fs[0] + ' 2>' + fs[1]
        else:
            fs = []
            q = ''
        if opt == 'J':
            opt = '-t Matlab -j ' + model_base + '.out'
        cmd = './translate.py --assume-valid --config-file=config.xml --Wu -a '\
              + opt + ' -o ' + model_base + suffix + ' ' + model_file \
              + options + q
        print cmd
        os.system(cmd)
        for f in fs:
            try:
                size = os.path.getsize(f)
                if size == 0:
                    os.remove(f)
            except OSError:
                pass
    return

for model in models:
    transform(model)

v.quit()

if gen_error:
    print "On a machine with maple, do the following to compute Jacobians:"
    print "for f in *.mpl; do maple -i $f > ${f:0:${#f}-4}.out; done"
    print "Then re-run this script with the -j option."
