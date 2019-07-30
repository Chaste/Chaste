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

"""
This part of PyCml deals with converting CellML models into programming language code.
It is a thin executable wrapper around translators.py.
"""

import os
import sys

# Make sure PyCml is on sys.path
pycml_path = os.path.dirname(os.path.realpath(__file__))
sys.path[0:0] = [pycml_path]

import translators # The actual code

if __name__ == '__main__':
    if '--profile' in sys.argv:
        import time, cProfile
        profile_name = '/tmp/pycml-profile-%f-%d' % (time.time(), os.getpid())
        cProfile.run('translators.run()', profile_name)
    else:
        translators.run()

    # For use in testing
    def euler(doc, t, nsteps=1000, dt=0.01):
        global tvar, state_vars, exprs
        tvar = t.free_vars[0]
        state_vars = t.state_vars
        for var in state_vars:
            var.set_value(float(var.initial_value))
        tvar.set_value(0.0)
        exprs = [e for e in doc.model.get_assignments()
                 if isinstance(e, translators.mathml_apply)]
        for _ in range(nsteps):
            for expr in exprs:
                expr.evaluate()
            tvar.set_value(tvar.get_value() + dt)
            for var in state_vars:
                var.set_value(var.get_value() +
                              dt * var.get_value(ode=tvar))
        return

    def writefile(doc, outfn='test.cml'):
        # Write out CellML file
        st = translators.open_output_stream(outfn)
        doc.xml(indent=1, stream=st)
        st.close()
        return

    def show_usage(doc):
        for comp in doc.model.component:
            for var in comp.variable:
                print var.fullname(), var._cml_usage_count


    def fix_divide_by_zero(doc):
        """
        Several models have equations of a form that may give rise to
        a divide by zero error on simulation, especially when lookup
        tables are used.  The general form is:
        
        (a * (V - v0)) / (exp(b * (V - v0)) - 1)
        
        When V = v0 this is undefined, however the limit of the
        function as V approaches v0 from either side is well-defined,
        and each limit is the same.  We approximate the limit by
        linear interpolation between values of the expression for
        (V-v0) = +/- 1e-10.
        """
        divides = [d.xml_parent
                   for d in doc.xml_xpath(u'//m:apply/m:divide')]
        for divide in divides:
            pass
        return
