#!/usr/bin/env python

"""Copyright (C) University of Oxford, 2005-2012

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
