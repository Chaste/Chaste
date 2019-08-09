# We want 1/2==0.5
from __future__ import division
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
This part of PyCml deals with converting CellML models into
programming language code, primarily C++ compatible with Chaste, but
supporting a few other languages also (and easily extensible).
"""

import optparse
import os
import re
import time
import sys
from cStringIO import StringIO

# Common CellML processing stuff
import pycml
from pycml import *  # Put contents in the local namespace as well
import optimize
import processors
import validator

__version__ = "$Revision$"[11:-2]

def version_comment(note_time=True):
    """Generate a version comment, with optional time info."""
    if note_time:
        t = '\non ' + time.asctime()
    else:
        t = ''
    text = """Processed by pycml - CellML Tools in Python
    (translators: %s, pycml: %s, optimize: %s)%s""" % (
        __version__, pycml.__version__, optimize.__version__, t)
    return text


def debugexpr(e):
    "For debugging."
    v = None
    if isinstance(e, cellml_variable):
        v = e
    elif isinstance(e, mathml_apply):
        v = e.assigned_variable()
    if v:
        r = (v==e, v.name, v.get_usage_count())
    else:
        r = (False, '', -1)
    return r


class TranslationError(RuntimeError):
    """Error thrown if CellML translation fails."""
    pass

class ConfigurationError(ValueError):
    """Error thrown if configuration file is invalid."""
    pass



class CellMLTranslator(object):
    """
    Base class for translators from CellML to programming languages.

    Provides various methods & attributes that can be overridden to
    achieve the desired output language and style.

    Also contains a registration system for subclasses, so the
    command-line client can know what translators are available.  See
    the register method for more information.
    """
    
    translators = {}
    class NameAlreadyRegistered(ValueError):
        pass
    @classmethod
    def register(cls, subclass, name):
        """Register a new translator subclass.

        Registers the subclass `subclass' with name `name' in the
        translators class attribute of CellMLTranslator.  If the name
        given is already in use, raises NameAlreadyRegistered.
        """
        if name in cls.translators:
            raise cls.NameAlreadyRegistered(name)
        cls.translators[name] = subclass
        return
    
    @staticmethod
    def generate_interface(doc, solver_info):
        """Generate an interface component connecting the model to whatever will use it.
        
        Stub method that subclasses can override to implement this functionality.
        """
        pass

    ###########################
    # Various language tokens #
    ###########################
    STMT_END = ';'                 # End of statement
    EQ_ASSIGN = ' = '              # Assignment operator
    COMMENT_START = '// '          # Start of a 1 line comment
    DOXYGEN_COMMENT_START = '//! ' # Start of a 1 line Doxygen comment
    LOGICAL_AND = ' && '
    LOGICAL_OR = ' || '
    LOGICAL_TRUE = ' true '
    ASSERT = 'EXCEPT_IF_NOT'       # Start of an assertion statement

    # Variable types
    TYPE_DOUBLE = 'double '
    TYPE_VOID = 'void '
    TYPE_CONST_DOUBLE = 'const double '
    TYPE_CONST_UNSIGNED = 'const unsigned '

    # Special constants
    TRUE = 'true'
    FALSE = 'false'
    PI = 'M_PI'
    E = 'M_E'
    NOT_A_NUMBER = 'NAN' # GNU extension, but fairly common
    
    # Whether the target language uses a subsidiary file, such as
    # a header file in C/C++
    USES_SUBSIDIARY_FILE = False
    # Mapping from primary file extension to subsidiary file extension
    FILE_EXTENSIONS = {'cpp': 'hpp',
                       'c': 'h',
                       'cxx': 'hxx'}

    def __init__(self, add_timestamp=True, options=None):
        """Create a translator."""
        self.options = options
        # Initially output should not be indented
        self.indent_level = 0
        # Character to indent with
        self.indent_char = ' '
        # No. of occurrences of indent_char per indent_level
        self.indent_factor = 4
        # Whether to use lookup tables where possible
        self.use_lookup_tables = True
        # Whether to add a timestamp comment to generated files
        self.add_timestamp = add_timestamp
        # Main output goes to the main file by default
        self._main_output_to_subsidiary = False

    def error(self, lines, xml=None):
        """Raise a translation error.

        lines is a list of strings describing what went wrong.
        A TranslationError with that message will be raised.

        If xml is given, it should be an element, which will be
        pretty-printed and included in the error.
        """
        if xml is not None:
            lines.extend(xml.xml(indent = u'yes',
                                 omitXmlDeclaration = u'yes').split('\n'))
        raise TranslationError('\n'.join(lines))

    @property
    def config(self):
        """Get the current document's configuration store."""
        return getattr(self.doc, '_cml_config', None)

    def translate(self, doc, model_filename, output_filename=None,
                  subsidiary_file_name=None,
                  class_name=None, v_variable=None,
                  continuation=None,
                  lookup_method_prefix='', row_lookup_method=False,
                  lt_index_uses_floor=True, constrain_table_indices=False):
        """Generate code for the given model.

        doc should be an instance of cellml_model representing a
        valid CellML model, such as might be produced from a call
        to
        >>> valid, doc = validator.CellMLValidator().validate(
        ...     model_filename, True)

        model_filename is the filename of the input model.
        The output program will by default be generated in the same
        folder, but with a different extension.  This can be
        overridden by supplying the output_filename keyword argument.

        By default the name of the class representing the model will
        be derived from the model name.  This can be overridden by
        passing an alternative as the class_name argument.

        The variable representing the transmembrane potential should
        be passed in using the v_variable argument.

        By default this method will perform some setup and then call
            self.output_top_boilerplate()
            self.output_mathematics()
            self.output_bottom_boilerplate()
        To alter this, pass a callable as the continuation parameter;
        this will then be called instead.

        lookup_method_prefix and row_lookup_method can be used to
        customise some aspects of lookup table usage.  The former is
        used by the Chaste translator to place lookup tables within a
        singleton class, while the latter can improve cache
        performance by looking up all tables in a single call, and
        returning an array of the results.
        
        lt_index_uses_floor specifies whether to use the floor()
        function to calculate the index into the lookup tables, or
        just cast to unsigned.

        constrain_table_indices specifies whether to throw an
        exception if lookup table index variables go outside the
        bounds specified (default), or just to cap them at the bounds.
        """
        self.doc = doc
        self.model = doc.model
        # Name of the class that will represent this model
        if class_name is None:
            self.class_name = u'CML_' + doc.model.name.replace('-', '_')
        else:
            self.class_name = class_name
        # Figure out the free & state vars in this model
        self.free_vars = doc.model.find_free_vars()
        self.state_vars = doc.model.find_state_vars()
        if len(self.free_vars) > 1:
            self.error(["Model has more than one free variable; exiting.",
                        "Free vars:" + str(self.free_vars)])
        if len(self.free_vars) == 0:
            if self.model.get_option('protocol'):
                # We may be dealing with an algebraic model; check for an Unknown variable
                for var in self.model.get_all_variables():
                    if var.get_type() == VarTypes.Unknown:
                        self.free_vars.append(var)
            if len(self.free_vars) != 1:
                self.error(["Model has no free variable; exiting."])
        # If only a single component, don't add it to variable names
        self.single_component = (len(getattr(self.model, u'component', [])) == 1)
        # Find the (index of the) transmembrane potential
        self.v_variable = v_variable
        if self.v_variable:
            self.v_variable_name = (v_variable.component.name, v_variable.name)
        else:
            self.v_variable = None
        for i, var in enumerate(self.state_vars):
            if var is v_variable:
                self.v_index = i
                break
        else:
            self.v_index = -1
        # Check to see if we're using lookup tables
        self.lookup_method_prefix = lookup_method_prefix
        self.row_lookup_method = row_lookup_method
        self.lt_index_uses_floor = lt_index_uses_floor
        self.constrain_table_indices = constrain_table_indices
        self.scan_for_lookup_tables()
        if not doc.lookup_tables:
            # No tables found
            self.use_lookup_tables = False

        # Extra configuration hook
        self.final_configuration_hook()

        # Open the output file(s)
        if output_filename is None:
            output_filename = self.output_file_name(model_filename)
        if self.USES_SUBSIDIARY_FILE:
            output_filename, self.subsidiary_filename = self.subsidiary_file_name(output_filename)
        self.output_filename = output_filename
        self.out = open_output_stream(output_filename)
        if self.USES_SUBSIDIARY_FILE:
            self.out2 = open_output_stream(self.subsidiary_filename)

        # Translate
        if continuation:
            continuation()
        else:
            self.output_top_boilerplate()
            self.output_mathematics()
            self.output_bottom_boilerplate()
        close_output_stream(self.out)
        if self.USES_SUBSIDIARY_FILE:
            close_output_stream(self.out2)
        return

    def final_configuration_hook(self):
        """A hook for subclasses to do some final configuration."""
        return

    def output_file_name(self, model_filename):
        """Generate a name for our output file, based on the input file."""
        return os.path.splitext(model_filename)[0] + '.cpp'
    
    def subsidiary_file_name(self, output_filename):
        """Generate a name for the subsidiary output file, based on the main one.
        
        Returns a pair (main_output_file_name, subsidiary_file_name).  This is in
        case the user specifies (e.g.) a .hpp file as the main output - we consider
        the main output to be the .cpp file.
        """
        base, ext = os.path.splitext(output_filename)
        ext = ext[1:] # Remove the '.'
        try:
            new_ext = self.FILE_EXTENSIONS[ext]
            swap = False
        except KeyError:
            swap = True
            for key, val in self.FILE_EXTENSIONS.iteritems():
                if val == ext:
                    new_ext = key
                    break
            else:
                # Turn off usage of subsidiary file
                self.USES_SUBSIDIARY_FILE = False
                return output_filename, None
        subsidiary_filename = base + '.' + new_ext
        if swap:
            output_filename, subsidiary_filename = subsidiary_filename, output_filename
        return output_filename, subsidiary_filename

    def send_main_output_to_subsidiary(self, to_subsidiary=True):
        """Set subsequent main-file writes to go to the subsidiary file instead.
        
        Supplying a False argument reverts this setting.
        
        Has no effect if not using a subsidiary file.
        """
        self._main_output_to_subsidiary = to_subsidiary

    def writeln(self, *args, **kwargs):
        """Write a line to our output file.

        Takes any number of strings as input, and writes them out to file.

        Unless the keyword argument indent is given as False, then the
        output will be indented to the level set by self.set_indent().
        Setting indent_level will override this value.
        Setting indent_offset will adjust the current level temporarily.

        If nl is set to False then a newline character will not be
        appended to the output.
        
        If subsidiary=True, then the line will be written to the subsidiary
        output file instead of the main one.  An error will be raised if
        there is no subsidiary output file.
        """
        if kwargs.get('subsidiary', False) or self._main_output_to_subsidiary:
            if not self.USES_SUBSIDIARY_FILE:
                self.error('Tried to write to non-existent subsidiary file')
            else:
                target = self.out2
        else:
            target = self.out
        indent = kwargs.get('indent', True)
        nl = kwargs.get('nl', True)
        if indent:
            level = kwargs.get('indent_level', self.indent_level)
            level += kwargs.get('indent_offset', 0)
            target.write(self.indent_char * self.indent_factor * level)
        target.write(''.join(map(str, args)))
        if nl:
            target.write('\n')

    def write(self, *args):
        """Write to our output file.

        This variant does not indent the output, or add a newline.
        """
        self.writeln(indent=False, nl=False, *args)
    
    def capture_output(self):
        """Make subsequent output operations write to a string buffer."""
        self._original_out = self.out
        self.out = StringIO()
    
    def get_captured_output(self):
        """Stop capturing output, and return what was captured as a string."""
        output = self.out.getvalue()
        self.out = self._original_out
        return output

    def output_comment(self, *args, **kwargs):
        """Output a (multi-line) string as a comment."""
        start = kwargs.get('start', self.COMMENT_START)
        if kwargs.get('pad', False):
            start = ' ' + start
        comment = ''.join(map(str, args))
        lines = comment.split('\n')
        for line in lines:
            self.writeln(start, line, **kwargs)

    def output_doxygen(self, *args, **kwargs):
        """Output a (multi-line) string as a Doxygen comment."""
        kwargs['start'] = self.DOXYGEN_COMMENT_START
        self.output_comment(*args, **kwargs)

    def set_indent(self, level=None, offset=None):
        """Set the indentation level for subsequent writes.

        If offset is given, adjust the level by that amount, otherwise
        set it to an absolute value.
        """
        if offset is not None:
            self.indent_level += offset
        else:
            self.indent_level = level

    def code_name(self, var, ode=False, prefix=None):
        """
        Return the full name of var in a form suitable for inclusion in a
        source file.
        
        If ode is True then return the name of the derivative of var
        instead.  We go directly to the source variable in this case,
        rather than including intermediate assignment statements as is
        done for connections.
        """
        if prefix is None:
            prefix = ['var_', 'd_dt_'][ode]
        if ode:
            var = var.get_source_variable(recurse=True)
        name = prefix + var.fullname(cellml=True)
        return name
    
    def varobj(self, varname):
        """Return the variable object that has code_name varname."""
        return cellml_variable.get_variable_object(self.model, varname)

    def var_display_name(self, var):
        """Return a display name for the given variable.
        
        If it has an oxmeta name, uses that.  Otherwise, looks first for a bqbiol:is annotation,
        or uses the cmeta:id if present, or the name attribute if not.  If there is an interface
        component, strip the name of it out of the display name.
        """
        if var.oxmeta_name:
            name = var.oxmeta_name
        else:
            for uri in var.get_rdf_annotations(('bqbiol:is', NSS['bqbiol'])):
                if '#' in uri:
                    name = uri[1 + uri.rfind('#'):]
                    break
            else:
                if hasattr(var, u'id') and var.id:
                    name = var.id
                else:
                    name = var.name
        iface = getattr(self.model, 'interface_component_name', '#N/A#')
        if name.startswith(iface):
            name = name[len(iface)+2:]
        return name

    @property
    def include_guard(self):
        """
        Get the include guard (for C/C++ output) for this cell model,
        based on the class name.
        """
        return self.class_name.upper() + '_HPP_'
    
    def output_top_boilerplate(self):
        """Output top boilerplate."""
        self.writeln('#ifndef _', self.include_guard, '_')
        self.writeln('#define _', self.include_guard, '_\n')
        self.output_comment('Model: ', self.model.name)
        self.output_comment(version_comment(self.add_timestamp))
        self.writeln()
        self.writeln('#include <cmath>')
        self.writeln('#include "AbstractOdeSystem.hpp"')
        self.writeln('#include "Exception.hpp"')
        self.writeln('#include "AbstractStimulusFunction.hpp"\n')
        self.writeln('class ', self.class_name, ' : public AbstractOdeSystem')
        self.writeln('{')
        self.writeln('private:')
        self.writeln('AbstractStimulusFunction *mpStimulus;\n',
                     indent_offset=1)
        self.writeln('public:')
        self.set_indent(1)
        self.writeln('const static unsigned _NB_OF_STATE_VARIABLES_ = ',
                   str(len(self.state_vars)), ';\n')
        self.writeln('//', ('-'*66))
        self.writeln('// Methods')
        self.writeln('//', ('-'*66), '\n')
        # Constructor
        self.writeln('', self.class_name,
                   '(AbstractStimulusFunction *stim)')
        self.writeln('    : AbstractOdeSystem(_NB_OF_STATE_VARIABLES_, ',
                     self.v_index, ')')
        self.open_block()
        self.writeln('mpStimulus = stim;\n')
        
        for i, var in enumerate(self.state_vars):
            self.writeln('// Y[', str(i), ']:')
            self.writeln('mVariableNames.push_back("', var.name, '");')
            self.writeln('mVariableUnits.push_back("', var.units, '");')
            init_val = getattr(var, u'initial_value', None)
            if init_val is None:
                init_comm = ' // Value not given in model'
                # Don't want compiler error, but shouldn't be a real number
                init_val = self.NOT_A_NUMBER
            else:
                init_comm = ''
            self.writeln('mInitialConditions.push_back(', init_val, ');',
                       init_comm, '\n')
        if self.use_lookup_tables: self.output_lut_generation()
        self.close_block()
        # Destructor
        self.writeln('~', self.class_name, '(void)')
        self.open_block()
        if self.use_lookup_tables: self.output_lut_deletion()
        self.close_block()
        # Lookup table declarations & methods
        if self.use_lookup_tables:
            self.output_lut_declarations()
            self.output_lut_methods()
        # Evaluation function
        self.writeln('void EvaluateYDerivatives (')
        self.writeln('        double ', self.code_name(self.free_vars[0]), ',')
        self.writeln('        const std::vector<double> &rY,')
        self.writeln('        std::vector<double> &rDY)')
        self.open_block()
        self.writeln('// Inputs:')
        self.writeln('// Time units: ', self.free_vars[0].units)
        for i, var in enumerate(self.state_vars):
            self.writeln('double ', self.code_name(var),
                         ' = rY[', str(i), '];')
            self.writeln('// Units: ', var.units, '; Initial value: ',
                         getattr(var, u'initial_value', 'Unknown'))
        self.writeln()
        if self.use_lookup_tables:
            self.output_table_index_generation()
        return

    def output_mathematics(self):
        """Output the mathematics in this model."""
        self.writeln(self.COMMENT_START, 'Mathematics')
        for expr in self.model.get_assignments():
            # Check this expression is actually used; don't output if not
            var = None
            if isinstance(expr, mathml_apply) and expr.is_assignment():
                var = expr.assigned_variable()
            elif isinstance(expr, cellml_variable):
                var = expr
            if not (var and var.get_usage_count() == 0):
                self.output_assignment(expr)
        return

    def output_bottom_boilerplate(self):
        """Output bottom boilerplate"""
        self.writeln('\n')
        for i, var in enumerate(self.state_vars):
            self.writeln('rDY[', str(i), '] = ', self.code_name(var, True),
                         ';')
        self.close_block()
        self.set_indent(offset=-1)
        self.writeln('};\n')
        self.writeln('#endif')
        return

    def output_assignment(self, expr):
        """Output an assignment expression."""
        if isinstance(expr, cellml_variable):
            # This may be the assignment of a mapped variable, or a constant
            t = expr.get_type()
            if t == VarTypes.Mapped:
                self.writeln(self.TYPE_CONST_DOUBLE, self.code_name(expr),
                             self.EQ_ASSIGN,
                             self.code_name(expr.get_source_variable()),
                             self.STMT_END, nl=False)
                self.output_comment(expr.units, indent=False, pad=True)
            elif t == VarTypes.Constant:
                self.writeln(self.TYPE_CONST_DOUBLE, self.code_name(expr),
                             self.EQ_ASSIGN, nl=False)
                self.output_number(expr.initial_value)
                self.writeln(self.STMT_END, indent=False, nl=False)
                self.output_comment(expr.units, indent=False, pad=True)
        else:
            # This is a mathematical expression
            self.writeln(self.TYPE_CONST_DOUBLE, nl=False)
            opers = expr.operands()
            self.output_lhs(opers.next())
            self.write(self.EQ_ASSIGN)
            self.output_expr(opers.next(), False)
            self.writeln(self.STMT_END, indent=False, nl=False)
            #1365: add a comment with the LHS units
            self.output_comment(expr._get_element_units(expr.eq.lhs, return_set=False).description(),
                                indent=False, pad=True)

    def output_lhs(self, expr):
        """Output the left hand side of an assignment expression."""
        if expr.localName == 'ci':
            self.output_variable(expr)
        elif expr.operator().localName == 'diff':
            self.write(self.code_name(expr.operator().dependent_variable, ode=True))

    def output_variable(self, ci_elt, ode=False):
        """Output a ci element, i.e. a variable lookup."""
        self.write(self.code_name(ci_elt.variable, ode=ode))

    def output_expr(self, expr, paren):
        """Output the expression expr.
        
        If paren is True then the context has requested parentheses around the
        output; if expr requires them then they will be added.
        """
        if self.use_lookup_tables and self.is_lookup_table(expr):
            self.output_table_lookup(expr, paren)
        elif isinstance(expr, mathml_apply):
            self.output_apply(expr, paren)
        elif isinstance(expr, mathml_piecewise):
            self.output_piecewise(expr, paren)
        elif isinstance(expr, mathml_ci):
            self.output_variable(expr)
        elif expr.localName == u'cn':
            self.output_number(expr)
        elif expr.localName == u'degree':
            # <degree> is just a wrapper around an expression
            self.output_expr(child_i(expr, 1), paren)
        elif expr.localName == u'logbase':
            # <logbase> is just a wrapper around an expression
            self.output_expr(child_i(expr, 1), paren)
        elif expr.localName == u'true':
            self.write(self.TRUE)
        elif expr.localName == u'false':
            self.write(self.FALSE)
        elif expr.localName == u'pi':
            self.write(self.PI)
        elif expr.localName == u'exponentiale':
            self.write(self.E)
        else:
            self.error(["Unsupported expression element " + expr.localName],
                       xml=expr)

    def output_number(self, expr):
        """Output the plain number expr.
        
        We make all constants parse as doubles to avoid problems with
        integer division or numbers too large for the int type.
        
        Negative numbers will be prefixed by a space to avoid unwanted
        decrement operations.
        """
        n = self.eval_number(expr)
        num = "%.17g" % n
        if num[0] == '-':
            num = ' ' + num
        if not '.' in num and not 'e' in num:
            num = num + '.0'
        self.write(num)

    def eval_number(self, expr):
        """Evaluate a number.

        If a (unicode) string, convert to float.
        If a cn element, call its evaluate method.
        """
        if isinstance(expr, mathml_cn):
            return expr.evaluate()
        else:
            return float(unicode(expr))

    # Map from operator element names to C++ function names,
    # where the translation is straightforward.
    function_map = {'power': 'pow', 'abs': 'fabs', 'ln': 'log', 'exp': 'exp',
                    'floor': 'floor', 'ceiling': 'ceil',
                    'factorial': 'factorial', # Needs external definition
                    'not': '!', 'rem': 'fmod',
                    'sin': 'sin', 'cos': 'cos', 'tan': 'tan',
                    'sec': '1/cos', 'csc': '1/sin', 'cot': '1/tan',
                    'sinh': 'sinh', 'cosh': 'cosh', 'tanh': 'tanh',
                    'sech': '1/cosh', 'csch': '1/sinh', 'coth': '1/tanh',
                    'arcsin': 'asin', 'arccos': 'acos', 'arctan': 'atan',
                    'arcsinh': 'asinh', 'arccosh': 'acosh', 'arctanh': 'atanh'}
    # Inverse reciprocal trig functions; these are represented by
    # key(x) = function_map[val](1/x)
    recip_trig = {'arcsec': 'arccos', 'arccsc': 'arcsin', 'arccot': 'arctan',
                  'arcsech': 'arccosh', 'arccsch': 'arcsinh', 'arccoth': 'arctanh'}
    # Operators
    nary_ops   = {'plus': '+', 'times': '*',
                  'and': '&&', 'or': '||'}
    binary_ops = {'divide': '/',
                  'xor': '!=', 'eq': '==', 'neq': '!=',
                  'geq': '>=','leq': '<=','gt': '>','lt': '<'}

    def output_apply(self, expr, paren):
        """Output an <apply> expression.
        
        paren is True if the context has requested parentheses.
        """
        op = expr.operator()
        if op.localName in self.function_map:
            self.output_function(self.function_map[op.localName],
                                 expr.operands(), paren)
        elif op.localName in self.recip_trig:
            self.output_function(self.function_map[self.recip_trig[op.localName]],
                                 expr.operands(), paren, reciprocal=True)
        elif op.localName == u'root':
            self.output_root(expr, paren)
        elif op.localName == u'log':
            self.output_log(expr, paren)
        elif op.localName in self.nary_ops:
            self.output_nary_operator(self.nary_ops[op.localName],
                                      expr.operands(), paren)
        elif op.localName in self.binary_ops:
            self.output_binary_operator(self.binary_ops[op.localName],
                                        expr.operands(), paren, expr)
        elif op.localName == u'minus':
            self.output_minus(expr, paren)
        elif op.localName == u'diff':
            # ODE occuring on the RHS
            self.write(self.code_name(op.dependent_variable, ode=True))
        elif op.localName == u'csymbol':
            self.output_special(expr, paren)
        else:
            # Unrecognised operator
            self.error(["Unsupported operator element " + str(op.localName)], xml=expr)

    def output_special(self, expr, paren):
        """Output a special-case operation, represented by a csymbol element.
        
        This needs to be implemented by subclasses if supported.
        """
        self.error(["Special-case operators are not supported by this translator."], xml=expr)

    def output_function(self, func_name, args, paren, reciprocal=False):
        """Output a function call with name func_name and arguments args.
        
        Parentheses are not required so paren is ignored.
        If reciprocal is True then pass the reciprocal of each arg to
        func_name.
        """
        self.write(func_name + '(')
        comma = False
        for arg in args:
            if comma: self.write(', ')
            else: comma = True
            if reciprocal:
                self.write('1/')
                self.output_expr(arg, True)
            else:
                self.output_expr(arg, False)
        self.write(')')

    def output_nary_operator(self, operator, operands, paren):
        """Output an n-ary operator (using infix notation).
        
        If paren is True, enclose the output in parentheses.
        """
        # TODO: Optimise - to use expm1(x) for computing exp(x)-1
        self.open_paren(paren)
        op = False
        for operand in operands:
            if op: self.write(' ' + operator + ' ')
            else: op = True
            self.output_expr(operand, True)
        self.close_paren(paren)

    def output_unary_operator(self, operator, operand, paren):
        """Output a unary operator (using prefix notation)."""
        self.open_paren(paren)
        self.write(operator)
        self.output_expr(operand, True)
        self.close_paren(paren)

    def output_binary_operator(self, operator, operands, paren, expr):
        """Output a binary operator.
        
        As output_nary_operator, but checks that len(list(operands)) == 2.
        """
        operands = list(operands)
        if len(operands) != 2:
            self.error(["Binary operator" + operator +
                        "does not have 2 operands."], xml=expr)
        self.output_nary_operator(operator, operands, paren)

    special_roots = {2: 'sqrt', 3: 'cbrt'}
    def output_root(self, expr, paren):
        """Output a root taken to some degree.

        If a degree qualifier element is not provided, uses default 2.
        """
        if hasattr(expr, u'degree'):
            # A degree is given.  Compute x^(1/b)
            # TODO: Optimise for when b==2 (sqrt) or b==3 (cbrt)
            # Try to evaluate expr.degree, and if the result is a key
            # of self.special_roots, use the value as the function to call.
            self.write('pow(')
            self.output_expr(expr.operands().next(), False)
            self.write(', 1/')
            self.output_expr(expr.degree, True)
            self.write(')')
        else:
            # Compute square root
            self.output_function('sqrt', expr.operands(), paren)

    def output_log(self, expr, paren):
        """Output a logarithm to the given base, which defaults to base 10."""
        if hasattr(expr, u'logbase'):
            # A base is provided.  Use the identity log_b(x) = log(x)/log(b)
            # TODO: Optimise for log2(x)
            self.open_paren(paren)
            self.output_function('log', expr.operands(), paren)
            self.write('/log(')
            self.output_expr(expr.logbase, False)
            self.write(')')
            self.close_paren(paren)
        else:
            # Use base 10
            self.output_function('log10', expr.operands(), paren)

    def output_minus(self, expr, paren):
        """Output either a unary or binary minus.

        Which is chosen depends on the number of operands.
        """
        operands = list(expr.operands())
        if len(operands) == 1:
            self.output_unary_operator('-', operands[0], paren)
        else:
            self.output_binary_operator('-', operands, paren, expr)

    def output_piecewise(self, expr, paren):
        """Output the piecewise expression expr.

        We use a cascading ternary if expression for simplicity.
        """
        self.open_paren(paren)
        for piece in getattr(expr, u'piece', []):
            self.output_expr(child_i(piece, 2), True) # Condition
            self.write(' ? ')
            self.output_expr(child_i(piece, 1), True) # Result
            self.write(' : ')
        if hasattr(expr, u'otherwise'):
            self.output_expr(child_i(expr.otherwise, 1), True) # Default case
        else:
            self.write(self.NOT_A_NUMBER)
        self.close_paren(paren)

    def open_paren(self, paren):
        if paren: self.write('(')
    def close_paren(self, paren):
        if paren: self.write(')')

    def open_block(self, **kwargs):
        """Open a new code block and increase indent."""
        self.writeln('{', **kwargs)
        self.set_indent(offset=1)
    def close_block(self, blank_line=True, **kwargs):
        """Close a code block and decrease indent."""
        self.set_indent(offset=-1)
        self.writeln('}', **kwargs)
        if blank_line:
            self.writeln(**kwargs)
        return

    ##############################
    # Dependency related methods #
    ##############################

    # These methods allow us to calculate which equations must be
    # output in order to compute a given set of quantities.
    def calculate_extended_dependencies(self, nodes, prune=[], prune_deps=[]):
        """Method moved to cellml_model."""
        return self.model.calculate_extended_dependencies(nodes, prune, prune_deps)

    def output_equations(self, nodeset):
        """Output the mathematics described by nodeset.

        nodeset represents a subset of the assignments in the model.
        Output assignments in the order given by a topological sort,
        but only include those in nodeset.

        Since a set of assignments is given, this method does not
        check whether variables are used - it is assumed that only
        assignments that are actually wanted are given in nodeset.
        """
        for expr in (e for e in self.model.get_assignments() if e in nodeset):
            self.output_assignment(expr)
        return

    def _vars_in(self, expr):
        """Return a set of variable objects used in the given expression.

        Will include state variables.  If the expression includes a derivative, the defining equation
        for that derivative will be included in the set.  Also if an expression is being
        replaced by a lookup table, this will only include the table key variable.
        """
        res = set()
        if self.use_lookup_tables and isinstance(expr, mathml) and self.is_lookup_table(expr):
            key_var = self.varobj(expr.getAttributeNS(NSS['lut'], u'var'))
            key_var = key_var.get_source_variable(recurse=True)
            res.add(key_var)
        elif isinstance(expr, mathml_ci):
            varobj = getattr(expr, '_cml_variable', None)
            if not varobj:
                varname = unicode(expr)
                varobj = self.varobj(varname.strip())
            if varobj:
                res.add(varobj)
        elif isinstance(expr, mathml_apply) and expr.operator().localName == u'diff':
            dep_varname = unicode(expr.ci)
            varobj = self.varobj(dep_varname.strip())
            res.add(varobj.get_ode_dependency(self.free_vars[0]))
        elif hasattr(expr, 'xml_children'):
            for child in expr.xml_children:
                res.update(self._vars_in(child))
        return res


    ########################
    # Lookup table methods #
    ########################

    # Lookup tables should be done in a cache- and memory-
    # efficient manner.
    #
    # Cache: Have one block of memory allocated for all tables with a
    # given index variable, such that entries are found at i*N+j,
    # where N is the no. of tables in the block, i is the index into a
    # table, and j is the table to read.  Change how lookups are done,
    # such that the lookup method is called once and returns a pointer
    # to the (i*N)'th entry.  Places where we now call the method then
    # index this pointer by j.
    #    The 'one block' part is done by default.
    #    The better lookup method is activated by --row-lookup-method.
    #
    # Memory: Extract the lookup tables into a separate class (in the
    # same .cpp file though).  This can then be made a singleton class
    # in a multi-cellular context.
    #    Chaste code generation has the option to do this, enabled by
    #    default.  Use --no-separate-lut-class to disable.

    def scan_for_lookup_tables(self):
        """Search for lookup tables used in this document.

        Store a list of suitable expressions on the document root.
        Generate a dictionary mapping tables to their index variables.
        """
        doc = self.doc
        # Remove xml:base to work around Amara bug!
        for elt in [doc, doc.model]:
            if u'base' in getattr(elt, 'xml_attributes', {}):
                print 'Delete base from', repr(elt)
                del elt.xml_attributes[u'base']
        # Get list of suitable expressions
        doc.lookup_tables = doc.xml_xpath(u"//*[@lut:possible='yes']")
        doc.lookup_tables.sort(cmp=element_path_cmp)
        # Map table keys (min, max, step, var) to an index variable
        doc.lookup_table_indexes = {}
        # Count the no. of lookup tables with each index variable
        doc.lookup_tables_num_per_index = {}
        if not doc.lookup_tables:
            # If no suitable expressions, we're done
            return
        # Search for table index variables already assigned
        table_indexes = [int(getattr(expr, u'table_index', -1))
                         for expr in doc.lookup_tables]
        tidx = max(table_indexes) + 1
        # Search for table names already assigned
        table_numbers = {}
        for expr in doc.lookup_tables:
            if hasattr(expr, u'table_name'):
                idx = expr.table_index
                table_numbers[idx] = max(table_numbers.get(idx, 0), 1 + int(expr.table_name))
        # Now assign new names, and construct mapping from tables to index variables
        for expr in doc.lookup_tables:
            # Get a suitable table index variable
            comp = expr.get_component()
            var = comp.get_variable_by_name(expr.var)
            var = var.get_source_variable(recurse=True)
            key = (expr.min, expr.max, expr.step, var)
            if not key in doc.lookup_table_indexes:
                var._cml_modifiable = True # Table index variables shouldn't be const, in case we constrain to table bounds
                if hasattr(expr, u'table_index'):
                    doc.lookup_table_indexes[key] = expr.table_index
                else:
                    doc.lookup_table_indexes[key] = unicode(tidx)
                    tidx += 1
                    expr.xml_set_attribute((u'lut:table_index', NSS['lut']),
                                           doc.lookup_table_indexes[key])
            # Get a table name, unique over all tables with this index variable
            if not hasattr(expr, u'table_name'):
                tnum = table_numbers.get(doc.lookup_table_indexes[key], 0)
                expr.xml_set_attribute((u'lut:table_name', NSS['lut']), unicode(tnum))
                table_numbers[doc.lookup_table_indexes[key]] = tnum + 1
        # Re-number table indices so they are contiguous starting from 0.
        table_index_map = {}
        table_name_map = {}
        tidx = 0
        for key in sorted(doc.lookup_table_indexes.keys()):
            idx = unicode(tidx)
            table_index_map[doc.lookup_table_indexes[key]] = idx
            table_name_map[idx] = {}
            doc.lookup_table_indexes[key] = idx
            doc.lookup_tables_num_per_index[idx] = 0
            tidx += 1
        # Make sure each lookup table is only listed once in doc.lookup_tables,
        # so we don't get 2 tables for the same expression!
        # Also re-number table names so they are contiguous starting from 0 for each table index.
        candidates = doc.lookup_tables[:]
        doc.lookup_tables = []
        listed = set()
        for expr in candidates:
            tid = (expr.table_index, expr.table_name)
            if not tid in listed:
                listed.add(tid)
                doc.lookup_tables.append(expr)
                # Renumber
                expr.table_index = table_index_map[expr.table_index]
                table_name_map[expr.table_index][expr.table_name] = unicode(doc.lookup_tables_num_per_index[expr.table_index])
                expr.table_name = table_name_map[expr.table_index][expr.table_name]
                # Increment count for this index variable
                doc.lookup_tables_num_per_index[expr.table_index] += 1
            else:
                # Just renumber to match the new id for this expression
                expr.table_index = table_index_map[expr.table_index]
                expr.table_name = table_name_map[expr.table_index][expr.table_name]
        return

    def lut_access_code(self, table_index, table_name, i):
        """Get the code for accessing the i'th element of the given table.
        """
        return '_lookup_table_%s[%s][%s]' % (table_index, i, table_name)
    
    def lut_parameters(self, key):
        """Get the bounds and step size for a particular table.
        
        key should be a key into self.lookup_table_indices.
        Returns (min, max, step, step_inverse) suitable for putting in generated code.
        """
        return key[0:3] + [unicode(1 / float(key[2]))]
    
    def lut_size_calculation(self, min, max, step):
        """Return the equivalent of '1 + (unsigned)((max-min)/step+0.5)'."""
        return '1 + (unsigned)((%s-%s)/%s+0.5)' % (max, min, step)

    def output_lut_generation(self, only_index=None):
        """Output code to generate lookup tables.

        There should be a list of suitable expressions available as self.doc.lookup_tables,
        to save having to search the whole model.
        
        If only_index is given, only generate tables using the given table index key.
        """
        # Don't use table lookups to generate the tables!
        self.use_lookup_tables = False
        # Allocate memory for tables
        for key, idx in self.doc.lookup_table_indexes.iteritems():
            if only_index is None or only_index == idx:
                min, max, step, _ = self.lut_parameters(key)
                self.writeln(self.TYPE_CONST_UNSIGNED, '_table_size_', idx, self.EQ_ASSIGN,
                             self.lut_size_calculation(min, max, step), self.STMT_END)
                self.writeln('_lookup_table_', idx, self.EQ_ASSIGN, 'new double[_table_size_', idx,
                             '][', self.doc.lookup_tables_num_per_index[idx], ']', self.STMT_END)
        # Generate each table in a separate loop
        for expr in self.doc.lookup_tables:
            var = expr.component.get_variable_by_name(expr.var)
            key = (expr.min, expr.max, expr.step, var.get_source_variable(recurse=True))
            idx = self.doc.lookup_table_indexes[key]
            if only_index is not None and only_index != idx:
                continue
            min, max, step, _ = self.lut_parameters(key)
            j = expr.table_name
            self.writeln('for (unsigned i=0 ; i<_table_size_', idx, '; i++)')
            self.open_block()
            self.writeln(self.TYPE_CONST_DOUBLE, self.code_name(var), self.EQ_ASSIGN, min,
                         ' + i*', step, self.STMT_END)
            self.writeln(self.lut_access_code(idx, j, 'i'), self.EQ_ASSIGN, nl=False)
            self.output_expr(expr, False)
            self.writeln(self.STMT_END, indent=False)
            self.close_block()
        self.use_lookup_tables = True

    def output_lut_deletion(self, only_index=None):
        """Output code to delete memory allocated for lookup tables."""
        for idx in self.doc.lookup_table_indexes.itervalues():
            if only_index is None or only_index == idx:
                self.writeln('if (_lookup_table_', idx, ')')
                self.open_block()
                self.writeln('delete[] _lookup_table_', idx, self.STMT_END)
                self.writeln('_lookup_table_', idx, self.EQ_ASSIGN, 'NULL', self.STMT_END)
                self.close_block(blank_line=False)

    def output_lut_declarations(self):
        """Output declarations for the lookup tables."""
        self.output_comment('Lookup tables')
        # Allocate memory, per index variable for cache efficiency
        for idx in self.doc.lookup_table_indexes.itervalues():
            num_tables = unicode(self.doc.lookup_tables_num_per_index[idx])
            self.writeln(self.TYPE_DOUBLE, '(*_lookup_table_', idx, ')[', num_tables, ']', self.STMT_END)
        self.writeln()

    def output_lut_index_declarations(self, idx):
        """Output declarations the variables used to index this table."""
        self.writeln('unsigned _table_index_', idx, self.STMT_END)
        factor = self.lut_factor(idx, include_type=True)
        if factor:
            self.writeln(factor, self.STMT_END)
        if self.row_lookup_method:
            self.writeln('double* _lt_', idx, '_row', self.STMT_END)

    def output_lut_indices(self):
        """Output declarations for the lookup table indices."""
        self.output_comment('Lookup table indices')
        for idx in self.doc.lookup_table_indexes.itervalues():
            self.output_lut_index_declarations(idx)
        self.writeln()
    
    def output_lut_methods(self):
        """Output the methods which look up values from lookup tables."""
        if self.row_lookup_method:
            self.output_lut_row_lookup_methods()
            return
        self.output_comment('Methods to look up values from lookup tables')
        self.output_comment('using ', self.config.options.lookup_type)
        for expr in self.doc.lookup_tables:
            j = expr.table_name
            idx = expr.table_index
            self.writeln('inline double _lookup_', j, '(unsigned i',
                         self.lut_factor('', include_type=True, include_comma=True), ')')
            self.open_block()
            self.output_single_lookup(idx, j, 'return ')
            self.close_block()
        self.writeln()
        return

    def output_single_lookup(self, tidx, tname, result):
        """Write the lookup calculation for a single entry.
        
        Used by output_lut_row_lookup_methods and output_lut_methods.
        """
        self.writeln(self.TYPE_CONST_DOUBLE, 'y1', self.EQ_ASSIGN,
                     self.lut_access_code(tidx, tname, 'i'), self.STMT_END)
        if self.config.options.lookup_type == 'linear-interpolation':
            self.writeln(self.TYPE_CONST_DOUBLE, 'y2', self.EQ_ASSIGN,
                         self.lut_access_code(tidx, tname, 'i+1'), self.STMT_END)
            self.writeln(result, 'y1 + (y2-y1)*', self.lut_factor(''), self.STMT_END)
        else:
            self.writeln(result, 'y1', self.STMT_END)

    def output_lut_row_lookup_methods(self):
        """Write methods that return a whole row of a lookup table.

        Note: assumes that table names are numbered sequentially from 0.
        """
        self.output_comment('Row lookup methods')
        self.output_comment('using ', self.config.options.lookup_type)
        for key, idx in self.doc.lookup_table_indexes.iteritems():
            num_tables = unicode(self.doc.lookup_tables_num_per_index[idx])
            self.writeln('double* _lookup_', idx, '_row(unsigned i',
                         self.lut_factor('', include_type=True, include_comma=True), ')')
            self.open_block()
            self.writeln('for (unsigned j=0; j<', num_tables, '; j++)')
            self.open_block()
            self.output_single_lookup(idx, 'j', '_lookup_table_%s_row[j] = ' % idx)
            self.close_block(False)
            self.writeln('return _lookup_table_', idx, '_row;')
            self.close_block()
        self.writeln()
        return
    
    def output_lut_row_lookup_memory(self):
        """Output declarations for the memory used by the row lookup methods."""
        self.output_comment('Row lookup methods memory')
        for key, idx in self.doc.lookup_table_indexes.iteritems():
            min, max, step, var = key
            num_tables = unicode(self.doc.lookup_tables_num_per_index[idx])
            self.writeln('double _lookup_table_', idx, '_row[', num_tables, '];')
        self.writeln()
        return

    def is_lookup_table(self, expr):
        """Return True iff expr can be replaced by a lookup table.

        Uses annotations from a previous analysis."""
        return expr.getAttributeNS(NSS['lut'], u'possible', '') == u'yes'
    
    def contained_table_indices(self, node):
        """Return all lookup tables used directly in computing this node.

        If this is an expression node, checks all its children for table
        lookups, and returns the set of table indices used.
        """
        result = set()
        if isinstance(node, amara.bindery.element_base):
            if self.is_lookup_table(node):
                result.add(node.table_index)
            else:
                for child in node.xml_children:
                    result.update(self.contained_table_indices(child))
        return result
    
    def lut_factor(self, idx, include_comma=False, include_type=False):
        """Return code for any extra factor needed to do a table lookup.
        
        Will return the empty string unless linear interpolation is being used.
        """
        if self.config.options.lookup_type == 'linear-interpolation':
            factor = '_factor_' + str(idx)
            if include_type: factor = self.TYPE_DOUBLE + factor
            if include_comma: factor = ', ' + factor
        else:
            factor = ''
        return factor

    def output_table_lookup(self, expr, paren):
        """Output code to look up expr in the appropriate table."""
        i = expr.table_index
        if self.row_lookup_method:
            self.write('_lt_', i, '_row[', expr.table_name, ']')
        else:
            self.write(self.lookup_method_prefix, '_lookup_', expr.table_name,
                       '(_table_index_', i, self.lut_factor(i, include_comma=True), ')')
        return

    def output_table_index_generation(self, time_name, nodeset=set()):
        """Output code to calculate indexes into any lookup tables.
        
        If time_name is given and table bounds are being checked, the time value will be included in the
        error message.  Note that we need to pass it in, since in some contexts the free variable is not
        defined.
        
        If nodeset is given, then filter the table indices calculated so
        that only those needed to compute the nodes in nodeset are defined.
        
        A nodeset is required if any table indices are computed variables rather than state variables.
        In this case, we use the equations within nodeset to calculate the values of the indices, and
        return a set containing just those nodes used, so we can avoid recalculating them later.
        """
        tables_to_index = set()
        nodes_used = set()
        for node in nodeset:
            tables_to_index.update(self.contained_table_indices(node))
        if tables_to_index or not nodeset:
            self.output_comment('Lookup table indexing')
        for key, idx in self.doc.lookup_table_indexes.iteritems():
            if not nodeset or idx in tables_to_index:
                var = key[-1]
                if var.get_type() is VarTypes.Computed:
                    if not nodeset:
                        raise TranslationError('Unable to generate lookup table indexed on', var, 'as it is a computed variable')
                    var_nodes = self.calculate_extended_dependencies([var]) & nodeset
                    self.output_equations(var_nodes)
                    nodes_used.update(var_nodes)
                self.output_table_index_checking(key, idx)
                if self.config.options.check_lt_bounds:
                    self.writeln('// LCOV_EXCL_START', indent=False)
                    self.writeln('if (_oob_', idx, ')')
                    if time_name is None:
                        dump_state_args = 'rY'
                    else:
                        dump_state_args = 'rY, ' + time_name
                    self.writeln('EXCEPTION(DumpState("', self.var_display_name(key[-1]),
                                 ' outside lookup table range", ', dump_state_args,'));', indent_offset=1)
                    self.writeln('// LCOV_EXCL_STOP', indent=False)
                self.output_table_index_generation_code(key, idx)
        self.writeln()
        return nodes_used
        
    def output_table_index_checking(self, key, idx):
        """Check whether a table index is out of bounds."""
        if self.config.options.check_lt_bounds:
            var = key[-1]
            min, max, _, _ = self.lut_parameters(key)
            varname = self.code_name(var)
            self.writeln('bool _oob_', idx, self.EQ_ASSIGN, 'false', self.STMT_END)
            self.writeln('if (', varname, '>', max, ' || ', varname, '<', min, ')')
            self.open_block()
            self.writeln('// LCOV_EXCL_START', indent=False)
            if self.constrain_table_indices:
                self.writeln('if (', varname, '>', max, ') ', varname, self.EQ_ASSIGN, max, self.STMT_END)
                self.writeln('else ', varname, self.EQ_ASSIGN, min, self.STMT_END)
            else:
                self.writeln('_oob_', idx, self.EQ_ASSIGN, 'true', self.STMT_END)
            self.writeln('// LCOV_EXCL_STOP', indent=False)
            self.close_block(blank_line=False)
    
    def output_table_index_generation_code(self, key, idx):
        """Method called by output_table_index_generation to output the code for a single table."""
        index_type = 'const unsigned '
        factor_type = 'const double '
        row_type = 'const double* const '
        var = key[-1]
        min, max, _, step_inverse = self.lut_parameters(key)
        offset = '_offset_' + idx
        offset_over_step = offset + '_over_table_step'
        varname = self.code_name(var)
        self.writeln(self.TYPE_CONST_DOUBLE, offset, self.EQ_ASSIGN, varname, ' - ', min, self.STMT_END)
        self.writeln(self.TYPE_CONST_DOUBLE, offset_over_step, self.EQ_ASSIGN,
                     offset, ' * ', step_inverse, self.STMT_END)
        idx_var = '_table_index_' + str(idx)
        if self.config.options.lookup_type == 'nearest-neighbour':
            if self.lt_index_uses_floor:
                self.writeln(index_type, idx_var, ' = (unsigned) round(', offset_over_step, ');')
            else:
                self.writeln(index_type, idx_var, ' = (unsigned) (', offset_over_step, '+0.5);')
        else:
            if self.lt_index_uses_floor:
                self.writeln(index_type, idx_var, ' = (unsigned) floor(', offset_over_step, ');')
            else:
                self.writeln(index_type, idx_var, ' = (unsigned)(', offset_over_step, ');')
            factor = self.lut_factor(idx)
            if factor:
                self.writeln(factor_type, factor, ' = ', offset_over_step, ' - ', idx_var, self.STMT_END)
        if self.row_lookup_method:
            self.writeln(row_type, '_lt_', idx, '_row = ', self.lookup_method_prefix, '_lookup_', idx,
                         '_row(', idx_var, self.lut_factor(idx, include_comma=True), ');')

class CellMLToChasteTranslator(CellMLTranslator):
    """
    As CellMLTranslator, but targets more recent Chaste style.

    Includes the ability to output a cell that can solve itself using
    backward Euler, if the appropriate analyses have been done on the
    model.  (See the -J and -j options to translate.py.)
    """

    # We want separate .cpp/.hpp files
    USES_SUBSIDIARY_FILE = True

    # Type of (a reference to) the state variable vector
    TYPE_VECTOR = 'std::vector<double> '
    TYPE_VECTOR_REF = 'std::vector<double>& '
    
    def writeln_hpp(self, *args, **kwargs):
        """Convenience wrapper for writing to the header file."""
        kwargs['subsidiary'] = True
        self.writeln(*args, **kwargs)

    def translate(self, *args, **kwargs):
        """Generate code for the given model."""
        our_kwargs = {'use_chaste_stimulus': False,
                      'separate_lut_class': True,
                      'convert_interfaces': False,
                      'use_modifiers': False,
                      'use_data_clamp': False,
                      'dynamically_loadable': False,
                      'use_protocol': False
                      }
        for key, default in our_kwargs.iteritems():
            setattr(self, key, kwargs.get(key, default))
            if key in kwargs:
                del kwargs[key]
        # Some other default settings
        self.use_backward_euler = False
        self.include_serialization = False
        # Last method's access specification
        self._last_method_access = 'private'
        return super(CellMLToChasteTranslator, self).translate(*args, **kwargs)

    def final_configuration_hook(self):
        """Set the LT method prefix (requires self.class_name to be set)."""
        if self.separate_lut_class:
            self.lt_class_name = self.class_name + '_LookupTables'
            self.lookup_method_prefix = self.lt_class_name + '::Instance()->'
        return super(CellMLToChasteTranslator, self).final_configuration_hook()
        
    def output_includes(self, base_class=None):
        """Output the start of each output file.
        
        As well as the #include lines, it also outputs the include guard for
        the .hpp file, and doxygen comment.
        
        If base_class is not None (and self.use_backward_euler isn't set)
        then includes that class' header instead of AbstractCardiacCell.
        
        If self.dynamically_loadable is set, includes extra headers needed
        for that case.
        
        Reads self.include_serialization and self.use_backward_euler.
        Sets self.base_class_name and self.class_inheritance.
        """
        self.writeln_hpp('#ifndef ', self.include_guard)
        self.writeln_hpp('#define ', self.include_guard, '\n')
        for sub in [False, True]:
            self.output_doxygen('@file\n\n',
                                'This source file was generated from CellML.\n\n',
                                'Model: ', self.model.name, '\n\n',
                                version_comment(self.add_timestamp),
                                '\n\n<autogenerated>',
                                subsidiary=sub)
            self.writeln(subsidiary=sub)
        # .cpp should include .hpp
        self.writeln('#include "', os.path.basename(self.subsidiary_filename), '"')
        if self.include_serialization:
            self.writeln_hpp('#include "ChasteSerialization.hpp"')
            self.writeln_hpp('#include <boost/serialization/base_object.hpp>')
        self.writeln('#include <cmath>')
        self.writeln('#include <cassert>')
        self.writeln('#include <memory>')
        if self.use_backward_euler:
            self.writeln_hpp('#include "AbstractBackwardEulerCardiacCell.hpp"')
            self.writeln('#include "CardiacNewtonSolver.hpp"')
            self.base_class_name = 'AbstractBackwardEulerCardiacCell<' + \
                str(self.nonlinear_system_size) + '>'
        elif self.options.rush_larsen:
            self.base_class_name = 'AbstractRushLarsenCardiacCell'
            self.writeln_hpp('#include "' + self.base_class_name + '.hpp"')
            if not self.doc._cml_rush_larsen:
                self.writeln('#include "Warnings.hpp"')
        elif self.options.grl1:
            self.base_class_name = 'AbstractGeneralizedRushLarsenCardiacCell'
            self.writeln_hpp('#include "' + self.base_class_name + '.hpp"')
        elif self.options.grl2: #1992 TODO: merge with above case
            self.base_class_name = 'AbstractGeneralizedRushLarsenCardiacCell'
            self.writeln_hpp('#include "' + self.base_class_name + '.hpp"')
        elif base_class:
            self.base_class_name = base_class
            self.writeln_hpp('#include "' + self.base_class_name + '.hpp"')
        else:
            self.base_class_name = 'AbstractCardiacCell'
            self.writeln_hpp('#include "' + self.base_class_name + '.hpp"')
        if self.use_modifiers:
            self.writeln_hpp('#include "AbstractCardiacCellWithModifiers.hpp"')
            self.writeln_hpp('#include "AbstractModifier.hpp"')
            # Modify the base class name
            self.base_class_name = 'AbstractCardiacCellWithModifiers<' + self.base_class_name + ' >'
        self.class_inheritance = ' : public ' + self.base_class_name
        if self.dynamically_loadable:
            self.writeln_hpp('#include "AbstractDynamicallyLoadableEntity.hpp"')
            self.class_inheritance += ', public AbstractDynamicallyLoadableEntity'
        if self.use_protocol:
            self.writeln_hpp('#include "AbstractTemplatedSystemWithOutputs.hpp"')
            self.class_inheritance += ', public AbstractTemplatedSystemWithOutputs<' + self.TYPE_VECTOR + '>'
        self.writeln('#include "Exception.hpp"')
        self.writeln('#include "OdeSystemInformation.hpp"')
        self.writeln('#include "RegularStimulus.hpp"')
        self.writeln_hpp('#include "AbstractStimulusFunction.hpp"')
        self.writeln('#include "HeartConfig.hpp"')
        self.writeln('#include "IsNan.hpp"')
        self.writeln('#include "MathsCustomFunctions.hpp"')
        self.writeln()
        self.writeln_hpp()
        
    def set_access(self, access):
        """Set the access specification for subsequent output.
        
        We keep track of the last access set, either via this method or
        output_method_start, and only output a new declaration to the
        header file if it changes.
        """
        if access != self._last_method_access:
            self._last_method_access = access
            self.writeln_hpp()
            self.writeln_hpp(access, ':', indent_offset=-1)

    def output_method_start(self, method_name, args, ret_type, access=None, defaults=[]):
        """Output the start of a method declaration/definition.
        
        Will write to both the .hpp and .cpp file.
        
        We keep track of the access of the last method, and only output a new
        declaration to the header file if it changes.  The default is to use
        the same access specification as last time.
        """
        DEBUG('translator', 'Generating code for method', method_name)
        if access:
            self.set_access(access)
        if ret_type:
            if ret_type[-1] != ' ':
                ret_type = ret_type + ' '
        else:
            ret_type = ''
        args_string_cpp = ', '.join(filter(None, map(str, args)))
        if defaults:
            assert len(defaults) == len(args)
            args_with_default = []
            for (arg, default) in zip(map(str, args), map(str, defaults)):
                if arg:
                    if default:
                        args_with_default.append(arg + '=' + default)
                    else:
                        args_with_default.append(arg)
            args_string_hpp = ', '.join(args_with_default)
        else:
            args_string_hpp = args_string_cpp
        self.writeln_hpp(ret_type, method_name, '(', args_string_hpp, ')', self.STMT_END)
        self.writeln(ret_type, self.class_name, '::', method_name, '(', args_string_cpp, ')')

    def output_derived_quantities(self):
        """Output a ComputeDerivedQuantities method if any such quantities exist.
        
        Looks for variables annotated with pycml:derived-quantity=yes, and generates
        a method to compute all these variables from a given state.
        """
        dqs = self.derived_quantities
        if dqs:
            self.output_method_start('ComputeDerivedQuantities',
                                     [self.TYPE_DOUBLE + self.code_name(self.free_vars[0]),
                                      'const ' + self.TYPE_VECTOR + '& rY'], # We need it to really be a reference
                                     self.TYPE_VECTOR, access='public')
            self.open_block()
            self.output_comment('Inputs:')
            self.output_comment('Time units: ', self.free_vars[0].units)
            # Work out what equations are needed
            if self.use_chaste_stimulus:
                i_stim = [self.doc._cml_config.i_stim_var]
            else:
                i_stim = []
            if self.use_data_clamp:
                prune = [self.config.i_data_clamp_data]
            else:
                prune = []
            nodeset = self.calculate_extended_dependencies(dqs, prune_deps=i_stim, prune=prune)
            # State variable inputs
            self.output_state_assignments(assign_rY=False, nodeset=nodeset)
            self.writeln()
            table_index_nodes_used = self.calculate_lookup_table_indices(nodeset, self.code_name(self.free_vars[0]))
            table_index_nodes_used.update(self.output_data_table_lookups(nodeset - table_index_nodes_used))
            # Output equations
            self.output_comment('Mathematics')
            self.output_equations(nodeset - table_index_nodes_used)
            self.writeln()
            # Assign to results vector
            self.writeln(self.vector_create('dqs', len(dqs)))
            for i, var in enumerate(dqs):
                self.writeln(self.vector_index('dqs', i), self.EQ_ASSIGN, self.code_name(var), self.STMT_END)
            self.writeln('return dqs', self.STMT_END)
            self.close_block(blank_line=True)
    
    def output_serialize_method(self):
        """This method outputs the boost serialize method for the 
        header files that need it."""
        # Serialization
        if self.include_serialization:
            self.writeln_hpp('friend class boost::serialization::access;')
            self.writeln_hpp('template<class Archive>')
            self.writeln_hpp('void serialize(Archive & archive, const unsigned int version)')
            self.open_block(subsidiary=True)
            self.writeln_hpp('archive & boost::serialization::base_object<', self.base_class_name,
                             ' >(*this);')
            if self.dynamically_loadable:
                self.writeln_hpp('archive & boost::serialization::base_object<AbstractDynamicallyLoadableEntity>(*this);')
            if self.use_modifiers:
                self.output_comment('Despite this class having modifier member variables, they are all added to the', subsidiary=True)
                self.output_comment('abstract class by the constructor, and archived via that, instead of here.', subsidiary=True) 
            self.close_block(subsidiary=True)
       
    def output_cell_parameters(self):
        """Output declarations, set & get methods for cell parameters.
        
        Sets self.cell_parameters to be those constant variables annotated with
        pycml:modifiable-parameter.  These use the mParameters functionality in
        Chaste.
        
        Also collects any variables annotated with an RDF oxmeta name into
        self.metadata_vars. Only constants and state variables are included.
        """
        # Find annotated parameters
        self.cell_parameters = filter(
            lambda v: v.is_modifiable_parameter,
            cellml_metadata.find_variables(self.model,
                                           ('pycml:modifiable-parameter', NSS['pycml']),
                                           'yes'))

        # Reduce intra-run variation
        self.cell_parameters.sort(key=self.var_display_name)

        for i, var in enumerate(self.cell_parameters):
            # Remember the var's index
            var._cml_param_index = i

        # Create set of all oxmeta-annotated variables
        vars = cellml_metadata.find_variables(self.model, ('bqbiol:is', NSS[u'bqbiol']))
        # Keep only the variables with an oxmeta name
        vars = filter(lambda v: v.oxmeta_name, vars)
        # We're interested in anything that isn't time or the stimulus
        self.metadata_vars = set([v for v in vars if v.get_type() != VarTypes.Free])
        self.metadata_vars.discard(self.doc._cml_config.i_stim_var)
        self.metadata_vars = list(self.metadata_vars)
        self.metadata_vars.sort(key=self.var_display_name)
        
        # #1464 Create a set of metadata variables that will have modifiers
        # We want to avoid writing out metadata for stimulus current as it is used once and then discarded.
        # \todo - use protocol information to put only the required modifiers into this list.
        self.modifier_vars = [v for v in self.metadata_vars if v.oxmeta_name not in cellml_metadata.STIMULUS_NAMES and v.oxmeta_name != 'membrane_capacitance']
        self.modifier_vars.sort(key=self.var_display_name)

        # Generate member variable declarations
        self.set_access('private')
        if self.metadata_vars:
            self.output_comment('\nSettable parameters and readable variables\n', subsidiary=True)
        
        # Write out the modifier member variables. 
        if self.use_modifiers:
            for var in self.modifier_vars:
                self.writeln_hpp('boost::shared_ptr<AbstractModifier> mp_' + var.oxmeta_name + '_modifier', self.STMT_END)    
        
        # Methods associated with oxmeta annotated variables
        # Don't use LT & modifiers for the const methods
        use_modifiers = self.use_modifiers
        self.use_modifiers = False
        use_lt = self.use_lookup_tables
        self.use_lookup_tables = False
        for var in self.metadata_vars:
            if var.is_statically_const(ignore_annotations=True):
#                 self.output_method_start('Get_' + var.oxmeta_name + '_constant', [], self.TYPE_DOUBLE)
#                 self.open_block()
#                 self.output_comment('Constant value given in CellML')
#                 nodeset = self.calculate_extended_dependencies([var])
#                 self.output_equations(nodeset)
#                 self.writeln('return ', self.code_name(var), self.STMT_END)
#                 self.close_block()
#                 self.writeln()
                if var in self.cell_parameters and var in self.modifier_vars:
                    # 'Forget' its index, so normal code generation occurs (#1647)
                    var._cml_has_modifier = True
        self.use_lookup_tables = use_lt
        self.use_modifiers = use_modifiers
        self.output_default_stimulus()
        self.output_intracellular_calcium()
        
        # Find & store derived quantities, for use elsewhere
        self.derived_quantities = cellml_metadata.find_variables(self.model,
                                                                 ('pycml:derived-quantity', NSS['pycml']),
                                                                 'yes')
        # Reduce intra-run variation
        self.derived_quantities.sort(key=self.var_display_name)
                
    def output_default_stimulus(self):
        """
        Output a default cell stimulus from the metadata specification
        as long as the following metadata exists:
         * membrane_stimulus_current_amplitude
         * membrane_stimulus_current_duration
         * membrane_stimulus_current_period
        and optionally:
         * membrane_stimulus_current_offset
         * membrane_stimulus_current_end
        
        Ensures that the amplitude of the generated RegularStimulus is negative.
        """
        vars = dict()
        for n in ['duration', 'amplitude', 'period', 'offset', 'end']:
            vars[n] = self.model.get_variable_by_oxmeta_name('membrane_stimulus_current_'+n, throw=False)
        if not (vars['duration'] and vars['amplitude'] and vars['period']):
            self.has_default_stimulus = False
            return
        self.has_default_stimulus = True
        nodeset = self.calculate_extended_dependencies(filter(None, vars.values()))

        self.output_method_start('UseCellMLDefaultStimulus', [], 'boost::shared_ptr<RegularStimulus>', 'public')
        self.open_block()
        self.output_comment('Use the default stimulus specified by CellML metadata')
        self.output_equations(nodeset)
        self.writeln('boost::shared_ptr<RegularStimulus> p_cellml_stim(new RegularStimulus(')
        self.writeln('        -fabs(', self.code_name(vars['amplitude']), '),')
        self.writeln('        ', self.code_name(vars['duration']), ',')
        self.writeln('        ', self.code_name(vars['period']), ',')
        if vars['offset']:
            self.writeln('        ', self.code_name(vars['offset']))
        else:
            self.writeln('        0.0')
        if vars['end']:
            self.writeln('      , ', self.code_name(vars['end']))
        self.writeln('        ))', self.STMT_END)
        self.writeln('mpIntracellularStimulus = p_cellml_stim', self.STMT_END)
        self.writeln('return p_cellml_stim', self.STMT_END)
        self.close_block(blank_line=True)
    
    def output_intracellular_calcium(self):
        """
        If a (state) variable has been annotated as cytosolic_calcium_concentration,
        generate a GetIntracellularCalciumConcentration method.
        """
        # Find cytosolic_calcium_concentration
        cai = self.doc.model.get_variable_by_oxmeta_name('cytosolic_calcium_concentration', throw=False)
        if cai and cai in self.state_vars:
            i = self.state_vars.index(cai[0])
            self.output_method_start('GetIntracellularCalciumConcentration', [], self.TYPE_DOUBLE, 'public')
            self.open_block()
            self.writeln('return ', self.vector_index('mStateVariables', i), self.STMT_END)
            self.close_block(blank_line=True)
        
    def code_name(self, var, *args, **kwargs):
        """
        Return the full name of var in a form suitable for inclusion in a source file.
        
        Overrides the base class version to access mParameters for parameters.
        """
        if hasattr(var, '_cml_param_index') and not (self.use_modifiers and getattr(var, '_cml_has_modifier', False)):
            return self.vector_index('mParameters', var._cml_param_index)
        elif var is getattr(self.model, u'_cml_Chaste_Cm', None):
            return 'HeartConfig::Instance()->GetCapacitance()'
        elif hasattr(var, '_cml_code_name'):
            return var._cml_code_name % {'time': self.code_name(self.free_vars[0])}
        else:
            return super(CellMLToChasteTranslator, self).code_name(var, *args, **kwargs)

    def output_top_boilerplate(self):
        """Output top boilerplate.
        
        This method outputs the constructor and destructor of the cell
        class, and also lookup table declarations and lookup methods.
        It also calls output_verify_state_variables.
        """
        self.include_serialization = True

        # Check if we're generating a Backward Euler model
        self.use_backward_euler = self.model.get_option('backward_euler')
        self.use_analytic_jacobian = (self.model.get_option('maple_output') and hasattr(self.model.solver_info, u'jacobian'))
        if self.use_backward_euler:
            assert hasattr(self.model, u'solver_info')
            # Find the size of the nonlinear system
            num_linear_odes = len(self.model.solver_info.xml_xpath(u'solver:linear_odes/m:math/m:apply'))
            self.nonlinear_system_size = len(self.state_vars) - 1 - num_linear_odes
            nonlinear_entries = self.model.solver_info.xml_xpath(u'solver:jacobian/solver:entry/@var_j')
            self.nonlinear_system_vars = map(self.varobj, nonlinear_entries[:self.nonlinear_system_size])
        # Start output
        self.output_includes()
        
        if self.use_backward_euler or self.options.rush_larsen or self.options.grl1 or self.options.grl2:
            # Keep the same signature as forward cell models, but note that the solver isn't used
            solver1 = 'boost::shared_ptr<AbstractIvpOdeSolver> /* unused; should be empty */'
            solver2 = ''
            #solver1 = solver2 = ''
        else:
            solver1 = 'boost::shared_ptr<AbstractIvpOdeSolver> pSolver'
            solver2 = 'pSolver'

        if self.use_lookup_tables and self.separate_lut_class:
            self.output_lut_class()

        # Cell model class
        self.writeln_hpp('class ', self.class_name, self.class_inheritance)
        self.open_block(subsidiary=True)
        
        # Put the boost serialize() method in if requested.
        self.output_serialize_method()

        # Parameter declarations, and set & get methods (#666)
        self.output_cell_parameters()
        # Constructor
        self.set_access('public')
        self.output_constructor([solver1, 'boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus'],
                                [solver2, len(self.state_vars), self.unsigned_v_index, 'pIntracellularStimulus'])
        # Destructor
        self.output_method_start('~'+self.class_name, [], '')
        self.open_block()
        self.close_block()
        # Other declarations & methods
        self.output_chaste_lut_methods()
        self.output_verify_state_variables()
        return
    
    @property
    def unsigned_v_index(self):
        if self.v_index == -1:
            return 'UNSIGNED_UNSET'
        else:
            return str(self.v_index)
    
    def output_verify_state_variables(self):
        """Output the VerifyStateVariables method.
        
        This will look for state variables annotated with pycml:range-low and/or pycml:range-high,
        which specify allowable ranges for these variables.  The generated method will check that
        they are within the range.  Both limits are included, i.e. they specify a closed interval.
        """
        
        # First work out if there are any constraints on state variables
        low_prop = ('pycml:range-low', NSS['pycml'])
        high_prop = ('pycml:range-high', NSS['pycml'])
        low_range_vars = filter(
            lambda v: v.get_type() == VarTypes.State,
            cellml_metadata.find_variables(self.model, low_prop))
        high_range_vars = filter(
            lambda v: v.get_type() == VarTypes.State,
            cellml_metadata.find_variables(self.model, high_prop))
        nodeset = set(low_range_vars + high_range_vars)
        
        # If not, don't bother writing the method, an empty implementation is in the abstract classes.
        if nodeset:
            # It's not appropriate to apply modifiers here - we want to check the actual values of the state
            use_modifiers = self.use_modifiers
            self.use_modifiers = False
        
            self.output_method_start('VerifyStateVariables', [], 'void')
            self.open_block()
            
            using_cvode = (self.TYPE_VECTOR_REF == CellMLToCvodeTranslator.TYPE_VECTOR_REF)
            if using_cvode:
                self.writeln('/* We only expect CVODE to keep state variables to within its tolerances,')
                self.writeln(' * not exactly the bounds prescribed to each variable that are checked here.')
                self.writeln(' *')
                self.writeln(' * For 99.99% of paces this->mAbsTol works,')
                self.writeln(' * For 99.999% of paces 10*this->mAbsTol is fine,')
                self.writeln(' * but unfortunately 100x seems to be required on rare occasions for upstrokes.')
                self.writeln(' * This sounds bad, but is probably typically only 1e-5 or 1e-6.')
                self.writeln(' */')
                self.writeln('const double tol = 100*this->mAbsTol;')
                            
            self.output_state_assignments(nodeset=nodeset)
            error_template = 'EXCEPTION(DumpState("State variable {0} has gone out of range. Check numerical parameters, for example time and space stepsizes, and/or solver tolerances"));'
            additional_tolerance_adjustment = ''
            for var in low_range_vars:
                if using_cvode:
                    additional_tolerance_adjustment = ' - tol'
                self.writeln('if (', self.code_name(var), ' < ', var.get_rdf_annotation(low_prop), additional_tolerance_adjustment, ')')
                self.open_block()
                #self.writeln('std::cout << "Too small: ', self.code_name(var), ' = " << ', self.code_name(var) , ' << std::endl << std::flush;')
                self.writeln(error_template.format(self.var_display_name(var)))
                self.close_block(False)
            for var in high_range_vars:
                if using_cvode:
                    additional_tolerance_adjustment = ' + tol'            
                self.writeln('if (', self.code_name(var), ' > ', var.get_rdf_annotation(high_prop), additional_tolerance_adjustment, ')')
                self.open_block()
                #self.writeln('std::cout << "Too large: ', self.code_name(var), ' = " << ', self.code_name(var) , ' << std::endl << std::flush;')
                self.writeln(error_template.format(self.var_display_name(var)))
                self.close_block(False)
            self.close_block(True)
            
            self.use_modifiers = use_modifiers
  
    def output_constructor(self, params, base_class_params):
        """Output a cell constructor.
        
        params is a list of constructor parameters, entries of which should be strings
        including both type and parameter name, which will be included verbatim in the
        generated code.
        
        base_class_params is a list of parameters to be supplied to the base class
        constructor.  Entries will be converted to strings.
        """
        self.output_method_start(self.class_name, params, '', access='public')
        self.writeln('    : ', self.base_class_name, '(')
        # Filter out empty params, to make backward Euler happy
        base_class_params = filter(None, map(str, base_class_params))
        for i, param in enumerate(base_class_params):
            if i == len(base_class_params)-1: comma = ')'
            else: comma = ','
            self.writeln(param, comma, indent_offset=3)
        self.open_block()
        self.output_comment('Time units: ', self.free_vars[0].units, '\n')
        self.writeln('this->mpSystemInfo = OdeSystemInformation<',
                     self.class_name, '>::Instance();')
        if self.v_index == -1 and self.v_variable:
            self.writeln('this->mVoltageIndex = GetAnyVariableIndex("',
                         self.var_display_name(self.v_variable), '");')
        if self.config.options.include_dt_in_tables:
            self.writeln(self.lt_class_name, '::Instance()->SetTimestep(mDt);')
        self.writeln('Init();\n')
        
        #1861 - Rush-Larsen
        if self.options.rush_larsen and not self.doc._cml_rush_larsen:
            self.writeln('WARNING("No elligible gating variables found for this Rush-Larsen cell model; using normal forward Euler.");')
        
        #1463 - default cellML stimulus
        if self.has_default_stimulus:
            self.output_comment('We have a default stimulus specified in the CellML file metadata')
            self.writeln('this->mHasDefaultStimulusFromCellML = true', self.STMT_END)
            
        #1464 - cleverer modifiers...
        if self.use_modifiers and self.modifier_vars:
            self.output_comment('These will get initialised to DummyModifiers in the base class method.')
            for var in self.modifier_vars:
                self.writeln('this->AddModifier("' + var.oxmeta_name + '",')
                self.writeln('                  mp_' + var.oxmeta_name + '_modifier)', self.STMT_END)        
        
        #666 - initialise parameters
        for var in self.cell_parameters:
            if var.get_type() == VarTypes.Constant:
                self.writeln(self.vector_index('this->mParameters', var._cml_param_index),
                             self.EQ_ASSIGN, var.initial_value, self.STMT_END, ' ',
                             self.COMMENT_START, var.fullname(), ' [', var.units, ']')
        #1354 - specify protocol outputs
        if self.use_protocol:
            outputs = cellml_metadata.find_variables(self.model,
                                                     ('pycml:output-variable', NSS['pycml']),
                                                     'yes')
            def write_output_info(output):
                if output.get_type() in [VarTypes.Free, VarTypes.Unknown]:
                    self.writeln('UNSIGNED_UNSET, FREE', indent=False, nl=False)
                elif output.get_type() == VarTypes.State:
                    self.writeln(self.state_vars.index(output), ', STATE', indent=False, nl=False)
                elif output.is_derived_quantity:
                    self.writeln(self.derived_quantities.index(output), ', DERIVED', indent=False, nl=False)
                elif output.is_modifiable_parameter:
                    self.writeln(self.cell_parameters.index(output), ', PARAMETER', indent=False, nl=False)
                else:
                    raise ValueError('Unexpected protocol output: ' + str(output))
            if outputs:
                outputs.sort(key=lambda v: self.var_display_name(v))
                self.output_comment('Protocol outputs')
                self.writeln('this->mOutputsInfo.resize(', len(outputs), ');')
                for i, output in enumerate(outputs):
                    self.writeln('this->mOutputsInfo[', i, ']', self.EQ_ASSIGN,
                                 'std::make_pair(', nl=False)
                    write_output_info(output)
                    self.writeln(')', self.STMT_END, indent=False)
                self.writeln()
            outputs = set(outputs)
            #1925 - outputs that are vectors
            prop = ('pycml:output-vector', NSS['pycml'])
            vector_names = set(cellml_metadata.get_targets(self.model, None,
                                                           cellml_metadata.create_rdf_node(prop)))
            self.writeln('this->mVectorOutputsInfo.resize(', len(vector_names), ');')
            self.writeln('this->mVectorOutputNames.resize(', len(vector_names), ');')
            for i, name in enumerate(sorted(vector_names)):
                self.writeln('this->mVectorOutputNames[', i, ']', self.EQ_ASSIGN, '"', name, '"', self.STMT_END)
                vector_outputs = cellml_metadata.find_variables(self.model, prop, name)
                assert len(vector_outputs) > 0
                if name == 'state_variable':
                    # Special case to ensure the ordering as an output matches the state vector in the ODE system
                    def get_state_index(v):
                        """Find the index of the state variable corresponding to this variable, which may be units converted."""
                        v = v.get_source_variable(recurse=True)
                        if v.get_type() is VarTypes.Computed:
                            v = v.get_dependencies()[0].get_dependencies()[0]
                        return self.state_vars.index(v)
                    vector_outputs.sort(key=get_state_index)
                else:
                    vector_outputs.sort(key=lambda v: self.var_display_name(v))
                self.writeln('this->mVectorOutputsInfo[', i, '].resize(', len(vector_outputs), ');')
                for j, output in enumerate(vector_outputs):
                    self.writeln('this->mVectorOutputsInfo[', i, '][', j, ']', self.EQ_ASSIGN,
                                 'std::make_pair(', nl=False)
                    write_output_info(output)
                    self.writeln(')', self.STMT_END, indent=False)
                self.writeln()
                outputs.update(vector_outputs)
            #1910 - SED-ML name mappings
            prop = ('pycml:alias', NSS['pycml'])
            aliased_vars = cellml_metadata.find_variables(self.model, prop, None)
            prop = cellml_metadata.create_rdf_node(prop)
            for var in aliased_vars:
                assert var in outputs
                source = cellml_metadata.create_rdf_node(fragment_id=var.cmeta_id)
                for alias in cellml_metadata.get_targets(self.model, source, prop):
                    name = self.var_display_name(var)
                    self.writeln('this->mNameMap["', alias, '"] = "', name, '";')
            #2178 - set up model outputs environment from above info
            self.writeln()
            self.writeln('ProcessOutputsInfo();')
            self.writeln()
            #2428 - also record protocol inputs
            inputs = cellml_metadata.find_variables(self.model, ('pycml:input-variable', NSS['pycml']), 'yes')
            if inputs:
                inputs.sort(key=lambda v: self.var_display_name(v))
                self.writeln('this->mInputNames.reserve(', len(inputs), ');')
                for input in inputs:
                    self.writeln('this->mInputNames.push_back("', self.var_display_name(input), '");')

        # Lookup table generation, if not in a singleton
        if self.use_lookup_tables and not self.separate_lut_class:
            self.output_lut_generation()
        self.output_extra_constructor_content()
        self.close_block()
        return
    
    def output_extra_constructor_content(self):
        """Hook for subclasses to add further content to the constructor."""
        pass

    # Methods for interpolating on data files, used by Functional Curation

    def output_data_tables(self):
        """Output the data for interpolated lookups by output_special."""
        if hasattr(self.model, '_cml_interp_exprs'):
            # We do have data tables, but PE may have moved them, so rescan the whole model just in case!
            # This also ensures that self.model._cml_interp_exprs lists expressions in the order in which they (may) need to
            # be evaluated, in case one table uses the value looked up from another.
            self.model._cml_interp_exprs[:] = []
            def find_tables(expr):
                if hasattr(expr, '_cml_interp_data'):
                    self.model._cml_interp_exprs.append(expr)
                for elt in expr.xml_element_children():
                    find_tables(elt)
            for expr in self.model.get_assignments():
                if isinstance(expr, mathml_apply):
                    find_tables(expr)
        for i, expr in enumerate(getattr(self.model, '_cml_interp_exprs', [])):
            name = expr._cml_interp_table_name = '_data_table_' + str(i)
            data = expr._cml_interp_data
            self.output_array_definition(name, data[1])
            # Set variable names for use in other data-table methods
            expr._cml_interp_table_index = '_data_table_index_' + str(i)
            expr._cml_interp_table_factor = '_data_table_factor_' + str(i)
            # Calculate and cache the inverse of the table step
            expr._cml_interp_step_inverse = "%.17g" % (1.0 / (data[0,1] - data[0,0]))

    def output_data_table_lookups(self, nodeset):
        """Output the index and factor generation code for any data tables used in the given nodeset.
        
        We may also use the equations within nodeset to calculate the value(s) of the independent variable(s),
        and return a set containing just those nodes so used, so we can avoid recalculating them later.
        """
        nodes_used = set()
        for expr in getattr(self.model, '_cml_interp_exprs', []):
            top_expr = expr
            piecewises = []
            while not (isinstance(top_expr, mathml_apply) and top_expr.is_top_level()):
                if isinstance(top_expr.xml_parent, mathml_piecewise):
                    piecewises.append((top_expr.xml_parent, top_expr))
                top_expr = top_expr.xml_parent
            if top_expr in nodeset:
                # Code to calculate the independent variable (it may have a units conversion for instance)
                indep_expr = expr.operands().next()
                vars_used = self._vars_in(indep_expr)
                for piecewise, child in piecewises:
                    # We may also need to calculate variables used in the bounds check below
                    for piece in getattr(piecewise, u'piece', []):
                        vars_used.update(self._vars_in(child_i(piece, 2)))
                indep_nodes_used = self.calculate_extended_dependencies(vars_used, prune=nodes_used)
                self.output_equations(indep_nodes_used)
                nodes_used.update(indep_nodes_used)
                indep_var = expr._cml_interp_table_index + '_raw'
                self.writeln(self.TYPE_CONST_DOUBLE, indep_var, self.EQ_ASSIGN, nl=False)
                self.output_expr(indep_expr, paren=False)
                self.writeln(self.STMT_END, indent=False)
                # Code to check we're within the table bounds.
                # We need to take account of any piecewise expressions we are within.
                data = expr._cml_interp_data
                def process_piecewise(piecewises):
                    if len(piecewises) == 0:
                        # Base case
                        self.write('(%s >= %.17g' % (indep_var, data[0,0]), self.LOGICAL_AND,
                                   '%s <= %.17g)' % (indep_var, data[0,-1]))
                    else:
                        piecewise, piece_or_otherwise = piecewises.pop()
                        needs_operator = False
                        for piece in getattr(piecewise, u'piece', []):
                            if needs_operator:
                                self.write(self.LOGICAL_OR)
                            if piece is piece_or_otherwise:
                                # This is the branch that the interpolate occurs in
                                self.write('(')
                                self.output_expr(child_i(piece, 2), paren=True)
                                self.write(self.LOGICAL_AND)
                                process_piecewise(piecewises)
                                self.write(')')
                            else:
                                self.output_expr(child_i(piece, 2), paren=True)
                            needs_operator = True
                        if hasattr(piecewise, u'otherwise'):
                            if needs_operator:
                                self.write(self.LOGICAL_OR)
                            if piecewise.otherwise is piece_or_otherwise:
                                process_piecewise(piecewises)
                            else:
                                self.write(self.LOGICAL_TRUE) # TODO: Think about this some more!
                self.writeln(self.ASSERT, '(', nl=False)
                process_piecewise(piecewises)
                self.writeln(')', self.STMT_END)
                # Output the index & factor generation code
                offset_over_step = "(%s - %.17g)*%s" % (indep_var, data[0,0], expr._cml_interp_step_inverse)
                self.writeln(self.TYPE_CONST_UNSIGNED, expr._cml_interp_table_index, self.EQ_ASSIGN,
                             self.fast_floor(offset_over_step), self.STMT_END)
                self.writeln(self.TYPE_CONST_DOUBLE, expr._cml_interp_table_factor, self.EQ_ASSIGN,
                             offset_over_step, ' - ', expr._cml_interp_table_index, self.STMT_END)
        return nodes_used

    def output_special(self, expr, paren):
        """Output a special-case operation, represented by a csymbol element.
        
        The only operation currently supported is a table lookup on hardcoded data.
        """
        assert expr.operator().definitionURL == u'https://chaste.cs.ox.ac.uk/nss/protocol/interp'
        self.open_paren(paren)
        self.write(expr._cml_interp_table_name, '[', expr._cml_interp_table_index, '] + ',
                   expr._cml_interp_table_factor, ' * (',
                   expr._cml_interp_table_name, '[1+'+expr._cml_interp_table_index, ']-',
                   expr._cml_interp_table_name, '[',expr._cml_interp_table_index, '])')
        self.close_paren(paren)

    # The following three will need overriding in non-C++ translators

    def output_array_definition(self, array_name, array_data):
        """Output code to create and fill a fixed-size 1d array."""
        self.writeln('const static double ', array_name, '[] = {', ', '.join(map(lambda f: "%.17g" % f, array_data)), '};')

    def fast_floor(self, arg):
        """Return code to compute the floor of an argument as an integer quickly, typically by casting."""
        return "(unsigned)(%s)" % arg

    # Normal lookup table methods

    def output_chaste_lut_methods(self):
        """
        Output lookup table declarations & methods, if not using a separate class,
        or output the method to get a pointer to the lookup table collection.
        """
        if self.use_lookup_tables:
            if self.separate_lut_class:
                self.output_method_start('GetLookupTableCollection', [], 'AbstractLookupTableCollection*')
                self.open_block()
                self.writeln('return ', self.lt_class_name, '::Instance();')
                self.close_block()
            else:
                self.send_main_output_to_subsidiary()
                self.output_lut_declarations()
                self.output_lut_row_lookup_memory()
                self.output_lut_methods()
                self.send_main_output_to_subsidiary(False)
    
    def lut_parameters(self, key):
        """Get the bounds and step size for a particular table.
        
        key should be a key into self.lookup_table_indices.
        Returns (min, max, step) suitable for putting in generated code.
        """
        if self.separate_lut_class:
            idx = self.doc.lookup_table_indexes[key]
            return map(lambda s: 'mTable%s[%s]' % (s, idx), ['Mins', 'Maxs', 'Steps', 'StepInverses'])
        else:
            return super(CellMLToChasteTranslator, self).lut_parameters(key)
    
    def output_lut_indexing_methods(self):
        """Output methods in the LT class for indexing the tables, and checking index bounds.
        
        These will be methods like
            const double * const IndexTable0(double index_var);
        if self.row_lookup_method, or like
            void IndexTable0(double index_var, unsigned& index, double& factor);
        otherwise, with
            bool CheckIndex0(double& index_var);
        for checking the bounds.
        """
        for key, idx in self.doc.lookup_table_indexes.iteritems():
            varname = self.code_name(key[-1])
            method_name = 'IndexTable' + str(idx)
            if self.row_lookup_method:
                method = 'const double * %s(double %s)' % (method_name, varname)
            else:
                factor = self.lut_factor(idx)
                idx_var = '_table_index_' + str(idx)
                if factor:
                    factor = ', double& ' + factor
                method = 'void %s(double %s, unsigned& %s%s)' % (method_name, varname, idx_var, factor)
            self.writeln(method)
            self.open_block()
            self.output_table_index_generation_code(key, idx, call_method=False)
            if self.row_lookup_method:
                self.writeln('return _lt_', idx, '_row;')
            self.close_block()
            # And check the indexes
            if self.config.options.check_lt_bounds:
                self.writeln('// LCOV_EXCL_START', indent=False)
                self.writeln('bool CheckIndex', idx, '(double& ', varname, ')')
                self.open_block()
                self.output_table_index_checking(key, idx, call_method=False)
                self.writeln('return _oob_', idx, self.STMT_END)
                self.close_block(blank_line=False)
                self.writeln('// LCOV_EXCL_STOP\n', indent=False)
    
    def output_table_index_checking(self, key, idx, call_method=True):
        """Override base class method to call the methods on the lookup table class if needed."""
        if self.separate_lut_class and call_method:
            if self.config.options.check_lt_bounds:
                var = key[-1]
                varname = self.code_name(var)
                self.writeln('const bool _oob_', idx, self.EQ_ASSIGN, self.lt_class_name,
                             '::Instance()->CheckIndex', idx, '(', varname, ')', self.STMT_END)
        else:
            super(CellMLToChasteTranslator, self).output_table_index_checking(key, idx)
    
    def output_table_index_generation_code(self, key, idx, call_method=True):
        """Override base class method to call the methods on the lookup table class if needed."""
        if self.separate_lut_class and call_method:
            var = key[-1]
            varname = self.code_name(var)
            method_name = self.lt_class_name + '::Instance()->IndexTable' + str(idx)
            if self.row_lookup_method:
                self.writeln('const double* const _lt_', idx, '_row = ', method_name, '(', varname, ');')
            else:
                factor = self.lut_factor(idx, include_comma=True)
                idx_var = '_table_index_' + str(idx)
                self.writeln(method_name, '(', varname, ', ', idx_var, factor, ');')
        else:
            super(CellMLToChasteTranslator, self).output_table_index_generation_code(key, idx)

    def output_lut_class(self):
        """Output a separate class for lookup tables.
        
        This will live entirely in the .cpp file."""
        # Lookup tables class
        self.writeln('class ', self.lt_class_name, ' : public AbstractLookupTableCollection')
        self.writeln('{')
        self.writeln('public:')
        self.set_indent(1)
        # Method to get the table instance object
        self.writeln('static ', self.lt_class_name, '* Instance()')
        self.open_block()
        self.writeln('if (mpInstance.get() == NULL)')
        self.writeln('{')
        self.writeln('mpInstance.reset(new ', self.lt_class_name, ');', indent_offset=1)
        self.writeln('}')
        self.writeln('return mpInstance.get();')
        self.close_block()
        # Method to free the table memory
        self.writeln('void FreeMemory()')
        self.open_block()
        self.output_lut_deletion()
        self.writeln('mNeedsRegeneration.assign(mNeedsRegeneration.size(), true);')
        self.close_block()
        # Table lookup methods
        self.output_lut_methods()
        self.output_lut_indexing_methods()
        # Destructor
        self.writeln('~', self.lt_class_name, '()')
        self.open_block()
        self.output_lut_deletion()
        self.close_block()
        # Make the class a singleton
        self.writeln('protected:', indent_level=0)
        self.writeln(self.lt_class_name, '(const ', self.lt_class_name, '&);')
        self.writeln(self.lt_class_name, '& operator= (const ', self.lt_class_name, '&);')
        # Constructor
        self.writeln(self.lt_class_name, '()')
        self.open_block()
        self.writeln('assert(mpInstance.get() == NULL);')
        if self.config.options.include_dt_in_tables:
            self.writeln('mDt = HeartConfig::Instance()->GetOdeTimeStep();')
            self.writeln('assert(mDt > 0.0);')
        num_indexes = len(self.doc.lookup_table_indexes)
        self.writeln('mKeyingVariableNames.resize(', num_indexes, ');')
        self.writeln('mNumberOfTables.resize(', num_indexes, ');')
        self.writeln('mTableMins.resize(', num_indexes, ');')
        self.writeln('mTableSteps.resize(', num_indexes, ');')
        self.writeln('mTableStepInverses.resize(', num_indexes, ');')
        self.writeln('mTableMaxs.resize(', num_indexes, ');')
        self.writeln('mNeedsRegeneration.resize(', num_indexes, ');')
        for key, idx in self.doc.lookup_table_indexes.iteritems():
            min, max, step, var = key
            num_tables = unicode(self.doc.lookup_tables_num_per_index[idx])
            self.writeln('mKeyingVariableNames[', idx, '] = "', self.var_display_name(var), '";')
            self.writeln('mNumberOfTables[', idx, '] = ', num_tables, self.STMT_END)
            self.writeln('mTableMins[', idx, '] = ', min, self.STMT_END)
            self.writeln('mTableSteps[', idx, '] = ', step, self.STMT_END)
            self.writeln('mTableStepInverses[', idx, '] = ', str(1/float(step)), self.STMT_END)
            self.writeln('mTableMaxs[', idx, '] = ', max, self.STMT_END)
            self.writeln('mNeedsRegeneration[', idx, '] = true;')
            self.writeln('_lookup_table_', idx, self.EQ_ASSIGN, 'NULL', self.STMT_END)
        self.writeln(self.lt_class_name, '::RegenerateTables();')
        self.close_block()
        # Table (re-)generation
        self.writeln('void RegenerateTables()')
        self.open_block()
        event_handler = 'AbstractLookupTableCollection::EventHandler::'
        self.writeln(event_handler, 'BeginEvent(', event_handler, 'GENERATE_TABLES);')
        if self.config.options.include_dt_in_tables:
            self.writeln(self.TYPE_CONST_DOUBLE, self.code_name(self.config.dt_variable), ' = mDt;')
            # Hack: avoid unused variable warning
            self.writeln('double _unused = ', self.code_name(self.config.dt_variable), ';')
            self.writeln('_unused = _unused;\n')
        for idx in self.doc.lookup_table_indexes.itervalues():
            self.writeln('if (mNeedsRegeneration[', idx, '])')
            self.open_block()
            self.output_lut_deletion(only_index=idx)
            self.output_lut_generation(only_index=idx)
            self.writeln('mNeedsRegeneration[', idx, '] = false;')
            self.close_block(blank_line=True)
        self.writeln(event_handler, 'EndEvent(', event_handler, 'GENERATE_TABLES);')
        self.close_block()
        # Private data
        self.writeln('private:', indent_level=0)
        self.writeln('/** The single instance of the class */')
        self.writeln('static std::shared_ptr<', self.lt_class_name, '> mpInstance;\n')
        if self.row_lookup_method:
            self.output_lut_row_lookup_memory()
        self.output_lut_declarations()
        # Close the class
        self.set_indent(0)
        self.writeln('};\n')
        # Define the instance pointer
        self.writeln('std::shared_ptr<', self.lt_class_name, '> ', self.lt_class_name, '::mpInstance;')
        self.writeln()
        return

    def output_state_assignments(self, exclude_nonlinear=False,
                                 assign_rY=True,
                                 nodeset=None,
                                 pointer=''):
        """Output statements extracting state variables from their vector.

        If exclude_nonlinear is set to true, state variables appearing
        in the nonlinear system will not be included.

        If nodeset is given, only state variables appearing in nodeset
        will be included.
        
        If pointer is given, then the state variables actually appear in the
        variable given by pointer, which is of type const std::vector<double>*.
        """
        used_vars = set()
        for var in self.state_vars:
            if ((not exclude_nonlinear or var not in self.nonlinear_system_vars)
                and (nodeset is None or var in nodeset)):
                used_vars.add(var)
        if assign_rY and used_vars:
            if pointer:
                self.output_comment('For state variable interpolation (SVI) we read in interpolated state variables,')
                self.output_comment('otherwise for ionic current interpolation (ICI) we use the state variables of this model (node).')
                if self.TYPE_VECTOR_REF == CellMLToChasteTranslator.TYPE_VECTOR_REF:
                    self.writeln('if (!%s) %s = &rGetStateVariables();' % (pointer, pointer))
                    self.writeln('const ', self.TYPE_VECTOR_REF, 'rY = *', pointer, self.STMT_END)
                else:
                    self.writeln(self.TYPE_VECTOR_REF, 'rY;')
                    self.writeln('bool made_new_cvode_vector = false;')
                    self.writeln('if (!%s)' % (pointer))
                    self.open_block()
                    self.writeln('rY = rGetStateVariables();')
                    self.close_block(False)
                    self.writeln('else')
                    self.open_block()
                    self.writeln('made_new_cvode_vector = true;')
                    self.writeln('rY = MakeNVector(*%s);' % (pointer))
                    self.close_block()
            else:
                self.writeln(self.TYPE_VECTOR_REF, 'rY = rGetStateVariables();')
        if self.options.protocol:
            low_prop = ('pycml:range-low', NSS['pycml'])
            high_prop = ('pycml:range-high', NSS['pycml'])
            def check_bound(prop, reln, var, value):
                prop_value = var.get_rdf_annotation(prop)
                if prop_value:
                    value = '(%s %s %s ? %s : %s)' % (value, reln, prop_value, prop_value, value)
                return value
        for i, var in enumerate(self.state_vars):
            if var in used_vars:
                if self.use_modifiers and var in self.modifier_vars:
                    value = self.modifier_call(var, self.vector_index('rY', i))
                else:
                    value = self.vector_index('rY', i)
                if self.options.protocol:
                    value = check_bound(low_prop, '<', var, value)
                    value = check_bound(high_prop, '>', var, value)
                #2116 - use supplied fixed voltage if we're clamping
                if var is self.v_variable:
                    value = '(mSetVoltageDerivativeToZero ? this->mFixedVoltage : %s)' % value
                self.writeln(self.TYPE_DOUBLE, self.code_name(var),
                             self.EQ_ASSIGN, value, self.STMT_END)
                self.writeln(self.COMMENT_START, 'Units: ', var.units,
                             '; Initial value: ',
                             getattr(var, u'initial_value', 'Unknown'))
                #621 TODO: maybe convert if state var dimensions include time
        self.writeln()
        return
    
    def modifier_call(self, var, current_value):
        """Return code for a call to a modifier function for an oxmeta-annotated variable.
        
        The modifier function takes 2 parameters: the current value of the variable,
        and the current time.  It returns a modified value for the variable.
        """
        return ('mp_' + var.oxmeta_name + '_modifier->Calc(' +
                current_value + ', ' + self.code_name(self.free_vars[0]) + ')')
    
    def vector_index(self, vector, i):
        """Return code for accessing the i'th index of vector."""
        return vector + '[' + str(i) + ']'
    
    def vector_create(self, vector, size):
        """Return code for creating a new vector with the given size."""
        return ''.join(map(str, [self.TYPE_VECTOR, vector, '(', size, ')', self.STMT_END]))
    
    def vector_initialise(self, vector, size):
        """Return code for creating an already-declared vector with the given size."""
        return ''.join(map(str, [vector, '.resize(', size, ')', self.STMT_END]))
    
    def output_nonlinear_state_assignments(self, nodeset=None):
        """Output assignments for nonlinear state variables."""
#        processed=[]
        for i, var in enumerate(self.nonlinear_system_vars):
            if not nodeset or var in nodeset and self.code_name(var):# and not self.code_name(var) in processed:
                self.writeln(self.TYPE_DOUBLE, self.code_name(var), self.EQ_ASSIGN,
                             self.vector_index('rCurrentGuess', i), self.STMT_END)
#                processed.append(self.code_name(var))
                #621 TODO: maybe convert if state var dimensions include time
        self.writeln()
        return
    
    def get_stimulus_assignment(self):
        """Return code for getting Chaste's stimulus current."""
        expr = self.doc._cml_config.i_stim_var
        output = self.code_name(expr) + self.EQ_ASSIGN
        get_stim = 'GetIntracellularAreaStimulus(' + self.code_name(self.free_vars[0]) + ')'
        if self.doc._cml_config.i_stim_negated:
            get_stim = '-' + get_stim
        return output + get_stim + self.STMT_END

    def output_equations(self, nodeset, zero_stimulus=False):
        """Output the mathematics described by nodeset.

        nodeset represents a subset of the assignments in the model.
        Output assignments in the order given by a topological sort,
        but only include those in nodeset.
        """
        # Special case for the stimulus current
        if self.doc._cml_config.i_stim_var in nodeset:
            if zero_stimulus:
                i_stim = self.doc._cml_config.i_stim_var
                stim_assignment = self.code_name(i_stim) + self.EQ_ASSIGN + '0.0' + self.STMT_END
            else:
                stim_assignment = self.get_stimulus_assignment()
        for expr in (e for e in self.model.get_assignments() if e in nodeset):
            # Special-case the stimulus current
            if self.use_chaste_stimulus or zero_stimulus:
                if isinstance(expr, cellml_variable) and expr is self.doc._cml_config.i_stim_var:
                    self.writeln(self.TYPE_CONST_DOUBLE, stim_assignment)
                elif not (isinstance(expr, mathml_apply) and
                          isinstance(expr.operator(), mathml_eq) and
                          isinstance(expr.eq.lhs, mathml_ci) and
                          expr.eq.lhs.variable is self.doc._cml_config.i_stim_var):
                    self.output_assignment(expr)
            else:
                self.output_assignment(expr)
        return

    def output_assignment(self, expr):
        """Output an assignment statement.

        Has overrides for various special cases.
        """
        clear_type = False
        writing_data_clamp_current = False
        # Figure out what is being assigned to
        if isinstance(expr, cellml_variable):
            assigned_var = expr
        else:
            if expr.eq.lhs.localName == 'ci':
                assigned_var = expr.eq.lhs.variable
                if assigned_var is self.config.i_data_clamp_current:
                    writing_data_clamp_current = True
                    self.output_comment('Special handling of data clamp current here (see #2708)')
                    self.output_comment('(we want to save expense of calling the interpolation method if possible.)')
                    self.writeln(self.TYPE_DOUBLE, self.code_name(assigned_var), self.EQ_ASSIGN, '0.0' , self.STMT_END)
                    self.writeln('if (mDataClampIsOn)')
                    self.open_block()
                    clear_type = True
            else:
                assigned_var = None # We don't store derivatives as members
                #907: Check if this is the derivative of the transmembrane potential
                if not self.use_backward_euler and expr.eq.lhs.diff.dependent_variable == self.v_variable:
                    clear_type = True
                    
        # Parameters don't need assigning
        has_modifier = self.use_modifiers and getattr(assigned_var, '_cml_has_modifier', False)
        if assigned_var in self.cell_parameters and not has_modifier:
            return
        # Is the variable declared elsewhere?
        if clear_type:
            self.TYPE_DOUBLE = self.TYPE_CONST_DOUBLE = ''
        elif getattr(assigned_var, '_cml_modifiable', False):
            # Override const-ness, e.g. for a lookup table index
            self.TYPE_CONST_DOUBLE = self.TYPE_DOUBLE
        if (assigned_var and self.use_modifiers and assigned_var in self.modifier_vars
            and assigned_var.get_type() != VarTypes.State):
            # "Constant" oxmeta-annotated parameters may be modified at run-time
            if has_modifier:
                # Turn off the modifier to figure out the base value
                assigned_var._cml_has_modifier = False
                rhs = self.code_name(assigned_var)
                assigned_var._cml_has_modifier = True
            else:
                self.capture_output()
                super(CellMLToChasteTranslator, self).output_assignment(expr)
                assignment = self.get_captured_output()
                eq_pos = assignment.find(self.EQ_ASSIGN)
                end_pos = assignment.find(self.STMT_END)
                rhs = assignment[eq_pos+len(self.EQ_ASSIGN):end_pos]
            if rhs:
                # If assigned_var is computed, it'll 'appear' twice - once with expr==assigned_var,
                # and once for the assignment mathml_apply.  The former will result in an empty rhs.
                self.writeln(self.TYPE_CONST_DOUBLE, self.code_name(assigned_var), self.EQ_ASSIGN,
                             self.modifier_call(assigned_var, rhs), self.STMT_END, nl=False)
                self.output_comment(assigned_var.units, indent=False, pad=True)
        else:
            super(CellMLToChasteTranslator, self).output_assignment(expr)
#        if assigned_var:
#            # Debug
#            self.writeln('EXCEPT_IF_NOT(!std::isinf(', self.code_name(assigned_var), '));')
#            self.writeln('EXCEPT_IF_NOT(!std::isnan(', self.code_name(assigned_var), '));')
        if clear_type:
            # Remove the instance attributes, thus reverting to the class members
            del self.TYPE_DOUBLE
            del self.TYPE_CONST_DOUBLE
        elif getattr(assigned_var, '_cml_modifiable', False):
            del self.TYPE_CONST_DOUBLE
            
        if writing_data_clamp_current:
            self.close_block(False)
            
        return

    def output_mathematics(self):
        """Output the mathematics in this model.

        When backward Euler is used, we do so in 5 methods:
         * UpdateTransmembranePotential  does a forward Euler step for V
         * ComputeOneStepExceptVoltage  co-ordinates a backward Euler step
         * ComputeResidual and ComputeJacobian are used in the Newton iteration
         * GetIIonic returns the total ionic current
        
        Rush-Larsen is implemented similarly, with:
         * EvaluateEquations  evaluate the model derivatives and alpha/beta terms
         * ComputeOneStepExceptVoltage  does a Rush-Larsen update for eligible variables,
           and a forward Euler step for other non-V state variables
        Generalised Rush-Larsen methods also have specialised handling; see the
        individual methods for details.

        For other solvers, only 2 methods are needed:
         * EvaluateYDerivatives computes the RHS of the ODE system
         * GetIIonic is as above
        
        Where derived-quantity annotations are present, we also generate a
        ComputeDerivedQuantities method.
        """
        self.output_get_i_ionic()
        if self.options.rush_larsen:
            self.output_rush_larsen_mathematics()
        elif self.use_backward_euler:
            self.output_backward_euler_mathematics()
        elif self.options.grl1:
            self.output_grl1_mathematics()
        elif self.options.grl2:
            self.output_grl2_mathematics()
        else:
            self.output_evaluate_y_derivatives()
        self.output_derived_quantities()
    
    def calculate_lookup_table_indices(self, nodeset, time_name=None):
        """Output the lookup table index calculations needed for the given equations, if tables are enabled.
        
        If time_name is given, it may be used in exception messages for tables out of bounds.
        Note that it is needed to be passed in, since GetIIonic does not know the time.
        
        Returns the subset of nodeset used in calculating the indices.
        """
        if self.use_lookup_tables:
            nodes_used = self.output_table_index_generation(time_name, nodeset=nodeset)
        else:
            nodes_used = set()
        return nodes_used

    def output_get_i_ionic(self):
        """Output the GetIIonic method."""
        use_modifiers = self.use_modifiers
        self.use_modifiers = False
        self.output_method_start('GetIIonic', ['const std::vector<double>* pStateVariables'],
                                 self.TYPE_DOUBLE, access='public', defaults=['NULL'])
        self.open_block()
        # Output mathematics to calculate ionic current, using solver_info.ionic_current.
        if (hasattr(self.model, u'solver_info') and hasattr(self.model.solver_info, u'ionic_current')):
            if not hasattr(self.model.solver_info.ionic_current, u'var'):
                raise ValueError('No ionic currents found; check your configuration file')
            nodes = map(lambda elt: self.varobj(unicode(elt)),
                        self.model.solver_info.ionic_current.var)
            # GetIIonic must not include the stimulus current
            i_stim = self.doc._cml_config.i_stim_var
            nodeset = self.calculate_extended_dependencies(nodes, prune_deps=[i_stim])
            #print map(lambda v: v.fullname(), nodes)
            #print filter(lambda p: p[2]>0, map(debugexpr, nodeset))
            # Output main part of maths
            self.output_state_assignments(nodeset=nodeset, pointer='pStateVariables')
            table_index_nodes_used = self.calculate_lookup_table_indices(nodeset)
            self.output_equations(nodeset - table_index_nodes_used, zero_stimulus=True)
            self.writeln()
            # Assign the total current to a temporary so we can check for NaN
            self.writeln(self.TYPE_CONST_DOUBLE, 'i_ionic', self.EQ_ASSIGN, nl=False)
            if self.doc._cml_config.i_ionic_negated:
                self.writeln('-(', nl=False, indent=False)
            plus = False
            for varelt in self.model.solver_info.ionic_current.var:
                if plus: self.write('+')
                else: plus = True
                self.output_variable(varelt)
            if self.doc._cml_config.i_ionic_negated:
                self.writeln(')', nl=False, indent=False)
            self.writeln(self.STMT_END, indent=False)
            if self.TYPE_VECTOR_REF == CellMLToCvodeTranslator.TYPE_VECTOR_REF:
                self.writeln('if (made_new_cvode_vector)')
                self.open_block()
                self.writeln('DeleteVector(rY);')
                self.close_block(False)
            self.writeln('EXCEPT_IF_NOT(!std::isnan(i_ionic));')
            self.writeln('return i_ionic', self.STMT_END)
        else:
            self.writeln('return 0.0;')
        self.close_block()
        self.use_modifiers = use_modifiers

    def output_evaluate_y_derivatives(self, method_name='EvaluateYDerivatives'):
        """Output the EvaluateYDerivatives method."""
        # Start code output
        self.output_method_start(method_name,
                                 [self.TYPE_DOUBLE + self.code_name(self.free_vars[0]),
                                  'const ' + self.TYPE_VECTOR_REF + 'rY',
                                  self.TYPE_VECTOR_REF + 'rDY'],
                                 'void', access='public')
        self.open_block()
        if not self.state_vars:
            # This isn't an ODE model!
            self.close_block()
            return
        self.output_comment('Inputs:')
        self.output_comment('Time units: ', self.free_vars[0].units)
        self.output_derivative_calculations(self.state_vars)
        # Assign to derivatives vector
        for i, var in enumerate(self.state_vars):
            self.writeln(self.vector_index('rDY', i), self.EQ_ASSIGN, self.code_name(var, True), self.STMT_END)
        self.close_block()
        
    def output_derivative_calculations(self, state_vars, assign_rY=False, extra_nodes=set(),
                                       extra_table_nodes=set()):
        """
        This is used by self.output_evaluate_y_derivatives and self.output_rush_larsen_mathematics
        to compute the derivatives (and any extra nodes, if given).  It contains the special logic
        to obey the mSetVoltageDerivativeToZero member variable in the generated code.
        Returns a nodeset containing the equations output.
        """
        # Work out what equations are needed to compute the derivatives
        derivs = set(map(lambda v: (v, self.free_vars[0]), state_vars))
        if self.v_variable in state_vars:
            dvdt = (self.v_variable, self.free_vars[0])
            derivs.remove(dvdt) #907: Consider dV/dt separately
        else:
            dvdt = None
        if self.use_chaste_stimulus:
            i_stim = [self.doc._cml_config.i_stim_var]
        else:
            i_stim = []
        nonv_nodeset = self.calculate_extended_dependencies(derivs|extra_nodes, prune_deps=i_stim)
        if dvdt:
            if self.use_data_clamp:
                prune = set([self.config.i_data_clamp_data]) | nonv_nodeset
            else:
                prune = nonv_nodeset
            v_nodeset = self.calculate_extended_dependencies([dvdt], prune=prune, prune_deps=i_stim)
        else:
            v_nodeset = set()
        # State variable inputs
        all_nodes = nonv_nodeset|v_nodeset
        self.output_state_assignments(assign_rY=assign_rY, nodeset=all_nodes)
        self.writeln()
        table_index_nodes_used = self.calculate_lookup_table_indices(all_nodes|extra_table_nodes, self.code_name(self.free_vars[0]))
        table_index_nodes_used.update(self.output_data_table_lookups(all_nodes - table_index_nodes_used))
        self.output_comment('Mathematics')
        #907: Declare dV/dt
        if dvdt:
            self.writeln(self.TYPE_DOUBLE, self.code_name(self.v_variable, ode=True), self.STMT_END)
        # Output mathematics required for non-dV/dt derivatives (which may include dV/dt)
        self.output_equations(nonv_nodeset - table_index_nodes_used)
        self.writeln()
        
        #907: Calculation of dV/dt
        if dvdt:
            self.writeln('if (mSetVoltageDerivativeToZero)')
            self.open_block()
            self.writeln(self.code_name(self.v_variable, ode=True), self.EQ_ASSIGN, '0.0', self.STMT_END)
            self.close_block(blank_line=False)
            self.writeln('else')
            self.open_block()
            self.output_equations(v_nodeset - table_index_nodes_used)
            self.close_block()
            
        return all_nodes | table_index_nodes_used

    def output_backward_euler_mathematics(self):
        """Output the mathematics methods used in a backward Euler cell.

        Outputs ComputeResidual, ComputeJacobian,
        UpdateTransmembranePotential and ComputeOneStepExceptVoltage.
        """
        dt_name = 'mDt'
        #model_dt = self.varobj(self.model.solver_info.dt)
        if self.nonlinear_system_size > 0:
            # Residual
            ##########
            argsize = '[' + str(self.nonlinear_system_size) + ']'
            self.output_method_start('ComputeResidual',
                                     [self.TYPE_DOUBLE + self.code_name(self.free_vars[0]),
                                      self.TYPE_CONST_DOUBLE + 'rCurrentGuess' + argsize,
                                      self.TYPE_DOUBLE + 'rResidual' + argsize],
                                     'void', access='public')
            self.open_block()
            # Output mathematics for computing du/dt for each nonlinear state var u
            nodes = map(lambda u: (u, self.free_vars[0]), self.nonlinear_system_vars)
            nodeset = self.calculate_extended_dependencies(nodes, prune_deps=[self.doc._cml_config.i_stim_var])
 
            
            self.output_state_assignments(exclude_nonlinear=True, nodeset=nodeset)
            self.writeln("//output_nonlinear_state_assignments")
            self.output_nonlinear_state_assignments(nodeset=nodeset)
            table_index_nodes_used = self.calculate_lookup_table_indices(nodeset, self.code_name(self.free_vars[0]))
            self.writeln("//output_equations")
            self.output_equations(nodeset - table_index_nodes_used)
            self.writeln()
            # Fill in residual
            for i, var in enumerate(self.state_vars):
                try:
                    j = self.nonlinear_system_vars.index(var)
                except ValueError:
                    j = -1
                if j != -1:
                    self.writeln('rResidual[', j, '] = rCurrentGuess[', j, '] - rY[', i, '] - ',
                                 dt_name, '*', self.code_name(var, ode=True), self.STMT_END)
            self.close_block()
            
            # Jacobian
            ##########
            self.output_method_start('ComputeJacobian',
                                     [self.TYPE_DOUBLE + self.code_name(self.free_vars[0]),
                                      self.TYPE_CONST_DOUBLE + 'rCurrentGuess' + argsize,
                                      self.TYPE_DOUBLE + 'rJacobian' + argsize + argsize],
                                     'void', access='public')
            self.open_block()
            # Mathematics that the Jacobian depends on
            used_vars = set()
            for entry in self.model.solver_info.jacobian.entry:
                used_vars.update(self._vars_in(entry.math))
            nodeset = self.calculate_extended_dependencies(used_vars, prune_deps=[self.doc._cml_config.i_stim_var])
            self.output_state_assignments(exclude_nonlinear=True, nodeset=nodeset)
            self.output_nonlinear_state_assignments(nodeset=nodeset)
            self.writeln(self.TYPE_CONST_DOUBLE, self.code_name(self.config.dt_variable), self.EQ_ASSIGN, dt_name, self.STMT_END, '\n')
            table_index_nodes_used = self.calculate_lookup_table_indices(nodeset|set(map(lambda e: e.math, self.model.solver_info.jacobian.entry)), self.code_name(self.free_vars[0]))
            self.output_equations(nodeset - table_index_nodes_used)
            self.writeln()
            # Jacobian entries
            for entry in self.model.solver_info.jacobian.entry:
                var_i, var_j = entry.var_i, entry.var_j
                i = self.nonlinear_system_vars.index(self.varobj(var_i))
                j = self.nonlinear_system_vars.index(self.varobj(var_j))
                self.writeln('rJacobian[', i, '][', j, '] = ', nl=False)
                entry_content = list(entry.math.xml_element_children())
                assert len(entry_content) == 1, "Malformed Jacobian matrix entry: " + entry.xml()
                self.output_expr(entry_content[0], False)
                self.writeln(self.STMT_END, indent=False)
#            self.output_comment('Debugging')
#            self.writeln('#ifndef NDEBUG', indent=False)
#            self.writeln('for (unsigned i=0; i<', len(self.nonlinear_system_vars), '; i++)')
#            self.writeln('for (unsigned j=0; j<', len(self.nonlinear_system_vars), '; j++)', indent_offset=1)
#            self.writeln('EXCEPT_IF_NOT(!std::isnan(rJacobian[i][j]));', indent_offset=2)
#            self.writeln('//DumpJacobianToFile(', self.code_name(self.free_vars[0]),
#                         ', rCurrentGuess, rJacobian, rY);')
#            self.writeln('#endif // NDEBUG', indent=False)
            self.close_block()
        # The other methods are protected
        self.writeln_hpp('protected:', indent_offset=-1)
        
        # UpdateTransmembranePotential
        ##############################
        self.output_method_start('UpdateTransmembranePotential',
                                 [self.TYPE_DOUBLE + self.code_name(self.free_vars[0])],
                                 'void', access='public')
        self.open_block()
        self.output_comment('Time units: ', self.free_vars[0].units)
        # Output mathematics to compute dV/dt
        nodes = [(self.state_vars[self.v_index], self.free_vars[0])]
        nodeset = self.calculate_extended_dependencies(nodes, prune_deps=[self.doc._cml_config.i_stim_var])
        self.output_state_assignments(nodeset=nodeset)
        table_index_nodes_used = self.calculate_lookup_table_indices(nodeset, self.code_name(self.free_vars[0]))
        self.output_equations(nodeset - table_index_nodes_used)
        # Update V
        self.writeln()
        self.writeln('rY[', self.v_index, '] += ', dt_name, '*',
                     self.code_name(self.state_vars[self.v_index], ode=True), self.STMT_END)
        self.close_block()

        # ComputeOneStepExceptVoltage
        #############################
        self.output_method_start('ComputeOneStepExceptVoltage',
                                 [self.TYPE_DOUBLE + self.code_name(self.free_vars[0])],
                                 'void', access='public')
        self.open_block()
        self.writeln(self.COMMENT_START, 'Time units: ',
                     self.free_vars[0].units)
        # Output mathematics to update linear state variables, using solver_info.linear_odes.
        # Also need to use output_equations for variables used in the update equations.
        linear_vars, update_eqns = [], {}
        used_vars = set() # NB: Also contains update equation if is a mathml_apply so table index generation works
        for u, t, update_eqn in SolverInfo(self.model).get_linearised_odes():
            assert t == self.free_vars[0]
            assert len(update_eqn) == 1
            update_eqn = update_eqn[0]
            linear_vars.append(u)
            update_eqns[id(u)] = update_eqn
            if not isinstance(update_eqn, mathml_cn): used_vars.add(update_eqn)
            used_vars.update(self._vars_in(update_eqn))
        # Output required equations for used variables
        nodeset = self.calculate_extended_dependencies(used_vars, prune_deps=[self.doc._cml_config.i_stim_var])
        self.output_state_assignments(nodeset=nodeset)
        if self.config.dt_variable in nodeset:
            self.writeln(self.TYPE_CONST_DOUBLE, self.code_name(self.config.dt_variable), self.EQ_ASSIGN,
                         dt_name, self.STMT_END, '\n')
        table_index_nodes_used = self.calculate_lookup_table_indices(nodeset, self.code_name(self.free_vars[0]))
        self.output_equations(nodeset - table_index_nodes_used)
        # Update state variables:
        #   rY[i] = (rY[i] + _g_j*dt) / (1 - _h_j*dt)
        self.writeln()
        linear_vars.sort(key=lambda v: v.fullname())
        for i, u in enumerate(linear_vars):
            j = self.state_vars.index(u)
            self.writeln('rY[', j, ']', self.EQ_ASSIGN, nl=False)
            self.output_expr(update_eqns[id(u)], False)
            self.writeln(self.STMT_END, indent=False)
        # Set up the Newton iteration, if needed
        self.writeln()
        if self.nonlinear_system_size > 0:
            self.writeln('double _guess[', self.nonlinear_system_size, '] = {', nl=False)
            comma = False
            idx_map = [0] * self.nonlinear_system_size
            for i, var in enumerate(self.state_vars):
                try:
                    j = self.nonlinear_system_vars.index(var)
                    idx_map[j] = i
                except ValueError:
                    pass
            for i in idx_map:
                if comma: self.write(',')
                else: comma = True
                self.write('rY[', i, ']')
            self.writeln('};', indent=False)
            # Solve
            CNS = 'CardiacNewtonSolver<%d,%s>' % (self.nonlinear_system_size, self.class_name)
            self.writeln(CNS, '* _p_solver = ', CNS, '::Instance();')
            self.writeln('_p_solver->Solve(*this, ', self.code_name(self.free_vars[0]), ', _guess);')
            # Update state
            for j, i in enumerate(idx_map):
                self.writeln('rY[', i, '] = _guess[', j, '];')
        self.close_block()
    
    def output_rush_larsen_mathematics(self):
        """Output the special methods needed for Rush-Larsen style cell models.
        
        We generate:
         * EvaluateEquations  evaluate the model derivatives and alpha/beta terms
         * ComputeOneStepExceptVoltage  does a Rush-Larsen update for eligible variables,
           and a forward Euler step for other non-V state variables
        """
        rl_vars = self.doc._cml_rush_larsen
        # EvaluateEquations
        ###################
        self.output_method_start('EvaluateEquations',
                                 [self.TYPE_DOUBLE + self.code_name(self.free_vars[0]),
                                  'std::vector<double> &rDY',
                                  'std::vector<double> &rAlphaOrTau',
                                  'std::vector<double> &rBetaOrInf'],
                                 'void', access='public')
        self.open_block()
        normal_vars = [v for v in self.state_vars if not v in rl_vars]
        nodes, table_nodes = set(), set()
        for _, alpha_or_tau, beta_or_inf, _ in rl_vars.itervalues():
            table_nodes.add(alpha_or_tau)
            nodes.update(self._vars_in(alpha_or_tau))
            table_nodes.add(beta_or_inf)
            nodes.update(self._vars_in(beta_or_inf))
        self.output_derivative_calculations(normal_vars, True, nodes, table_nodes)
        # Now assign input vectors
        for i, var in enumerate(self.state_vars):
            if var in rl_vars:
                # Fill in rAlphaOrTau & rBetaOrInf
                self.writeln(self.vector_index('rAlphaOrTau', i), self.EQ_ASSIGN, nl=False)
                self.output_expr(rl_vars[var][1], False)
                self.writeln(self.STMT_END, indent=False)
                self.writeln(self.vector_index('rBetaOrInf', i), self.EQ_ASSIGN, nl=False)
                self.output_expr(rl_vars[var][2], False)
                self.writeln(self.STMT_END, indent=False)
            else:
                # Fill in rDY
                self.writeln(self.vector_index('rDY', i), self.EQ_ASSIGN, self.code_name(var, True), self.STMT_END)
        self.close_block()
        
        # ComputeOneStepExceptVoltage
        #############################
        self.output_method_start('ComputeOneStepExceptVoltage',
                                 ['const std::vector<double> &rDY',
                                  'const std::vector<double> &rAlphaOrTau',
                                  'const std::vector<double> &rBetaOrInf'],
                                 'void', access='public')
        self.open_block()
        self.writeln('std::vector<double>& rY = rGetStateVariables();')
        for i, var in enumerate(self.state_vars):
            if var in rl_vars:
                # Rush-Larsen update
                conv = rl_vars[var][3] or ''
                if conv: conv = '*' + str(conv)
                if rl_vars[var][0] == 'ab':
                    # Alpha & beta formulation
                    self.open_block()
                    self.writeln(self.TYPE_CONST_DOUBLE, 'tau_inv = rAlphaOrTau[', i, '] + rBetaOrInf[', i, '];')
                    self.writeln(self.TYPE_CONST_DOUBLE, 'y_inf = rAlphaOrTau[', i, '] / tau_inv;')
                    self.writeln('rY[', i, '] = y_inf + (rY[', i, '] - y_inf)*exp(-mDt', conv, '*tau_inv);')
                    self.close_block(blank_line=False)
                else:
                    # Tau & inf formulation
                    self.writeln('rY[', i, '] = rBetaOrInf[', i, '] + (rY[', i, '] - rBetaOrInf[', i, '])',
                                 '*exp(-mDt', conv, '/rAlphaOrTau[', i, ']);')
            elif var is not self.v_variable:
                # Forward Euler update
                self.writeln('rY[', i, '] += mDt * rDY[', i, '];')
        self.close_block()
    
    #Megan E. Marsh, Raymond J. Spiteri 
    #Numerical Simulation Laboratory 
    #University of Saskatchewan 
    #December 2011 
    #Partial support provided by research grants from 
    #the National Science and Engineering Research 
    #Council (NSERC) of Canada and the MITACS/Mprime 
    #Canadian Network of Centres of Excellence.
    def output_derivative_calculations_grl(self, var, assign_rY=False, extra_nodes=set(), extra_table_nodes=set()):
        """This is used by self.output_grl?_mathematics to get equations for each variable separately.

        Returns a node set with the equations output.
        """
        # Work out what equations are needed to compute the derivative of var
        if var in self.state_vars:
            dvardt = (var, self.free_vars[0])
            var_nodeset = self.calculate_extended_dependencies([dvardt])
        else:
            var_nodeset = set()
        # State variable inputs
        self.output_state_assignments(nodeset=var_nodeset, assign_rY=assign_rY)
        self.writeln()
        table_index_nodes_used = self.calculate_lookup_table_indices(var_nodeset, self.code_name(self.free_vars[0]))
        self.output_comment('Mathematics')
        self.output_equations(var_nodeset - table_index_nodes_used)
        return var_nodeset | table_index_nodes_used

    def find_grl_partial_derivatives(self):
        """If we have analytic Jacobian information available from Maple, find the terms needed for GRL methods.

        This caches where the diagonal entries are in the matrix, indexed by the state variable objects currently in use,
        since the entries in the matrix may reference non-partially-evaluated variables.
        """
        if not hasattr(self, 'jacobian_diagonal'):
            self.jacobian_diagonal = {}
        if self.use_analytic_jacobian and not self.jacobian_diagonal:
            for entry in self.model.solver_info.jacobian.entry:
                if entry.var_i == entry.var_j:
                    # It's a diagonal entry
                    var = self.varobj(entry.var_i).get_source_variable(recurse=True)
                    assert var in self.state_vars, "Jacobian diagonal entry is not in the state vector: " + entry.xml()
                    entry_content = list(entry.math.xml_element_children())
                    assert len(entry_content) == 1, "Malformed Jacobian entry: " + entry.xml()
                    self.jacobian_diagonal[var] = entry_content[0]

    def output_grl_compute_partial(self, i, var):
        """Compute the partial derivative of f(var) wrt var, the i'th variable in the state vector.

        This uses an analytic Jacobian if available; otherwise it approximates using finite differences.
        """
        self.output_method_start('EvaluatePartialDerivative'+str(i),
                                 [self.TYPE_DOUBLE + self.code_name(self.free_vars[0]),
                                  'std::vector<double>& rY', 'double delta', 'bool forceNumerical'],
                                 'double', access='public', defaults=['', '', '', 'false'])
        self.open_block()
        self.writeln('double partialF;')
        if self.jacobian_diagonal:
            # Work out what equations are needed to compute the analytic derivative
            self.writeln('if (!forceNumerical && this->mUseAnalyticJacobian)')
            self.open_block()
            entry = self.jacobian_diagonal[var]
            nodeset = self.calculate_extended_dependencies(self._vars_in(entry))
            self.output_state_assignments(nodeset=nodeset, assign_rY=False)
            table_index_nodes_used = self.calculate_lookup_table_indices(nodeset|set([entry]), self.code_name(self.free_vars[0]))
            self.output_equations(nodeset)
            # Calculate the derivative
            self.writeln('partialF = ', nl=False)
            self.output_expr(entry, paren=False)
            self.writeln(self.STMT_END, indent=False)
            self.close_block(blank_line=False)
            self.writeln('else')
            self.open_block()
        # Numerical approximation
        self.writeln('const double y_save = rY[', i, '];')
        self.writeln('rY[', i, '] += delta;')
        self.writeln('const double temp = EvaluateYDerivative', i, '(', self.code_name(self.free_vars[0]), ', rY);')
        self.writeln('partialF = (temp-mEvalF[', i, '])/delta;')
        self.writeln('rY[', i, '] = y_save;')
        if self.jacobian_diagonal:
            self.close_block(blank_line=False)
        self.writeln('return partialF;')
        self.close_block()

    #Megan E. Marsh, Raymond J. Spiteri 
    #Numerical Simulation Laboratory 
    #University of Saskatchewan 
    #December 2011 
    #Partial support provided by research grants from 
    #the National Science and Engineering Research 
    #Council (NSERC) of Canada and the MITACS/Mprime 
    #Canadian Network of Centres of Excellence.
    def output_grl1_mathematics(self):
        """Output the special methods needed for GRL1 style cell models.

        We generate:
         * UpdateTransmembranePotential update V_m
         * ComputeOneStepExceptVoltage  does a GRL1 update for variables except voltage
         * EvaluateYDerivativeI for each variable I
        """
        self.find_grl_partial_derivatives()
        ########################################################UpdateTransmembranePotential
        self.output_method_start('UpdateTransmembranePotential',
                                 [self.TYPE_DOUBLE + self.code_name(self.free_vars[0])],
                                 'void', access='public')
        self.open_block()
        self.writeln('std::vector<double>& rY = rGetStateVariables();')
        self.writeln('unsigned v_index = GetVoltageIndex();')
        self.writeln('const double delta = 1e-8;')
        self.writeln()
        # Compute partial derivative of dV wrt V
        self.writeln(self.TYPE_DOUBLE, self.code_name(self.v_variable, ode=True), self.STMT_END)
        self.output_derivative_calculations_grl(self.v_variable)
        self.writeln()
        self.writeln('double evalF = ', self.code_name(self.v_variable, ode=True), self.STMT_END)
        self.writeln('mEvalF[', self.v_index, '] = ', self.code_name(self.v_variable, ode=True), self.STMT_END)
        self.writeln('double partialF = EvaluatePartialDerivative', self.v_index, '(', self.code_name(self.free_vars[0]), ', rY, delta, true);')
        self.writeln('if (fabs(partialF) < delta)')
        self.open_block()
        self.writeln('rY[v_index] += evalF*mDt;')
        self.close_block(False)
        self.writeln('else')
        self.open_block()
        self.writeln('rY[v_index] += (evalF/partialF)*(exp(partialF*mDt)-1.0);')
        self.close_block()
        self.close_block()

        #########################################################ComputeOneStepExceptVoltage
        self.output_method_start('ComputeOneStepExceptVoltage',
                                 [self.TYPE_DOUBLE + self.code_name(self.free_vars[0])],
                                 'void', access='public')
        self.open_block()
        # Set up variables
        self.writeln('std::vector<double>& rY = rGetStateVariables();')
        self.writeln('const double delta = 1e-8;')
        self.writeln()

        # Evaluate RHS of equations (except dV/dt)
        non_v_vars = self.state_vars[:]
        if self.v_variable in non_v_vars:
            non_v_vars.remove(self.v_variable)
        self.output_derivative_calculations(non_v_vars)

        # Compute partial derivatives (for non-V)
        for i, var in enumerate(self.state_vars):
            if var is not self.v_variable:
                self.writeln('mEvalF[', i, '] = ', self.code_name(var, ode=True), self.STMT_END)
                self.writeln('mPartialF[', i, '] = EvaluatePartialDerivative', i, '(', self.code_name(self.free_vars[0]), ', rY, delta);')

        # Do the GRL updates
        for i, var in enumerate(self.state_vars):
            if var is not self.v_variable:
                self.open_block()
                self.writeln('if (fabs(mPartialF[', i, ']) < delta)')
                self.open_block()
                self.writeln('rY[', i, '] += mDt*', self.code_name(var, True), ';')
                self.close_block(False)
                self.writeln('else')
                self.open_block()
                self.writeln('rY[', i, '] += (', self.code_name(var, True), '/mPartialF[', i, '])*(exp(mPartialF[', i, ']*mDt)-1.0);')
                self.close_block()
                self.close_block()
        self.close_block()

        #########################################################Evaluate each equation
        for i, var in enumerate(self.state_vars):
            self.output_method_start('EvaluateYDerivative'+str(i),
                                     [self.TYPE_DOUBLE + self.code_name(self.free_vars[0]),
                                      'std::vector<double>& rY'],
                                     'double', access='public')
            self.open_block()
            if var is self.v_variable:
                self.writeln(self.TYPE_DOUBLE, self.code_name(self.v_variable, ode=True), self.STMT_END)
            self.output_derivative_calculations_grl(var)
            self.writeln()
            self.writeln('return ', self.code_name(var, True), ';')
            self.close_block()

            self.output_grl_compute_partial(i, var)

    #Megan E. Marsh, Raymond J. Spiteri 
    #Numerical Simulation Laboratory 
    #University of Saskatchewan 
    #December 2011 
    #Partial support provided by research grants from 
    #the National Science and Engineering Research 
    #Council (NSERC) of Canada and the MITACS/Mprime 
    #Canadian Network of Centres of Excellence.     
    def output_grl2_mathematics(self):
        """Output the special methods needed for GRL2 style cell models.

        We generate:
         * Update TransmembranePotential update V_m
         * ComputeOneStepExceptVoltage  does a GRL2 update for variables except voltage
         * EvaluateYDerivativeI for each variable I
        """
        self.find_grl_partial_derivatives()
        ########################################################UpdateTransmembranePotential
        self.output_method_start('UpdateTransmembranePotential',
                                 [self.TYPE_DOUBLE + self.code_name(self.free_vars[0])],
                                 'void', access='public')
        self.open_block()
        self.writeln('std::vector<double>& rY = rGetStateVariables();')
        self.writeln('const unsigned v_index = GetVoltageIndex();')
        self.writeln('const double delta = 1e-8;')
        self.writeln('const double yinit = rY[v_index];')
        self.writeln()

        # Do the first half step
        self.writeln(self.TYPE_DOUBLE, self.code_name(self.v_variable, ode=True), self.STMT_END)
        self.output_derivative_calculations_grl(self.v_variable)
        self.writeln()
        self.writeln('double evalF = ', self.code_name(self.v_variable, ode=True), self.STMT_END)
        self.writeln('mEvalF[', self.v_index, '] = ', self.code_name(self.v_variable, ode=True), self.STMT_END)
        self.writeln('double partialF = EvaluatePartialDerivative', self.v_index, '(', self.code_name(self.free_vars[0]), ', rY, delta, true);')
        self.writeln('if (fabs(partialF) < delta)')
        self.open_block()
        self.writeln('rY[v_index] += 0.5*evalF*mDt;')
        self.close_block(False)
        self.writeln('else')
        self.open_block()
        self.writeln('rY[v_index] += (evalF/partialF)*(exp(partialF*0.5*mDt)-1.0);')
        self.close_block()

        # Do the second half step
        self.writeln('rY[v_index] = yinit;')
        self.writeln('evalF = EvaluateYDerivative', self.v_index, '(', self.code_name(self.free_vars[0]), ', rY);')
        self.writeln('mEvalF[', self.v_index, '] = evalF;')
        self.writeln('partialF = EvaluatePartialDerivative', self.v_index, '(', self.code_name(self.free_vars[0]), ', rY, delta, true);')
        self.writeln('if (fabs(partialF) < delta)')
        self.open_block()
        self.writeln('rY[v_index] = yinit + evalF*mDt;')
        self.close_block(False)
        self.writeln('else')
        self.open_block()
        self.writeln('rY[v_index] = yinit + (evalF/partialF)*(exp(partialF*mDt)-1.0);')
        self.close_block()
        self.close_block() # End method

        #########################################################ComputeOneStepExceptVoltage
        self.output_method_start('ComputeOneStepExceptVoltage',
                                 [self.TYPE_DOUBLE + self.code_name(self.free_vars[0])],
                                 'void', access='public')
        self.open_block()
        # Set up variables
        self.writeln('std::vector<double>& rY = rGetStateVariables();')
        self.writeln('const double delta=1e-8;')
        self.writeln('const unsigned size = GetNumberOfStateVariables();')
        self.writeln('mYInit = rY;')
        self.writeln('double y_save;')
        self.writeln()

        # Calculate partial derivatives
        self.output_derivative_calculations(self.state_vars)
        for i, var in enumerate(self.state_vars):
            self.writeln(self.vector_index('mEvalF', i), self.EQ_ASSIGN, self.code_name(var, True), self.STMT_END)
        self.writeln()
        for i, var in enumerate(self.state_vars):
            if var is not self.v_variable:
                self.writeln('mPartialF[', i, '] = EvaluatePartialDerivative', i, '(', self.code_name(self.free_vars[0]), ', rY, delta);')

        # Update all variables
        self.writeln('for (unsigned var=0; var<size; var++)')
        self.open_block()
        self.writeln('if (var == ', self.v_index, ') continue;')
        self.writeln('if (fabs(mPartialF[var]) < delta)')
        self.open_block()
        self.writeln('rY[var] = mYInit[var] + 0.5*mDt*mEvalF[var];')
        self.close_block(False)
        self.writeln('else')
        self.open_block()
        self.writeln('rY[var] = mYInit[var] + (mEvalF[var]/mPartialF[var])*(exp(mPartialF[var]*0.5*mDt)-1.0);')
        self.close_block()
        self.close_block()
        self.writeln()

        # Determine new partial derivatives
        for i, var in enumerate(self.state_vars):
            if var is not self.v_variable:
                self.writeln()
                self.writeln('y_save = rY[', i, '];')
                self.writeln('rY[', i, '] = mYInit[', i, '];')
                self.writeln('mEvalF[', i, '] = EvaluateYDerivative', i, '(', self.code_name(self.free_vars[0]), ', rY);')
                self.writeln('mPartialF[', i, '] = EvaluatePartialDerivative', i, '(', self.code_name(self.free_vars[0]), ', rY, delta);')
                self.writeln('rY[', i, '] = y_save;')

        # Update all variables
        self.writeln('for (unsigned var=0; var<size; var++)')
        self.open_block()
        self.writeln('if (var == ', self.v_index, ') continue;')
        self.writeln('if (fabs(mPartialF[var]) < delta)')
        self.open_block()
        self.writeln('rY[var] = mYInit[var] + mDt*mEvalF[var];')
        self.close_block(False)
        self.writeln('else')
        self.open_block()
        self.writeln('rY[var] = mYInit[var] + (mEvalF[var]/mPartialF[var])*(exp(mPartialF[var]*mDt)-1.0);')
        self.close_block()
        self.close_block()
        self.writeln()
        self.close_block() # End method

        #########################################################Evaluate each equation
        for i, var in enumerate(self.state_vars):
            self.output_method_start('EvaluateYDerivative'+str(i),
                                     [self.TYPE_DOUBLE + self.code_name(self.free_vars[0]),
                                      'std::vector<double>& rY'],
                                     'double', access='public')
            self.open_block()
            if var is self.v_variable:
                self.writeln(self.TYPE_DOUBLE, self.code_name(self.v_variable, ode=True), self.STMT_END)
            self.output_derivative_calculations_grl(var)
            self.writeln()
            self.writeln('return '+self.code_name(var, True)+';')
            self.close_block()

            self.output_grl_compute_partial(i, var)


    def output_model_attributes(self):
        """Output any named model attributes defined in metadata.

        Such attributes are given by compound RDF annotations:
          model --pycml:named-attribute--> bnode
          bnode --pycml:name--> Literal(Attribute name, string)
          bnode --pycml:value--> Literal(Attribute value, double)
        """
        model = self.model
        meta_id = model.cmeta_id
        attrs = []
        if meta_id:
            property = cellml_metadata.create_rdf_node(('pycml:named-attribute', NSS['pycml']))
            name_prop = cellml_metadata.create_rdf_node(('pycml:name', NSS['pycml']))
            value_prop = cellml_metadata.create_rdf_node(('pycml:value', NSS['pycml']))
            source = cellml_metadata.create_rdf_node(fragment_id=meta_id)
            attr_nodes = cellml_metadata.get_targets(model, source, property)
            for node in attr_nodes:
                name = cellml_metadata.get_target(model, node, name_prop)
                value = cellml_metadata.get_target(model, node, value_prop)
                attrs.append((name, value))
        for name, value in attrs:
            self.writeln('this->mAttributes["', name, '"] = ', value, ';')
        if attrs:
            self.writeln()

    def output_bottom_boilerplate(self):
        """Output bottom boilerplate.

        End class definition, output ODE system information (to .cpp) and
        serialization code (to .hpp), and end the file.
        """
        # End main class
        self.set_indent(offset=-1)
        self.writeln_hpp('};\n\n')
        # ODE system information
        self.writeln('template<>')
        self.writeln('void OdeSystemInformation<', self.class_name,
                     '>::Initialise(void)')
        self.open_block()
        self.writeln('this->mSystemName', self.EQ_ASSIGN, '"', self.model.name, '"', self.STMT_END)
        self.writeln('this->mFreeVariableName', self.EQ_ASSIGN,
                     '"', self.var_display_name(self.free_vars[0]), '"', self.STMT_END)
        self.writeln('this->mFreeVariableUnits', self.EQ_ASSIGN,
                     '"', self.free_vars[0].units, '"', self.STMT_END)
        self.writeln()
        def output_var(vector, var):
            self.writeln('this->m', vector, 'Names.push_back("', self.var_display_name(var), '");')
            self.writeln('this->m', vector, 'Units.push_back("', var.units, '");')
        for i, var in enumerate(self.state_vars):
            self.output_comment('rY[', str(i) ,']:')
            output_var('Variable', var)
            init_val = getattr(var, u'initial_value', None)
            if init_val is None:
                init_comm = ' // Value not given in model'
                # Don't want compiler error, but shouldn't be a real number
                init_val = self.NOT_A_NUMBER
            else:
                init_comm = ''
            self.writeln('this->mInitialConditions.push_back(', init_val, ');',
                       init_comm, '\n')
        # Model parameters
        for i,var in enumerate(self.cell_parameters):
            if var.get_type() == VarTypes.Constant:
                self.output_comment('mParameters[', str(i), ']:')
                output_var('Parameter', var)
                self.writeln()
        # Derived quantities
        for i,var in enumerate(self.derived_quantities):
            self.output_comment('Derived Quantity index [', str(i), ']:')
            output_var('DerivedQuantity', var)
            self.writeln()
        self.output_model_attributes()
        self.writeln('this->mInitialised = true;')
        self.close_block()
        self.writeln()
        # Serialization
        if self.include_serialization:
            self.output_comment('Needs to be included last', subsidiary=True)
            self.writeln_hpp('#include "SerializationExportWrapper.hpp"')
            self.writeln_hpp('CHASTE_CLASS_EXPORT(', self.class_name, ')')
            self.output_comment('Serialization for Boost >= 1.36')
            self.writeln('#include "SerializationExportWrapperForCpp.hpp"')
            self.writeln('CHASTE_CLASS_EXPORT(', self.class_name, ')')
            self.writeln_hpp()
            self.writeln_hpp('namespace boost')
            self.open_block(subsidiary=True)
            self.writeln_hpp('namespace serialization')
            self.open_block(subsidiary=True)
            # Save
            self.writeln_hpp('template<class Archive>')
            self.writeln_hpp('inline void save_construct_data(')
            self.writeln_hpp('Archive & ar, const ', self.class_name,
                             ' * t, const unsigned int fileVersion)',
                             indent_offset=1)
            self.open_block(subsidiary=True)
            self.writeln_hpp('const boost::shared_ptr<AbstractIvpOdeSolver> p_solver = t->GetSolver();')
            self.writeln_hpp('const boost::shared_ptr<AbstractStimulusFunction> p_stimulus = t->GetStimulusFunction();')
            self.writeln_hpp('ar << p_solver;')
            self.writeln_hpp('ar << p_stimulus;')
            self.close_block(subsidiary=True)
            # Load
            self.writeln_hpp('template<class Archive>')
            self.writeln_hpp('inline void load_construct_data(')
            self.writeln_hpp('Archive & ar, ', self.class_name,
                             ' * t, const unsigned int fileVersion)',
                             indent_offset=1)
            self.open_block(subsidiary=True)
            self.writeln_hpp('boost::shared_ptr<AbstractIvpOdeSolver> p_solver;')
            self.writeln_hpp('boost::shared_ptr<AbstractStimulusFunction> p_stimulus;')
            self.writeln_hpp('ar >> p_solver;')
            self.writeln_hpp('ar >> p_stimulus;')
            self.writeln_hpp('::new(t)', self.class_name, '(p_solver, p_stimulus);')
            self.close_block(subsidiary=True)
            self.close_block(subsidiary=True)
            self.close_block(subsidiary=True)
        if self.dynamically_loadable:
            # Write the C function to create instances of this cell model
            self.writeln('extern "C"')
            self.open_block()
            self.writeln('AbstractCardiacCellInterface* MakeCardiacCell(')
            self.writeln('boost::shared_ptr<AbstractIvpOdeSolver> pSolver,', indent_offset=2)
            self.writeln('boost::shared_ptr<AbstractStimulusFunction> pStimulus)', indent_offset=2)
            self.open_block()
            self.writeln('return new ', self.class_name, '(pSolver, pStimulus);')
            self.close_block()
            self.close_block()
        # End file
        self.writeln_hpp('#endif // ', self.include_guard)
        return

    def output_lhs(self, expr):
        """Output the left hand side of an assignment expression."""
        if expr.localName == 'ci':
            self.output_variable(expr)
        elif expr.operator().localName == 'diff':
            ci_elt = expr.operands().next()
            self.output_variable(ci_elt, ode=True)
        return

    def output_variable(self, ci_elt, ode=False):
        """Output a ci element, i.e. a variable lookup."""
        if hasattr(ci_elt, '_cml_variable') and ci_elt._cml_variable:
            self.write(self.code_name(ci_elt.variable, ode=ode))
        else:
            # This ci element doesn't have all the extra annotations.  It is a fully
            # qualified name though.  This is typically because PE has been done.
            prefix = ['var_', 'd_dt_'][ode]
            varname = unicode(ci_elt)
            try:
                var = self.varobj(varname)
            except KeyError:
                var = None
            if var:
                self.write(self.code_name(var, ode=ode))
            else:
                # Assume it's a suitable name
                self.write(prefix + varname)
        return

    def output_function(self, func_name, args, *posargs, **kwargs):
        """Override base class method for special case of abs with 2 arguments.
        
        This comes from Maple's Jacobians, and should generate signum of the second argument.
        """
        args = list(args)
        if func_name == 'fabs' and len(args) == 2:
            super(CellMLToChasteTranslator, self).output_function('Signum', [args[1]], *posargs, **kwargs)
        else:
            super(CellMLToChasteTranslator, self).output_function(func_name, args, *posargs, **kwargs)
    
    @staticmethod
    def get_current_units_options(model):
        """
        Return a list of units objects that give the possibilities for the dimensions
        of transmembrane ionic currents.
        """
        chaste_units = cellml_units.create_new(
            model, 'uA_per_cm2',
            [{'units': 'ampere', 'prefix': 'micro'},
             {'units': 'metre', 'prefix': 'centi', 'exponent': '-2'}])
        microamps = cellml_units.create_new(model, u'microamps',
                                            [{'units':'ampere', 'prefix':'micro'}])
        A_per_F = cellml_units.create_new(model, 'A_per_F',
                                          [{'units': 'ampere'},
                                           {'units': 'farad', 'exponent': '-1'}])
        return [chaste_units, microamps, A_per_F]

    # Name in CellML for the variable representing Chaste's membrane capacitance
    MEMBRANE_CAPACITANCE_NAME = u'chaste_membrane_capacitance'
    
    # Name of the component added to interface the model to Chaste
    INTERFACE_COMPONENT_NAME = u'chaste_interface'

    @staticmethod
    def add_special_conversions(converter, comp):
        """Add special units conversions for ionic currents.
        
        Adds conversions for the two other common conventions to/from the units expected by Chaste,
        uA/cm^2.  The cases are:
        
        1. Current in amps/farads.
           In this case we convert to uA/uF then multiply by Chaste's value
           for the membrane capacitance (in uF/cm^2).
        2. Current in amps, capacitance in farads.
           We assume the cell model conceptually represents a cell, and hence
           that its membrane capacitance is supposed to represent the same
           thing as Chaste's.  Thus convert current to uA, capacitance to uF,
           and return current/capacitance * Chaste's capacitance.
        
        comp is a component to which we should add any necessary variables, i.e. Chaste's capacitance.
        """
        klass = CellMLToChasteTranslator
        model = converter.model
        # Variables needed by some conversions
        model_Cm = model.get_config('Cm_variable')
        uF_per_cm2 = cellml_units.create_new(model, 'uF_per_cm2',
                                             [{'units': 'farad', 'prefix': 'micro'},
                                              {'units': 'metre', 'prefix': 'centi', 'exponent': '-2'}])
        Chaste_Cm = converter.add_variable(comp, klass.MEMBRANE_CAPACITANCE_NAME, uF_per_cm2)
        model._cml_Chaste_Cm = Chaste_Cm # Record for use in code_name
        # Add the conversions
        chaste_units, microamps, A_per_F = klass.get_current_units_options(model)
        converter.add_special_conversion(A_per_F, chaste_units,
                                         lambda expr: converter.times_rhs_by(expr, Chaste_Cm))
        converter.add_special_conversion(chaste_units, A_per_F,
                                         lambda expr: converter.divide_rhs_by(expr, Chaste_Cm))
        if model_Cm:
            converter.add_special_conversion(microamps, chaste_units,
                    lambda expr: converter.times_rhs_by(converter.divide_rhs_by(expr, model_Cm),
                                                        Chaste_Cm))
            converter.add_special_conversion(chaste_units, microamps,
                    lambda expr: converter.divide_rhs_by(converter.times_rhs_by(expr, model_Cm),
                                                         Chaste_Cm))
        
    @staticmethod
    def generate_interface(doc, solver_info):
        """Generate an interface component connecting the model to Chaste.
        
        On return from this method, Chaste code will only need to interact with variables in
        the new interface component.  It will contain the transmembrane potential, the ionic
        and stimulus currents, the simulation time, and the derivatives.
        
        It may also contain other variables depending on the model, for example the intracellular
        calcium concentration (if annotated), modifiable parameters, and derived quantities.
        
        If the --convert-interfaces option has been supplied, units conversion will then be
        performed on this component, ensuring that all these variables are in the units expected
        by Chaste and linked by suitable conversions to the rest of the model.
        
        Note that if partial evaluation is then performed, the model will be collapsed into a
        single component.  However, the interface will still be preserved in the correct units.
        """
        model = doc.model
        config = doc._cml_config
        klass = CellMLToChasteTranslator
        # Create Chaste units definitions
        ms = cellml_units.create_new(model, 'millisecond',
                                     [{'units': 'second', 'prefix': 'milli'}])
        mV = cellml_units.create_new(model, 'millivolt',
                                     [{'units': 'volt', 'prefix': 'milli'}])

        milliMolar = cellml_units.create_new(model,'millimolar', 
                            [{'units': 'mole', 'prefix': 'milli'},
                             {'units': 'litre', 'exponent': '-1'}])


        current_units, microamps = klass.get_current_units_options(model)[0:2]
        # The interface generator
        generator = processors.InterfaceGenerator(model, name=klass.INTERFACE_COMPONENT_NAME)


        iface_comp = generator.get_interface_component()
        # In case we need to convert initial values, we create the units converter here
        if config.options.convert_interfaces:
            warn_only = not config.options.fully_automatic and config.options.warn_on_units_errors
            notifier = NotifyHandler(level=logging.WARNING)
            logging.getLogger('units-converter').addHandler(notifier)
            converter = processors.UnitsConverter(model, warn_only, show_xml_context_only=True)
            klass.add_special_conversions(converter, iface_comp)
            generator.set_units_converter(converter)
        # And now specify the interface
        t = model.find_free_vars()[0]
        if not ms.dimensionally_equivalent(t.get_units()):
            # Oops!
            raise TranslationError('Time does not have dimensions of time')
        generator.add_input(t, ms)

        if doc.model.get_option('backward_euler'):
            # Backward Euler code generation requires access to the time step
            model_dt = solver_info.create_dt(generator, t.component, t.get_units())
            config.dt_variable = generator.add_input(model_dt, ms)
            config.dt_variable.set_pe_keep(True)
        elif doc.model.get_option('maple_output'):
            # CVODE Jacobians need to be able to scale for time too
            fake_dt = generator.add_variable(t.component, 'fake_dt', ms, initial_value='1.0')
            fake_dt._set_type(VarTypes.Constant)
            config.dt_variable = generator.add_input(fake_dt, t.get_units())
            config.dt_variable.set_is_modifiable_parameter(False)
            config.dt_variable.set_pe_keep(True)

        if config.options.use_chaste_stimulus and config.i_stim_var:
            # We need to make it a constant so add_input doesn't complain, then make it computed
            # again so that exposing metadata-annotated variables doesn't make it a parameter!
            generator.make_var_constant(config.i_stim_var, 0)
            config.i_stim_var = generator.add_input(config.i_stim_var, current_units,
                                                    annotate=False, convert_initial_value=False)
            generator.make_var_computed_constant(config.i_stim_var, 0)
            # Also convert variables that make up the default stimulus
            # Note: we vary in/out-put primarily to test units conversion of initial values
            def add_oxmeta_ioput(oxmeta_name, units, inout):
                var = doc.model.get_variable_by_oxmeta_name(oxmeta_name, throw=False)
                if var:
                    meth = getattr(generator, 'add_%sput' % inout)
                    newvar = meth(var, units, annotate=False)
                    newvar.set_pe_keep(True)
            for n in ['duration', 'period', 'offset', 'end']:
                add_oxmeta_ioput('membrane_stimulus_current_'+n, ms, 'in')
            add_oxmeta_ioput('membrane_stimulus_current_amplitude', current_units, 'out')

        if config.V_variable:
            config.V_variable = generator.add_input(config.V_variable, mV)



        ionic_vars = config.i_ionic_vars
        if ionic_vars:
            i_ionic = generator.add_output_function('i_ionic', 'plus', ionic_vars, current_units)
            config.i_ionic_vars = [i_ionic]

        if doc.model.get_option('use_data_clamp'):
            assert config.V_variable and ionic_vars
            # Create g_clamp
            conductance_units = current_units.quotient(mV).simplify()
            i_data_clamp_conductance = generator.add_variable(iface_comp, 'membrane_data_clamp_current_conductance', conductance_units, initial_value='0.0')
            i_data_clamp_conductance._set_type(VarTypes.Constant)
            i_data_clamp_conductance.set_pe_keep(True) # This prevents it becoming 'chaste_interface__membrane_data_clamp_current_conductance'
            config.i_data_clamp_conductance = generator.add_input(i_data_clamp_conductance, conductance_units)
            # Create V_clamp
            data_var = config.i_data_clamp_data = generator.add_variable(iface_comp, 'experimental_data_voltage', mV, initial_value='0.0')
            data_var._set_type(VarTypes.Constant)
            data_var.set_pe_keep(True)
            data_var._cml_code_name = 'GetExperimentalVoltageAtTimeT(%(time)s)'
            # Create the current: I = g_clamp * (V - V_clamp)
            current_var = config.i_data_clamp_current = generator.add_variable(iface_comp, 'membrane_data_clamp_current', current_units)
            current_var._set_type(VarTypes.Computed)
            current_var.set_is_derived_quantity(True)
            sub = mathml_apply.create_new(model, u'minus', [config.V_variable.name, data_var.name])
            times = mathml_apply.create_new(model, u'times', [config.i_data_clamp_conductance.name, sub])
            assign = mathml_apply.create_new(model, u'eq', [current_var.name, times])
            generator.add_expr_to_comp(iface_comp, assign)
            # Make dV/dt use the new current
            def process_ci(elt):
                # Add reference to new current after first existing ionic current
                ref = mathml_ci.create_new(model, local_current_var.name)
                elt.xml_parent.xml_insert_after(elt, ref)
            if hasattr(ionic_vars[0], '_cml_ref_in_dvdt'):
                local_current_var = generator.connect_variables(current_var, (ionic_vars[0]._cml_ref_in_dvdt.component.name, current_var.name))
                process_ci(ionic_vars[0]._cml_ref_in_dvdt)
            else:
                dVdt = config.V_variable.get_all_expr_dependencies()[0]
                local_current_var = generator.connect_variables(current_var, (config.V_variable.component.name, current_var.name))
                def process_ci_elts(elt):
                    """Recursively process any ci elements in the tree rooted at elt."""
                    if isinstance(elt, mathml_ci):
                        if elt.variable is ionic_vars[0]:
                            process_ci(elt)
                    else:
                        for child in getattr(elt, 'xml_children', []):
                            process_ci_elts(child)
                process_ci_elts(dVdt)

        #Unit conversion for cytosolic_calcium_variable
        ##Try to check  if cytosolic_calcium_variable has a dimension of molar. If it fails its units may not be defined properly
        #Cannot just use generator.add_input as thsi may cause duplicates in BackwardsEuler odes
        try:
            if config.options.convert_interfaces and config.cytosolic_calcium_variable and milliMolar.dimensionally_equivalent(config.cytosolic_calcium_variable.get_units()):
                if config.cytosolic_calcium_variable.get_type()==VarTypes.Computed:
                    config.cytosolic_calcium_variable = generator.add_output(config.cytosolic_calcium_variable, milliMolar)
                else:
                    config.cytosolic_calcium_variable = generator.add_input(config.cytosolic_calcium_variable, milliMolar)
        except AttributeError:
            DEBUG('generate_interface', "Model has no cytosolic_calcium_variable")
        except:
             DEBUG('generate_interface', "Unit conversion for cytosolic_calcium_variable failed.")

        # Finish up
        def errh(errors):
            raise TranslationError("Creation of Chaste interface component failed:\n  " + str(errors))
        generator.finalize(errh, check_units=False)
        # Apply units conversions just to the interface, if desired
        if config.options.convert_interfaces:
            converter.add_conversions_for_component(iface_comp)
            converter.finalize(errh, check_units=False)
            notifier.flush()
            logging.getLogger('units-converter').removeHandler(notifier)
            if notifier.messages:
                msg = 'Problems occurred converting model variables to Chaste units.\n'
                if ionic_vars and ionic_vars[0].get_units().dimensionally_equivalent(microamps):
                    msg += 'To convert the ionic currents for this model, '\
                           'the model membrane capacitance needs to be identified.'
                if config.options.fully_automatic:
                    raise TranslationError(msg)
                else:
                    print >>sys.stderr, msg

class CellMLToCvodeTranslator(CellMLToChasteTranslator):
    """Translate a CellML model to C++ code for use with Chaste+CVODE."""

    # Type of (a reference to) the state variable vector
    TYPE_VECTOR = 'N_Vector '
    TYPE_VECTOR_REF = 'N_Vector ' # CVODE's vector is actually a pointer type

    def vector_index(self, vector, i):
        """Return code for accessing the i'th index of vector."""
        return 'NV_Ith_S(' + vector + ', ' + str(i) + ')'

    def vector_create(self, vector, size):
        """Return code for creating a new vector with the given size."""
        return ''.join(map(str, [self.TYPE_VECTOR, vector, self.EQ_ASSIGN,
                                 'N_VNew_Serial(', size, ')', self.STMT_END]))

    def vector_initialise(self, vector, size):
        """Return code for creating an already-declared vector with the given size."""
        return ''.join(map(str, [vector, self.EQ_ASSIGN, 'N_VNew_Serial(', size, ')', self.STMT_END]))

    def output_top_boilerplate(self):
        """Output top boilerplate code.

        This method outputs #includes, and the start of the cell class
        with constructor, destructor, and LT methods.
        """
        # CVODE is optional in Chaste
        self.writeln("#ifdef CHASTE_CVODE")
        self.writeln_hpp("#ifdef CHASTE_CVODE")

        self.include_serialization = True
        self.use_backward_euler = False
        
        if self.use_data_clamp:
            self.use_analytic_jacobian = False # Todo - data current not included in analytic jacobian yet.
            self.output_includes(base_class='AbstractCvodeCellWithDataClamp')
        else:
            self.use_analytic_jacobian = (self.model.get_option('maple_output') and hasattr(self.model.solver_info, u'jacobian'))
            self.output_includes(base_class='AbstractCvodeCell')
        
        # Separate class for lookup tables?
        if self.use_lookup_tables and self.separate_lut_class:
            self.output_lut_class()
        self.output_data_tables()
        # Start cell model class
        self.writeln_hpp('class ', self.class_name, self.class_inheritance)
        self.open_block(subsidiary=True)
        # Serialization
        self.output_serialize_method()
        # Parameter declarations, and set & get methods (#666)
        self.output_cell_parameters()
        # Constructor
        self.output_constructor(['boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver /* unused; should be empty */',
                                 'boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus'],
                                ['pOdeSolver', len(self.state_vars), self.unsigned_v_index, 'pIntracellularStimulus'])
        # Destructor
        self.output_method_start('~'+self.class_name, [], '', access='public')
        self.open_block()
        self.close_block()
        # Other declarations & methods
        self.output_chaste_lut_methods()
        self.output_verify_state_variables()

    def output_mathematics(self):
        """Output the mathematics in this model.
        
        Two methods are needed:
         * EvaluateYDerivatives computes the RHS of the ODE system
         * GetIIonic returns the total ionic current
        """
        self.output_get_i_ionic()
        self.output_evaluate_y_derivatives(method_name='EvaluateYDerivatives')
        if self.use_analytic_jacobian:
            self.output_jacobian()
        self.output_derived_quantities()
    
    def output_bottom_boilerplate(self):
        """Call superclass method, then end the CHASTE_CVODE guard."""
        super(CellMLToCvodeTranslator, self).output_bottom_boilerplate()
        # CVODE is optional in Chaste
        self.writeln("#endif // CHASTE_CVODE")
        self.writeln_hpp("#endif // CHASTE_CVODE")
    
    def output_extra_constructor_content(self):
        """Tell the base class if we have an analytic Jacobian."""
        if self.use_analytic_jacobian:
            self.writeln('mUseAnalyticJacobian = true;')
            self.writeln('mHasAnalyticJacobian = true;')
    
    def _count_operators(self, exprs, result=None):
        if result is None: result = {}
        for expr in exprs:
            if isinstance(expr, mathml_apply):
                op = expr.operator().localName
                result[op] = 1 + result.get(op, 0)
            children = expr.xml_element_children()
            if children:
                self._count_operators(children, result)
        return result
        
    def output_jacobian(self):
        """Output an analytic Jacobian for CVODE to use."""
        self.output_method_start('EvaluateAnalyticJacobian',
                                 [self.TYPE_DOUBLE + self.code_name(self.free_vars[0]),
                                  self.TYPE_VECTOR + 'rY', self.TYPE_VECTOR + 'rDY',
                                  'CHASTE_CVODE_DENSE_MATRIX rJacobian',
                                  self.TYPE_VECTOR + 'rTmp1', self.TYPE_VECTOR + 'rTmp2', self.TYPE_VECTOR + 'rTmp3'],
                                 'void', access='public')
        self.open_block()
        # Mathematics that the Jacobian depends on
        used_vars = set([self.config.dt_variable])
        for entry in self.model.solver_info.jacobian.entry:
            used_vars.update(self._vars_in(entry.math))
        nodeset = self.calculate_extended_dependencies(used_vars, prune_deps=[self.doc._cml_config.i_stim_var])
        self.output_state_assignments(nodeset=nodeset, assign_rY=False)
        table_index_nodes_used = self.calculate_lookup_table_indices(nodeset|set(map(lambda e: e.math, self.model.solver_info.jacobian.entry)), self.code_name(self.free_vars[0]))
        self.output_equations(nodeset - table_index_nodes_used)
        self.writeln()
        # Jacobian entries, sorted by index with rows varying fastest
        self.output_comment('Matrix entries')
        entries = []
        def gv(vn):
            return self.varobj(vn).get_source_variable(recurse=True)
        for entry in self.model.solver_info.jacobian.entry:
            var_i, var_j = gv(entry.var_i), gv(entry.var_j)
            i = self.state_vars.index(var_i)
            j = self.state_vars.index(var_j)
            entry_content = list(entry.math.xml_element_children())
            assert len(entry_content) == 1, "Malformed Jacobian entry: " + entry.xml()
            entry = entry_content[0]
            if not (isinstance(entry, mathml_cn) and entry.evaluate() == 0.0):
                entries.append((j, i, var_i is self.v_variable, entry))
        entries.sort()
        for j, i, is_V, entry in entries:
            self.writeln('IJth(rJacobian, ', i, ', ', j, ') = ', self.code_name(self.config.dt_variable), ' * (', nl=False)
            paren = False
            if is_V:
                self.write('mSetVoltageDerivativeToZero ? 0.0 : ')
                paren = True
            self.output_expr(entry, paren)
            self.writeln(')', self.STMT_END, indent=False)
#        self.output_comment('Debugging!')
#        self.writeln('#ifndef NDEBUG', indent=False)
#        self.writeln('for (long int j=0; j<N; j++)')
#        self.open_block()
#        self.writeln('for (long int i=0; i<N; i++)')
#        self.open_block()
#        self.writeln('if (std::isnan(IJth(rJacobian, i, j)))')
#        self.open_block()
#        self.writeln('std::cerr << "NAN at J(" << i << "," << j << ")" << DumpState("", rY);')
#        self.close_block(blank_line=False)
#        self.writeln('EXCEPT_IF_NOT(!std::isnan(IJth(rJacobian, i, j)));')
#        self.close_block(blank_line=False)
#        self.close_block(blank_line=False)
#        self.writeln('//CheckAnalyticJacobian(', self.code_name(self.free_vars[0]),
#                     ', rY, rDY, rJacobian, rTmp1, rTmp2, rTmp3);')
#        self.writeln('#endif // NDEBUG', indent=False)
        self.close_block()

class CellMLToMapleTranslator(CellMLTranslator):
    """Translate a CellML model to Maple code."""

    # Language tokens that differ from the default
    EQ_ASSIGN = ' := '   # Assignment operator
    COMMENT_START = '# ' # Start of a 1 line comment
    # Types are determined automatically by Maple
    TYPE_DOUBLE = ''
    TYPE_CONST_DOUBLE = ''
    # Some constants are different
    PI = 'Pi'
    E = 'exp(1)'

    def __init__(self, omit_constants=False, compute_full_jacobian=False,
                 **kwargs):
        """Create a Maple translator.

        If omit_constants is set to true, assignments will not be
        generated for constant variables.  This should be used if
        these values will be altered at runtime, in order to prevent
        derivatives being calculated incorrectly.

        Set compute_full_jacobian to True to make Maple compute the
        Jacobian of the whole ODE system, rather than just the
        nonlinear portion.

        Other keyword arguments are all passed to the base class.
        """
        super(CellMLToMapleTranslator, self).__init__(**kwargs)
        # Maple translation doesn't support lookup tables
        self.use_lookup_tables = False
        # Translation parameters
        self.omit_constants = omit_constants
        self.compute_full_jacobian = compute_full_jacobian
        # Update some function names
        self.function_map = CellMLTranslator.function_map.copy()
        del self.function_map['power']
        self.function_map.update(
            {'abs': 'abs', 'ln': 'ln', 'not': 'not',
             'sec': 'sec', 'csc': 'csc', 'cot': 'cot',
             'sech': 'sech', 'csch': 'csch', 'coth': 'coth',
             'arcsin': 'arcsin', 'arccos': 'arccos', 'arctan': 'arctan',
             'arcsec': 'arcsec', 'arccsc': 'arccsc', 'arccot': 'arccot',
             'arcsinh': 'arcsinh', 'arccosh': 'arccosh', 'arctanh': 'arctanh',
             'arcsech': 'arcsech', 'arccsch': 'arccsch', 'arccoth': 'arccoth'})
        self.recip_trig = {}
        self.nary_ops = CellMLTranslator.nary_ops.copy()
        self.nary_ops.update({'and': 'and', 'or': 'or'})
        self.binary_ops = CellMLTranslator.binary_ops.copy()
        self.binary_ops.update({'xor': 'xor', 'eq': '=', 'neq': '<>',
                                'power': '^'})
        self.special_roots = {}

    def output_file_name(self, model_filename):
        """Generate a name for our output file, based on the input file."""
        return os.path.splitext(model_filename)[0] + '.mpl'

    def output_top_boilerplate(self):
        """Output top boilerplate."""
        self.writeln('# Model: ', self.model.name)
        self.output_comment(version_comment(self.add_timestamp))
        self.writeln()
        self.writeln('interface(prettyprint=0);\n')
        if self.compute_full_jacobian:
            self.writeln('print("FULL JACOBIAN");')
        return
#
#    def output_mathematics(self):
#        if self.compute_full_jacobian:
#            pass
#        else:
#            super(CellMLToMapleTranslator, self).output_mathematics()

    def output_bottom_boilerplate(self):
        """Output bottom boilerplate."""
        self.output_comment('\nJacobian calculation\n')
        if self.compute_full_jacobian:
            state_vars = self.state_vars
            # Record the ordering of the state variables, since they're refered to by index below
            for i, var in enumerate(state_vars):
                self.writeln('print("--%d--%s--");' % (i+1, self.code_name(var)))
            # Jacobian calculation for the whole ODE system, i.e. each df_i/du_j
            state_var_names = map(self.code_name, state_vars)
            self.writeln('jacobian', self.EQ_ASSIGN, 'array(')
            self.write('[')
            for var_i in state_vars:
                if var_i is not state_vars[0]: self.write(',')
                self.write('[')
                for var_j in state_vars:
                    if var_j is not state_vars[0]: self.write(',')
                    self.write('diff(', self.code_name(var_i, ode=True), ',', self.code_name(var_j), ')')
                self.write(']')
            self.write(']);\n')
            self.writeln('with(codegen):')
            self.writeln('J', self.EQ_ASSIGN, 'optimize(jacobian);')
        elif hasattr(self.model, '_cml_nonlinear_system_variables'):
            # Jacobian calculation for Jon Whiteley's algorithm
            vars_text = self.model._cml_nonlinear_system_variables
            if type(vars_text) == type([]):
                var_objs = self.model._cml_nonlinear_system_variables
            else:
                # Get variable objects from names
                varnames = map(lambda s: s.split(','), vars_text.split(';'))
                var_objs = map(lambda (c, v):
                               self.model.get_variable_by_name(c, v),
                               varnames)
            # Output the Newton iteration expression for each variable
            for var in var_objs:
                self.writeln('g_', self.code_name(var), self.EQ_ASSIGN,
                             self.code_name(var), ' - ',
                             self.code_name(var), '_old - delta_t*',
                             self.code_name(var, ode=True), self.STMT_END)
            # Output the derivative calculations
            for var_i in var_objs:
                for var_j in var_objs:
                    self.writeln('print("--', self.code_name(var_i), '/',
                                 self.code_name(var_j), '--");')
                    self.writeln('diff(g_', self.code_name(var_i), ', ',
                                 self.code_name(var_j), ');')
        # Tell Maple to quit when done
        self.writeln()
        self.writeln('quit;')
        return

    def output_assignment(self, expr):
        """Output an assignment expression.

        Optionally, if this is an assignment of a constant, don't output,
        so that differentation doesn't optimise expressions away.
        """
        if isinstance(expr, cellml_variable):
            # This may be the assignment of a mapped variable, or a constant
            t = expr.get_type()
            if t == VarTypes.Mapped:
                self.writeln(self.TYPE_DOUBLE, self.code_name(expr),
                             self.EQ_ASSIGN,
                             self.code_name(expr.get_source_variable()),
                             self.STMT_END)
            elif t == VarTypes.Constant and not self.omit_constants:
                self.writeln(self.TYPE_CONST_DOUBLE, self.code_name(expr),
                             self.EQ_ASSIGN, nl=False)
                self.output_number(expr.initial_value)
                self.writeln(self.STMT_END, indent=False)
        else:
            # This is a mathematical expression
            self.writeln(self.TYPE_DOUBLE, nl=False)
            opers = expr.operands()
            self.output_lhs(opers.next())
            self.write(self.EQ_ASSIGN)
            self.output_expr(opers.next(), False)
            self.writeln(self.STMT_END, indent=False)
        return

    def output_number(self, expr):
        """Output the plain number expr.
        
        With Maple, there is no need to make all constants parse as
        doubles to avoid problems with integer division or numbers too
        large for the int type.
        
        Negative numbers will be prefixed by a space to avoid unwanted
        decrement operations.
        """
        n = self.eval_number(expr)
        num = "%.17g" % n
        if num[0] == '-':
            num = '(' + num + ')'
        self.write(num)

    def output_root(self, expr, paren):
        """Output a root taken to some degree.

        If a degree qualifier element is not provided, uses default 2.
        """
        if hasattr(expr, u'degree'):
            # A degree is given.  Compute x^(1/b)
            self.open_paren(paren)
            self.output_expr(expr.operands().next(), True)
            self.write('^(1/')
            self.output_expr(expr.degree, True)
            self.write(')')
            self.close_paren(paren)
        else:
            # Compute square root
            self.output_function('sqrt', expr.operands(), paren)

    def output_log(self, expr, paren):
        """Output a logarithm to the given base, which defaults to base 10."""
        if hasattr(expr, u'logbase'):
            # A base is provided.  Use the log[b](x) function.
            self.write('log[')
            self.output_expr(expr.logbase, False)
            self.write(']')
            self.output_function('', expr.operands(), paren)
        else:
            # Use base 10
            self.output_function('log10', expr.operands(), paren)

    def output_piecewise(self, expr, paren):
        """Output the piecewise expression expr.

        We use the if operator.
        """
        self.write('piecewise(')
        need_comma = False
        for piece in getattr(expr, u'piece', []):
            if need_comma:
                self.write(',')
            self.output_expr(child_i(piece, 2), False) # Condition
            self.write(',')
            self.output_expr(child_i(piece, 1), False) # Result
            need_comma = True
        if hasattr(expr, u'otherwise'):
            if need_comma:
                self.write(',')
            self.output_expr(child_i(expr.otherwise, 1), paren) # Default case
        self.write(')')

class CellMLToHaskellTranslator(CellMLTranslator):
    """Translate a CellML model to a Haskell version.

    This does not produce a 'runnable' version of the model, but
    rather an equivalent model in effectively an abstract syntax,
    which can then be interpreted by a suitable interpreter.
    
    This allows us to more easily specify an operational semantics for
    CellML, without having to worry about the XML issues in the
    interpreter itself.
    """

    STMT_END = ''
    COMMENT_START = '-- '
    TYPE_DOUBLE = ''
    TYPE_CONST_DOUBLE = ''
    E = '(exp 1)'
    PI = 'pi'
    TRUE = '(Bool True)'
    FALSE = '(Bool False)'

    def __init__(self, **kwargs):
        """Create a Haskell translator.

        Keyword arguments are all passed to the base class.
        """
        super(CellMLToHaskellTranslator, self).__init__(**kwargs)
        # We don't use lookup tables in Haskell code
        self.use_lookup_tables = False
        return

    def output_file_name(self, model_filename):
        """Generate a name for our output file, based on the input file."""
        return os.path.splitext(model_filename)[0] + '.hs'

    def stringify(self, s):
        """Quote a string."""
        return '"' + s + '"'

    def code_name(self, var, ode=False, full=False):
        """
        Return the full name of var in a form suitable for inclusion in a
        source file.
        
        The functionality of ode=True is implemented in output_apply
        rather than here, so this parameter must be False.

        If full is True, include the name of the owning component.

        If PE has been performed (there is only 1 component, and variable
        names have been munged) then transform the munged name to Haskell
        munged form.
        """
        if ode:
            raise NotImplementedError # Never used; see output_apply.
        if self.single_component and '__' in var.name:
            name = var.name.replace('__', ',')
        elif full:
            name = var.xml_parent.name + ',' + var.name
        else:
            name = var.name
        return self.stringify(name)

    def output_top_boilerplate(self):
        """Output top boilerplate.

        Outputs the imports and model-level units definitions.
        Components and connections are output by output_mathematics.
        """
        self.class_name = self.class_name.lower()
        self.module_name = self.class_name[0].upper() + self.class_name[1:]
        self.writeln('module ', self.module_name, ' where')
        self.writeln('-- Model: ', self.model.name)
        self.output_comment(version_comment(self.add_timestamp))
        self.writeln()
        self.writeln('import CellML')
        self.writeln('import Units')
        self.writeln('import Environment')
        self.writeln()
        # Model definition
        self.writeln(self.class_name, ' = Model ',
                     self.stringify(self.class_name),
                     ' units components connections')
        self.writeln('  where')
        self.set_indent(offset=1)
        # Model-level units definitions
        self.writeln('units =')
        self.output_units(self.model)
        return

    def output_unit(self, udef):
        """Output a single units definition, recursively."""
        def prefix(uref):
            """Return a prefix of a units reference,
            as an integer to which 10 may be raised."""
            prefix = getattr(uref, u'prefix_', 0)
            if prefix in uref.SI_PREFIXES:
                prefix = uref.SI_PREFIXES[prefix]
            else:
                prefix = int(prefix)
            return prefix
        def num(n):
            """Wrap numbers for Haskell."""
            if n < 0:
                return "(" + str(n) + ")"
            else:
                return n

        self.write('(')
        if udef.is_base_unit():
            # Base unit
            self.write('BaseUnits ', self.stringify(udef.name))
        elif udef.is_simple():
            # Simple units
            self.write('SimpleUnits ', num(udef.get_multiplier()), ' ',
                       num(prefix(udef.unit)), ' ')
            self.output_unit(udef.unit.get_units_element())
            self.write(' ', num(udef.get_offset()))
        else:
            # Complex units
            self.write('ComplexUnits [')
            uref_comma = False
            for uref in udef.unit:
                if uref_comma: self.write(',')
                else: uref_comma = True
                self.write('Unit ', num(uref.get_multiplier()), ' ',
                           num(prefix(uref)), ' ')
                self.output_unit(uref.get_units_element())
                self.write(' ', num(uref.get_exponent()))
            self.write(']')
        self.write(')')
        return

    def output_units(self, units_parent):
        """Output the units definitions in this model or component."""
        self.open_list()
        comma, output = False, False
        for udef in getattr(units_parent, u'units', []):
            output = True
            if comma: self.writeln(',', indent=False)
            else: comma = True
            # Output a single definition
            self.writeln('UDef ', self.stringify(udef.name), nl=False)
            self.output_unit(udef)
        if output:
            self.writeln('', indent=False)
        self.close_list()
        return

    def output_mathematics(self):
        """Output the mathematics in this model.

        This method outputs the components and connections."""
        # Components
        self.writeln('components =')
        self.open_list()
        comma = False
        for comp in getattr(self.model, u'component', []):
            if comma: self.writeln(',')
            else: comma = True
            self.output_component(comp)
        self.writeln('')
        self.close_list()
        # Connections
        self.writeln('connections =')
        self.open_list()
        comma = False
        for var in (v for v in self.model.get_assignments()
                    if isinstance(v, cellml_variable)
                    if v.get_type() == VarTypes.Mapped):
            if comma: self.writeln(',', indent=False)
            else: comma = True
            self.output_connection(var)
        self.writeln('', indent=False)
        self.close_list()
        return

    def output_component(self, comp):
        """Output a single component."""
        self.writeln('MkComp ', self.stringify(comp.name))
        self.output_units(comp)
        # Variable declarations, associating units with var names
        self.open_list()
        comma = False
        for var in getattr(comp, u'variable', []):
            if comma: self.writeln(',')
            else: comma = True
            self.writeln('VarDecl ', self.code_name(var), ' ',
                         self.stringify(var.units))
        self.close_list()
        # And now the mathematics
        self.open_list()
        comma = False
        # Constants
        for var in (v for v in getattr(comp, u'variable', [])
                    if v.get_type() == VarTypes.Constant):
            if comma: self.writeln(',')
            else: comma = True
            self.output_assignment(var)
        # Expressions
        for math in getattr(comp, u'math', []):
            for expr in getattr(math, u'apply', []):
                if comma: self.writeln(',')
                else: comma = True
                self.output_assignment(expr)
        self.close_list()
        return

    def output_connection(self, conn):
        """Output a single connection."""
        to_var = conn
        from_var = conn.get_source_variable()
        self.writeln('VarMap', nl=False)
        self.write(' (', self.stringify(from_var.xml_parent.name), ',',
                   self.stringify(from_var.name), ')')
        self.write(' (', self.stringify(to_var.xml_parent.name), ',',
                   self.stringify(to_var.name), ')')
        return

    def output_bottom_boilerplate(self):
        """Output bottom boilerplate."""
        self.set_indent(offset=-1)
        self.writeln()
        self.output_comment('Evaluate derivatives at the start of time.')
        self.writeln()
        # Initial environment
        self.writeln('initial_environment :: Env')
        self.writeln('initial_environment = make_env')
        self.open_list()
        self.writeln('  (Var ', self.code_name(self.free_vars[0], full=True),
                     ', Val (Number 0))')
        for sv in self.state_vars:
            self.writeln(', (Var ', self.code_name(sv, full=True),
                         ', Val ', nl=False)
            self.output_number(sv.initial_value, as_value=True)
            self.writeln(')', indent=False)
        self.close_list()
        self.writeln('results = run_cellml ', self.class_name,
                     ' initial_environment')
        self.writeln()
        # Dynamic environment for PE
        self.writeln('dynamic_environment :: Env')
        self.writeln('dynamic_environment = foldr def initial_environment')
        self.open_list()
        # Include all variables marked as pe:keep
        comma = False
        for comp in getattr(self.model, u'component', []):
            for var in getattr(comp, u'variable', []):
                if var.pe_keep:
                    self.writeln([' ', ','][comma], ' (Var ',
                                 self.code_name(var, full=True),
                                 ', Val DynamicMarker)')
                    if not comma: comma = True
        self.close_list()
        self.writeln('where def (k,v) env = define env k v', indent_offset=1)
        self.writeln('pe_results = reduce_and_run_cellml ', self.class_name,
                     ' dynamic_environment')
        return

    def open_list(self):
        """Open a multi-line list."""
        self.set_indent(offset=1)
        self.writeln('[')
        self.set_indent(offset=1)
        return

    def close_list(self):
        """Close a multi-line list."""
        self.set_indent(offset=-1)
        self.writeln(']')
        self.set_indent(offset=-1)
        return

    def output_assignment(self, expr):
        """Output an assignment expression."""
        if isinstance(expr, cellml_variable):
            # Assignment of a constant
            self.writeln('Assign (Var ', self.code_name(expr), ') ', nl=False)
            self.output_number(expr.initial_value, units=expr.units)
            self.writeln(indent=False)
        else:
            # This is a mathematical expression
            opers = expr.operands()
            self.writeln('Assign ', nl=False)
            self.output_lhs(opers.next())
            self.write(' ')
            self.output_expr(opers.next(), True)
            self.writeln(indent=False)
        return

    def output_lhs(self, expr):
        """Output the left hand side of an assignment expression."""
        if expr.localName == 'ci':
            self.output_variable(expr, lhs=True)
        elif expr.operator().localName == 'diff':
            v1 = expr.operator().dependent_variable
            v2 = expr.operator().independent_variable
            self.write('(Ode ', self.code_name(v1), ' ', self.code_name(v2),
                       ')')
        return

    def output_variable(self, ci_elt, lhs=False):
        """Output a ci element, i.e. a variable lookup."""
        type_constructor = ['Variable', 'Var'][lhs]
        self.write('(', type_constructor, ' ',
                   self.code_name(ci_elt.variable), ')')
        return
    
    def output_number(self, expr, as_value=False, units=None):
        """Output the plain number expr.
        
        With Haskell there is no need to force numbers to parse as doubles.
        We do need to bracket negative numbers.
        """
        n = self.eval_number(expr)
        num = "%.17g" % n
        if num[0] == '-':
            num = "(" + num + ")"
        tc = ['Num', 'Number'][as_value]
        if not as_value:
            if units is None:
                units = getattr(expr, u'units', '')
            units = '(Left ' + self.stringify(units) + ')'
        else:
            units = ''
        self.write("(", tc, " ", num, " ", units, ")")
        return
    
    def output_apply(self, expr, paren):
        """Output an <apply> expression.
        
        paren is True if the context has requested parentheses.
        """
        op = expr.operator()
        op_name = op.localName.capitalize()
        self.open_paren(paren)
        # Some operators are special-cased, but most map directly
        if op_name == u'Root':
            self.output_special_apply(expr, u'Root', u'degree', u'Sqrt')
        elif op_name == u'Log':
            self.output_special_apply(expr, u'Log', u'logbase', u'Ln')
        elif op_name == u'Diff':
            if self.single_component:
                # A bit of a hack, due to the odd way this case is
                # handled by other translators - in effect the
                # translator has to do some PE...
                self.write('Apply Diff [Variable ',
                           self.code_name(
                    op.dependent_variable.get_source_variable(recurse=True)),
                           ', Variable ',
                           self.code_name(
                    op.independent_variable.get_source_variable(recurse=True)),
                           ']')
            else:
                self.really_output_apply(op_name, list(expr.operands()) +
                                         [expr.bvar.ci])
        else:
            self.really_output_apply(op_name, expr.operands())
        self.close_paren(paren)
        return

    def really_output_apply(self, operator, operands):
        """Actually output code for the application of an operator."""
        self.write('Apply ', operator, ' [')
        comma = False
        for operand in operands:
            if comma: self.write(',')
            else: comma = True
            self.output_expr(operand, False)
        self.write(']')
        return

    def output_special_apply(self, expr, op_name, qual_name, special_name):
        """Output a special-cased apply expression.

        op_name is the name of the general case operator.  If the
        expression has a qualifier called qual_name, this will be
        used, with the qualifier's value as second operand.  Otherwise,
        the operator special_name will be used, with a single operand.
        """
        if hasattr(expr, qual_name):
            self.really_output_apply(op_name, [expr.operands().next(),
                                               getattr(self, qual_name)])
        else:
            self.really_output_apply(special_name, expr.operands())
        return

    def output_piecewise(self, expr, paren):
        """Output the piecewise expression expr."""
        self.open_paren(paren)
        self.write('Piecewise [')
        comma = False
        for piece in getattr(expr, u'piece', []):
            if comma: self.write(',')
            else: comma = True
            self.write('Case ')
            self.output_expr(child_i(piece, 2), True) # Condition
            self.write(' ')
            self.output_expr(child_i(piece, 1), True) # Result
        self.write('] ')
        if hasattr(expr, u'otherwise'):
            self.write('(Just ')
            self.output_expr(child_i(expr.otherwise, 1), True) # Default case
            self.write(')')
        else:
            self.write('Nothing')
        self.close_paren(paren)
        return

class CellMLToMatlabTranslator(CellMLTranslator):
    """Translate a CellML model to Matlab code.

    The normal case generates a .m file such as could be used with ode45.
    When lookup tables are used (TODO), the file generated represents a
    function that returns a function handle, suitable for use with ODE
    solvers.
    """

    # Language tokens that differ from the default
    COMMENT_START = '% ' # Start of a 1 line comment
    # Types are determined automatically by Matlab
    TYPE_DOUBLE = ''
    TYPE_CONST_DOUBLE = ''
    # Some constants are different
    PI = 'pi'
    E = 'exp(1)'
    NOT_A_NUMBER = 'NaN'

    def __init__(self, **kwargs):
        super(CellMLToMatlabTranslator, self).__init__(**kwargs)
        # Update some function, etc. names
        self.function_map = CellMLTranslator.function_map.copy()
        self.function_map.update(
            {'power': 'power', 'abs': 'abs',
             'xor': 'xor', 'not': '~', 'rem': 'rem',
             'sec': 'sec', 'csc': 'csc', 'cot': 'cot',
             'sech': 'sech', 'csch': 'csch', 'coth': 'coth',
             'arcsec': 'asec', 'arccsc': 'acsc', 'arccot': 'acot',
             'arcsech': 'asech', 'arccsch': 'acsch', 'arccoth': 'acoth'
             })
        self.recip_trig = {}
        self.binary_ops = CellMLTranslator.binary_ops.copy()
        del self.binary_ops['xor']
        self.binary_ops['neq'] = '~='
        self.binary_ops['divide'] = './'
        self.nary_ops = CellMLTranslator.nary_ops.copy()
        self.nary_ops['times'] = '.*'
        self.special_roots = {2: 'sqrt'}

    def output_file_name(self, model_filename):
        """Generate a name for our output file, based on the input file."""
        name = os.path.splitext(model_filename)[0] + '.m'
        # Matlab doesn't like long names :(
        if len(name) > 60:
            # Take end part so we get version/variant info if present
            name = name[-60:]
        return name

    def translate(self, doc, *args, **kwargs):
        """Generate code for the model or its Jacobian matrix."""
        self.variable_name_map = {}
        if hasattr(doc.model, u'solver_info') and \
               hasattr(doc.model.solver_info, u'jacobian'):
            kwargs['continuation'] = self.output_jacobian
        if 'output_filename' in kwargs and len(kwargs['output_filename'])>60:
            # Take end part so we get version/variant info if present
            kwargs['output_filename'] = kwargs['output_filename'][-60:]
        return super(CellMLToMatlabTranslator,
                     self).translate(doc, *args, **kwargs)

    def output_top_boilerplate(self):
        """Output top boilerplate."""
        self.output_comment(version_comment(self.add_timestamp))
        t = self.code_name(self.free_vars[0])
        # Matlab doesn't like long names :(
        if len(self.class_name) > 60:
            # Take end part so we get version/variant info if present
            self.class_name = self.class_name[-60:]
            # Strip leading underscores
            while self.class_name[0] == '_':
                self.class_name = self.class_name[1:]
        
        if self.use_lookup_tables:
            self.writeln('function dy_fun_ptr = ', self.class_name, '_lt(step)')
            self.output_comment('Generate a function to evaluate using '
                                'lookup tables the model ', self.model.name, '.')
            self.output_comment('The function returned is f, where dU/dt = f(t, U).')
            self.writeln()
            self.set_indent(offset=1)
            self.output_lut_generation()
            self.output_lut_lookups()
            self.writeln('tables = generate_tables(step);')
        else:
            self.writeln('function [dy_fun_ptr initial_values V_index t_units state_var_names] = ',
                         self.class_name, '()')
            self.output_comment('Get evaluation function and metadata for the model ',
                                self.model.name, '.')
            self.output_comment('\nReturns the function f (where dU/dt = f(t, U)),\n'
                                'suitable initial values for the system,\n'
                                'the index of the transmembrane potential within '
                                'the state variable vector,\n'
                                'the multiplicative factor of the time units,\n'
                                'and the names of the state variables.')
            self.set_indent(offset=1)
            self.writeln('V_index = ', self.v_index+1, ';')
            self.writeln('state_var_names = cell(1, ', len(self.state_vars), ');')
            self.writeln('initial_values = zeros(1, ', len(self.state_vars), ');')
            for i, var in enumerate(self.state_vars):
                self.writeln('state_var_names{', i+1, '}', self.EQ_ASSIGN,
                             "'", var.fullname(), "';")
                self.writeln('initial_values(', i+1, ')', self.EQ_ASSIGN,
                             getattr(var, u'initial_value', self.NOT_A_NUMBER), ';')
            t_var = self.free_vars[0]
            t_units = t_var.component.get_units_by_name(t_var.units)
            self.writeln('t_units = ', t_units.get_multiplicative_factor(), ';')
        self.writeln('function dy = dy_fun(',t,', y)')
        self.set_indent(offset=1)
        self.output_comment('Time units: ', self.free_vars[0].units)
        self.writeln()
        for i, var in enumerate(self.state_vars):
            self.writeln(self.code_name(var), self.EQ_ASSIGN, 'y(', i+1,
                         ');')
            self.output_comment('Units: ', var.units, '; Initial value: ',
                                getattr(var, u'initial_value', 'Unknown'))
        self.writeln()
        if self.use_lookup_tables:
            for key, i in self.doc.lookup_table_indexes.iteritems():
                i = int(i) + 1
                min, max, step, var = key
                varname = self.code_name(var)
                self.writeln('table_values{', i, '} = lookup_', i,
                             '(tables, ', varname, ', step);')
        self.writeln()
        return

    def output_bottom_boilerplate(self):
        """Output bottom boilerplate."""
        self.writeln()
        self.writeln('dy = zeros(size(y));')
        for i, var in enumerate(self.state_vars):
            self.writeln('dy(', str(i+1), ') = ',
                         self.code_name(var, ode=True), ';')
        self.set_indent(offset=-1)
        self.writeln('end')
        self.writeln()
        self.writeln('dy_fun_ptr = @dy_fun;')
        self.set_indent(offset=-1)
        self.writeln('end')

    def code_name(self, var, ode=False, shorten=True):
        """Matlab has an upper limit on the length of variable names!"""
        full_name = super(CellMLToMatlabTranslator, self).code_name(var, ode)
        if shorten:
            full_name = self.shorten_name(full_name)
        return full_name

    def shorten_name(self, var_name):
        """If the name is too long for Matlab, shorten it."""
        if len(var_name) > 60:
            # Actual bound is 63, but let's be cautious
            try:
                return self.variable_name_map[var_name]
            except KeyError:
                new_name = 'shortened_var_' + str(len(self.variable_name_map))
                self.variable_name_map[var_name] = new_name
                return new_name
        else:
            return var_name

    def output_number(self, expr):
        """Output the plain number expr.
        
        With Matlab, there is no need to make all constants parse as
        doubles to avoid problems with integer division or numbers too
        large for the int type.
        
        Negative numbers will be prefixed by a space to avoid unwanted
        decrement operations.
        """
        n = self.eval_number(expr)
        num = "%.17g" % n
        if num[0] == '-':
            num = ' ' + num
        self.write(num)

    def output_root(self, expr, paren):
        """Output a root taken to some degree.

        If a degree qualifier element is not provided, uses default 2.
        """
        if hasattr(expr, u'degree'):
            # A degree is given.  Compute nthroot(x, b)
            x = expr.operands().next()
            b = expr.degree
            self.output_function('nthroot', [x, b], paren)
        else:
            # Compute square root
            self.output_function('sqrt', expr.operands(), paren)

    def output_piecewise(self, expr, paren):
        """Output the piecewise expression expr.

        Uses an ifexpr.m file to code if expressions.
        """
        num_ifs = 0
        for piece in getattr(expr, u'piece', []):
            num_ifs += 1
            self.write('ifexpr(')
            self.output_expr(child_i(piece, 2), False) # Condition
            self.write(',')
            self.output_expr(child_i(piece, 1), False) # Result
            self.write(',')
        if hasattr(expr, u'otherwise'):
            self.output_expr(child_i(expr.otherwise, 1), paren) # Default case
        else:
            self.write(self.NOT_A_NUMBER) # If this is hit, things get ugly
        for i in range(num_ifs):
            self.close_paren(True)

    def output_lut_generation(self):
        """Output code to generate lookup tables.

        There should be a list of suitable expressions available as
        self.doc.lookup_tables, to save having to search the whole
        model.
        """
        self.writeln('function tables = generate_tables(step)')
        self.set_indent(offset=1)
        self.output_comment('Generate all the lookup tables for this model.\n'
                            'Returns a cell array containing matrices, each column of '
                            'which contain one table.')
        self.use_lookup_tables = False
        for key, idx in self.doc.lookup_table_indexes.iteritems():
            min, max, step, var = key
            i = int(idx) + 1
            table_extent = unicode(float(max) - float(min))
            num_tables = unicode(self.doc.lookup_tables_num_per_index[idx])
            self.writeln('tables{', i, '} = zeros(1+floor(', table_extent, '/step),',
                         num_tables, ');')
        for expr in self.doc.lookup_tables:
            j = int(expr.table_name) + 1
            i = int(expr.table_index) + 1
            var = expr.get_component().get_variable_by_name(expr.var)
            varname = self.code_name(var)
            self.writeln(varname, ' = [', expr.min, ':step:', expr.max, '];')
            self.writeln('tables{', i, '}(:,', j, ') = ', nl=False)
            self.output_expr(expr, False)
            self.writeln(';', indent=False)
        self.use_lookup_tables = True
        self.set_indent(offset=-1)
        self.writeln('end')
        self.writeln()

    def output_lut_lookups(self):
        """Output the functions that perform table lookups."""
        for key, idx in self.doc.lookup_table_indexes.iteritems():
            i = int(idx) + 1
            min, max, step, var = key
            self.writeln('function val = lookup_', i, '(tables, var, step)')
            self.set_indent(offset=1)
            self.output_comment('Lookup all tables for variable var')
            self.writeln('if ~isreal(var)')
            self.writeln("error(['Index variable value ' num2str(var) ' is not real'])",
                         indent_offset=1)
            self.writeln('end')
            self.writeln('table_lower = ', min, ';')
            self.writeln('table_upper = ', max, ';')
            self.writeln('if var < table_lower || var >= table_upper')
            self.writeln("error(['Index variable value ' num2str(var) ' outside table bounds'])",
                         indent_offset=1)
            self.writeln('end')
            self.writeln('i = 1 + floor((var - table_lower)/step);')
            self.writeln('y1 = tables{', i, '}(i, :);')
            self.writeln('y2 = tables{', i, '}(i+1, :);')
            self.writeln('var_i = table_lower + step*(i-1);')
            self.writeln('val = y1 + (y2-y1) .* (var-var_i) ./ step;')
            self.set_indent(offset=-1)
            self.writeln('end')
            self.writeln()

    def output_table_lookup(self, expr, paren):
        """Output code to look up expr in the appropriate table."""
        i = int(expr.table_index) + 1
        j = int(expr.table_name) + 1
        self.write('table_values{', i, '}(', j, ')')

    def output_jacobian(self):
        """Generate code to compute the Jacobian matrix for this model."""
        t = self.code_name(self.free_vars[0])
        self.writeln('function J = jacobian(',t,', y)')
        self.set_indent(offset=1)
        self.output_comment('Jacobian matrix for the model ', self.model.name)
        self.output_comment('Evaluates the matrix J, where J(j,i) = d f_i / d u_j')
        # State variable assignments
        for i, var in enumerate(self.state_vars):
            self.writeln(self.code_name(var), self.EQ_ASSIGN, 'y(', str(i+1),
                         ');')
            self.output_comment('Units: ', var.units, '; Initial value: ',
                                getattr(var, u'initial_value', 'Unknown'))
        # Mathematics that the Jacobian depends on
        used_vars = set()
        for entry in self.model.solver_info.jacobian.entry:
            used_vars.update(self._vars_in(entry.math))
        nodeset = self.calculate_extended_dependencies(used_vars)
        self.output_equations(nodeset)
        self.writeln()
        # Jacobian entries
        state_var_names = map(lambda v: self.code_name(v, shorten=False),
                              self.state_vars)
        self.writeln('J = zeros(length(y));')
        for entry in self.model.solver_info.jacobian.entry:
            var_i, var_j = entry.var_i, entry.var_j
            i = state_var_names.index(var_i) + 1
            j = state_var_names.index(var_j) + 1
            self.writeln('J(', j, ',', i, ') = ', nl=False)
            entry_content = list(entry.math.xml_element_children())
            assert len(entry_content) == 1
            self.output_expr(entry_content[0], False)
            self.writeln(self.STMT_END, indent=False)
        self.set_indent(offset=-1)
        self.writeln('end')

    def output_variable(self, ci_elt, ode=False):
        """Output a ci element, i.e. a variable lookup."""
        if hasattr(ci_elt, '_cml_variable') and ci_elt._cml_variable:
            self.write(self.code_name(ci_elt.variable, ode=ode))
        else:
            # This ci element is in the solver_info section, thus
            # doesn't have all the extra annotations.  It is a fully
            # qualified name though.
            prefix = ['var_', 'd_dt_'][ode]
            varname = unicode(ci_elt)
            if varname[0] == '(':
                # (compname,varname)
                cname, vname = varname[1:-1].split(u',')
                if self.single_component:
                    varname = vname
                else:
                    varname = cname + '__' + vname
            elif varname == u'delta_t':
                # Special case for the timestep in ComputeJacobian
                prefix = ''
                varname = 'mDt'
            else:
                # var_cname__vname
                varname = varname[4:]
            self.write(self.shorten_name(prefix + varname))
        return



class CellMLToPythonTranslator(CellMLToChasteTranslator):
    """Output code suitable for the Python implementation of Functional Curation."""
    
    STMT_END = ''
    COMMENT_START = '# '
    DOXYGEN_COMMENT_START = '## '
    LOGICAL_AND = ' and '
    LOGICAL_OR = ' or '
    LOGICAL_TRUE = ' True '
    ASSERT = 'assert '
    TYPE_DOUBLE = ''
    TYPE_CONST_DOUBLE = ''
    TYPE_VOID = ''
    TYPE_CONST_UNSIGNED = ''
    TYPE_VECTOR = ''
    TYPE_VECTOR_REF = ''
    TRUE = 'True'
    FALSE = 'False'
    M_PI = 'math.pi'
    M_E = 'math.e'
    NOT_A_NUMBER = 'float("nan")'
    USES_SUBSIDIARY_FILE = False
    
    binary_ops = CellMLToChasteTranslator.binary_ops.copy()
    binary_ops.update({'rem': '%'})
    nary_ops = CellMLToChasteTranslator.nary_ops.copy()
    nary_ops.update({'and': 'and', 'or': 'or'})
    function_map = {'power': 'math.pow', 'abs': 'abs', 'ln': 'math.log', 'log': 'math.log10', 'exp': 'math.exp',
                    'floor': 'math.floor', 'ceiling': 'math.ceil',
                    'factorial': 'factorial', # Needs external definition
                    'not': 'not',
                    'sin': 'math.sin', 'cos': 'math.cos', 'tan': 'math.tan',
                    'sec': '1/math.cos', 'csc': '1/math.sin', 'cot': '1/math.tan',
                    'sinh': 'math.sinh', 'cosh': 'math.cosh', 'tanh': 'math.tanh',
                    'sech': '1/math.cosh', 'csch': '1/math.sinh', 'coth': '1/math.tanh',
                    'arcsin': 'math.asin', 'arccos': 'math.acos', 'arctan': 'math.atan',
                    'arcsinh': 'math.asinh', 'arccosh': 'math.acosh', 'arctanh': 'math.atanh'}
    special_roots = {2: 'math.sqrt'}
    
    def output_file_name(self, model_filename):
        """Generate a name for our output file, based on the input file."""
        return os.path.splitext(model_filename)[0] + '.py'
    
    def open_block(self, **kwargs):
        """Just increase indent; we assume the previous line included a colon."""
        self.set_indent(offset=1)

    def close_block(self, blank_line=True, **kwargs):
        """Decrease indent, and optionally add an extra blank line."""
        self.set_indent(offset=-1)
        if blank_line:
            self.writeln(**kwargs)
        return

    def code_name(self, var, *args, **kwargs):
        """Return the full name of var in a form suitable for inclusion in a source file.
        
        Overrides the base class version to access self.parameters for parameters.
        """
        if hasattr(var, '_cml_param_index'):
            return self.vector_index(self.param_vector_name, var._cml_param_index)
        else:
            return super(CellMLToPythonTranslator, self).code_name(var, *args, **kwargs)

    def output_log(self, expr, paren):
        """Output a logarithm to the given base, which defaults to base 10."""
        if hasattr(expr, u'logbase'):
            # A base is provided.
            self.output_function('math.log', list(expr.operands()) + [expr.logbase], paren)
        else:
            # Use base 10
            self.output_function('math.log10', expr.operands(), paren)

    def output_root(self, expr, paren):
        """Output a root taken to some degree.

        If a degree qualifier element is not provided, uses default 2.
        """
        if hasattr(expr, u'degree'):
            # A degree is given.  Compute x^(1/b)
            self.write('math.pow(')
            self.output_expr(expr.operands().next(), False)
            self.write(', 1/')
            self.output_expr(expr.degree, True)
            self.write(')')
        else:
            # Compute square root
            self.output_function('math.sqrt', expr.operands(), paren)

    def output_piecewise(self, expr, paren):
        """Output the piecewise expression expr.

        We use a cascading ternary if expression for simplicity.
        """
        self.open_paren(paren)
        for piece in getattr(expr, u'piece', []):
            self.output_expr(child_i(piece, 1), True) # Result
            self.write(' if ')
            self.output_expr(child_i(piece, 2), True) # Condition
            self.write(' else ')
        if hasattr(expr, u'otherwise'):
            self.output_expr(child_i(expr.otherwise, 1), True) # Default case
        else:
            self.write(self.NOT_A_NUMBER)
        self.close_paren(paren)

    def vector_create(self, vector, size):
        """Return code for creating a new vector with the given size."""
        return ''.join(map(str, [vector, self.EQ_ASSIGN, 'np.zeros(', size, ')', self.STMT_END]))

    def vector_initialise(self, vector, size):
        """Return code for creating an already-declared vector with the given size."""
        return self.vector_create(vector, size)
    
    def analyse_model(self):
        """Figure out protocol inputs & outputs of interest, and record details as member variables."""
        assert self.use_protocol
        # Single-valued outputs
        self._outputs = cellml_metadata.find_variables(self.model,
                                                       ('pycml:output-variable', NSS['pycml']),
                                                       'yes')
        self._outputs.sort(key=lambda v: self.var_display_name(v))
        # Vector-valued outputs
        self._vector_outputs = {}
        prop = ('pycml:output-vector', NSS['pycml'])
        vector_names = set(cellml_metadata.get_targets(self.model, None,
                                                       cellml_metadata.create_rdf_node(prop)))
        for name in vector_names:
            vector_outputs = cellml_metadata.find_variables(self.model, prop, name)
            assert len(vector_outputs) > 0
            if name == 'state_variable':
                # Special case to ensure the ordering as an output matches the state vector in the ODE system
                def get_state_index(v):
                    """Find the index of the state variable corresponding to this variable, which may be units converted."""
                    v = v.get_source_variable(recurse=True)
                    if v.get_type() is VarTypes.Computed:
                        v = v.get_dependencies()[0].get_dependencies()[0]
                    return self.state_vars.index(v)
                vector_outputs.sort(key=get_state_index)
            else:
                vector_outputs.sort(key=lambda v: self.var_display_name(v))
            self._vector_outputs[name] = vector_outputs
        # Find model parameters that can be set from the protocol
        self.cell_parameters = filter(
            lambda v: v.is_modifiable_parameter,
            cellml_metadata.find_variables(self.model,
                                           ('pycml:modifiable-parameter', NSS['pycml']),
                                           'yes'))
        self.cell_parameters.sort(key=lambda v: self.var_display_name(v))
        for i, var in enumerate(self.cell_parameters):
            # Remember the var's index
            var._cml_param_index = i
        self.param_vector_name = 'self.parameters'
    
    def output_common_imports(self):
        """Output imports common to both Python and Cython code."""
        self.output_doxygen('@file\n\n',
                            'This source file was generated from CellML.\n\n',
                            'Model: ', self.model.name, '\n\n',
                            version_comment(self.add_timestamp),
                            '\n\n<autogenerated>')
        self.writeln()
        self.writeln('import numpy as np')
        self.writeln()
        self.writeln('import fc.simulations.model as Model')
        self.writeln('import fc.utility.environment as Env')
        self.writeln('import fc.language.values as V')
        self.writeln()
        
    def output_common_constructor_content(self):
        """Output __init__ content common to both Python and Cython code."""
        self.writeln('self.freeVariableName = "', self.var_display_name(self.free_vars[0]), '"')
        self.writeln('self.freeVariable = 0.0')
        self.writeln(self.vector_create('self.state', len(self.state_vars)))
        self.writeln('self.stateVarMap = {}')
        self.writeln(self.vector_create('self.initialState', len(self.state_vars)))
        input_names = set() # Check for duplicates
        for i, var in enumerate(self.state_vars):
            for name in self.get_ontology_names(var, no_names_ok=True):
                if name in input_names:
                    raise ValueError('Duplicate input variable name "' + name + '" found')
                input_names.add(name)
                self.writeln('self.stateVarMap["', name, '"] = ', i)
            init_val = getattr(var, u'initial_value', None)
            init_comm = ' # ' + var.fullname() + ' ' + var.units
            if init_val is None:
                init_comm += '; value not given in model'
                # Don't want compiler error, but shouldn't be a real number
                init_val = self.NOT_A_NUMBER
            self.writeln(self.vector_index('self.initialState', i), self.EQ_ASSIGN, init_val, init_comm)
        self.writeln()
        self.writeln('self.parameterMap = {}')
        self.writeln(self.vector_create('self.parameters', len(self.cell_parameters)))
        for var in self.cell_parameters:
            for name in self.get_ontology_names(var):
                if name in input_names:
                    raise ValueError('Duplicate input variable name "' + name + '" found')
                input_names.add(name)
                self.writeln('self.parameterMap["', name, '"] = ', var._cml_param_index)
            self.writeln(self.vector_index('self.parameters', var._cml_param_index),
                         self.EQ_ASSIGN, var.initial_value, self.STMT_END, ' ',
                         self.COMMENT_START, var.fullname(), ' ', var.units)
        # List outputs, and create objects for the GetOutputs method
        self.writeln()
        self.writeln('self.outputNames = []')
        self.writeln('outputs = self._outputs = []')
        output_names = set() # Check for duplicate local parts
        for var in self._outputs:
            # TODO: A later optimisation could look at which names the protocol actually uses, and only generate those.
            for name in self.get_ontology_names(var):
                if name in output_names:
                    raise ValueError('Duplicate output name "' + name + '" found')
                output_names.add(name)
                self.writeln('self.outputNames.append("', name, '")')
                self.writeln('outputs.append(np.array(0.0))')
        for name, vars in self._vector_outputs.iteritems():
            if name in output_names:
                raise ValueError('Duplicate output name "' + name + '" found')
            output_names.add(name)
            self.writeln('self.outputNames.append("', name, '")')
            self.writeln('outputs.append(np.array([', ', '.join(['0.0'] * len(vars)), ']))')
        self.writeln()

    def output_top_boilerplate(self):
        """Output file content occurring before the model equations."""
        self.analyse_model()
        # Start file output
        self.output_common_imports()
        self.writeln('import math')
        if self.options.numba:
            self.writeln('import numba')
            self.writeln('from numba import autojit, jit, void, double, object_')
        self.writeln()
        if self.options.numba:
            self.writeln('@jit')
        self.writeln('class ', self.class_name, '(Model.AbstractOdeModel):')
        self.open_block()
        # Constructor
        if self.options.numba:
            self.writeln('@void()')
        self.writeln('def __init__(self):')
        self.open_block()
        self.output_common_constructor_content()
        #2390 TODO: Units info
        self.writeln('Model.AbstractOdeModel.__init__(self)')
        self.close_block()
    
    def output_state_assignments(self, nodeset, stateVectorName):
        """Assign state variables used by nodeset to local names."""
        self.output_comment('State variables')
        for i, var in enumerate(self.state_vars):
            if var in nodeset:
                self.writeln(self.TYPE_CONST_DOUBLE, self.code_name(var), self.EQ_ASSIGN, self.vector_index(stateVectorName, i), self.STMT_END)
        self.writeln()

    def output_mathematics(self):
        """Output the mathematics in this model.

        This just generates the ODE right-hand side function, EvaluateRhs(self, t, y)
        """
        if self.options.numba:
            self.writeln('@jit(double[:](object_, double, double[:], double[:]))')
        self.writeln('def EvaluateRhs(self, ', self.code_name(self.free_vars[0]), ', y, ydot=np.empty(0)):')
        self.open_block()
        self.writeln('if ydot.size == 0:')
        self.writeln(self.vector_create('ydot', len(self.state_vars)), indent_offset=1)
        # Work out what equations are needed to compute the derivatives
        derivs = set(map(lambda v: (v, self.free_vars[0]), self.state_vars))
        nodeset = self.calculate_extended_dependencies(derivs)
        # Code to do the computation
        self.output_state_assignments(nodeset, 'y')
        self.output_comment('Mathematics')
        self.output_equations(nodeset)
        self.writeln()
        # Assign to derivatives vector
        for i, var in enumerate(self.state_vars):
            self.writeln(self.vector_index('ydot', i), self.EQ_ASSIGN, self.code_name(var, True), self.STMT_END)
        self.writeln('return ydot')
        self.close_block()
    
    def output_bottom_boilerplate(self):
        """Output file content occurring after the model equations, i.e. the GetOutputs method."""
        if self.options.numba:
            self.writeln('@object_()')
        self.writeln('def GetOutputs(self):')
        self.open_block()
        self.output_get_outputs_content()
        self.close_block()
    
    def get_ontology_names(self, var, no_names_ok=False):
        """Get the local names of this variable within any ontology annotations.
        
        We look at all annotations of this variable using bqbiol:is, and if any of them occur within namespaces mapped
        in the protocol, we extract the local part of the annotation URI, after the base defined by the protocol.
        Returns a list of such names, raising an error if none exist, unless no_names_ok is True.
        """
        names = []
        name_uris = var.get_rdf_annotations(('bqbiol:is', NSS['bqbiol']))
        for name_uri in name_uris:
            # Iterate through possible URI bases to find which this one is part of
            for uri_base in self.model._cml_protocol_namespaces.itervalues():
                local_part = cellml_metadata.namespace_member(name_uri, uri_base, wrong_ns_ok=True)
                if local_part:
                    names.append(local_part)
                    break # No other bases possible for this URI
        if not names and not no_names_ok:
            raise ValueError('No suitable name annotations found for variable ' + str(var))
        return names
    
    def output_get_outputs_content(self):
        """Output the content and open/close block for the GetOutputs method."""
        # Figure out what equations are needed to compute the outputs
        output_vars = set(self._outputs)
        for vars in self._vector_outputs.itervalues():
            output_vars.update(vars)
        nodeset = self.calculate_extended_dependencies(output_vars)
        # Do the calculations
        self.writeln(self.TYPE_CONST_DOUBLE, self.code_name(self.free_vars[0]), self.EQ_ASSIGN, 'self.freeVariable')
        self.output_state_assignments(nodeset, 'self.state')
        nodes_used = self.output_data_table_lookups(nodeset)
        self.output_comment('Mathematics computing outputs of interest')
        self.output_equations(nodeset - nodes_used)
        self.writeln()
        # Put the results in a list to be returned to the caller
        self.writeln('outputs = self._outputs')
        output_count = 0
        for var in self._outputs:
            # TODO: A later optimisation could look at which names the protocol actually uses, and only generate those.
            for name in self.get_ontology_names(var):
                self.writeln('outputs[', output_count, '][()] = ', self.code_name(var))
                output_count += 1
        for name, vars in self._vector_outputs.iteritems():
            for i, var in enumerate(vars):
                self.writeln('outputs[', output_count, '][', i, '] = ', self.code_name(var))
            output_count += 1
        self.writeln('return outputs')
        self.close_block()


class CellMLToCythonTranslator(CellMLToPythonTranslator):
    """Output Cython code suitable for the Python implementation of Functional Curation.
    
    Unlike the base class, code generated by this translator can't inherit from a pure Python base class.
    It also hardcodes using our Cython wrapper of CVODE as the solver.
    
    Note that we use 2 types of vector in the generated code: numpy arrays with the same names as for
    CellMLToPythonTranslator provide the same interface to the FC python code, and N_Vector views on the
    same memory provide fast access for the ODE solver.
    """

    USES_SUBSIDIARY_FILE = True
#     TYPE_VECTOR = 'cdef Sundials.N_Vector'
#     TYPE_VECTOR_REF = 'cdef Sundials.N_Vector'
    TYPE_DOUBLE = 'cdef double '
    TYPE_CONST_DOUBLE = 'cdef double '

    def output_file_name(self, model_filename):
        """Generate a name for our output file, based on the input file."""
        return os.path.splitext(model_filename)[0] + '.pyx'

    def subsidiary_file_name(self, output_filename):
        """Our subsidiary file is the setup.py used to build the extension."""
        return output_filename, os.path.join(os.path.dirname(output_filename), 'setup.py')

#     def vector_index(self, vector, i):
#         """Return code for accessing the i'th index of vector."""
#         return '(<Sundials.N_VectorContent_Serial>(' + vector + ').content).data[' + str(i) + ']'
# 
#     def vector_create(self, vector, size):
#         """Return code for creating a new vector with the given size."""
#         return ''.join(map(str, [self.TYPE_VECTOR, vector, self.EQ_ASSIGN,
#                                  'Sundials.N_VNew_Serial(', size, ')', self.STMT_END]))
# 
#     def vector_initialise(self, vector, size):
#         """Return code for creating an already-declared vector with the given size."""
#         return ''.join(map(str, [vector, self.EQ_ASSIGN, 'Sundials.N_VNew_Serial(', size, ')', self.STMT_END]))

    def output_assignment(self, expr):
        """Output an assignment statement.
        
        Avoids most of the magic in the Chaste version of this method, except for handling parameters specially.
        """
        if isinstance(expr, cellml_variable) and expr in self.cell_parameters:
            return
        return CellMLTranslator.output_assignment(self, expr)

    def output_top_boilerplate(self):
        """Output file content occurring before the model equations: basically just imports in this case.
        
        The main RHS 'method' is actually a plain function so we can use it as a C callback.
        """
        self.analyse_model()
        self.write_setup_py()
        # Start file output
        self.writeln('# cython: profile=True')
        self.output_common_imports()
        self.writeln('cimport libc.math as math')
        self.writeln('cimport numpy as np')
        self.writeln('import os')
        self.writeln('import shutil')
        self.writeln('import sys')
        self.writeln()
        self.writeln('from fc.sundials.solver cimport CvodeSolver')
        self.writeln('cimport fc.sundials.sundials as Sundials')
        self.writeln('from fc.utility.error_handling import ProtocolError')
        self.writeln()
        self.output_data_tables()

    def output_bottom_boilerplate(self):
        """Output file content occurring after the model equations, i.e. the model class."""
        base_class = 'CvodeSolver'
        self.writeln('cdef class ', self.class_name, '(', base_class, '):')
        self.open_block()
        # Declare member attributes. Note that state and _state come from the base class.
        self.writeln('cdef public char* freeVariableName')
        self.writeln('cdef public double freeVariable')
        self.writeln('cdef public object stateVarMap')
        self.writeln('cdef public np.ndarray initialState')
        self.writeln('cdef public object parameterMap')
        self.writeln('cdef public np.ndarray parameters')
        self.writeln('cdef public object outputNames')
        self.writeln()
        self.writeln('cdef public object savedStates')
        self.writeln('cdef public object env')
        self.writeln('cdef public bint dirty')
        self.writeln('cdef public char* outputPath')
        self.writeln('cdef public object indentLevel')
        self.writeln()
        self.writeln('cdef public object _module')
        self.writeln('cdef public object simEnv')
        self.writeln()
        self.writeln('cdef Sundials.N_Vector _parameters')
        self.writeln('cdef public object _outputs')
        self.writeln()
        # Constructor
        self.writeln('def __init__(self):')
        self.open_block()
        self.output_common_constructor_content()
        self.writeln('self.state = self.initialState.copy()')
        self.writeln('self.savedStates = {}')
        self.writeln('self.dirty = False')
        self.writeln('self.indentLevel = 0')
        self.writeln('self.AssociateWithModel(self)')
        self.writeln('self._parameters = Sundials.N_VMake_Serial(len(self.parameters), <Sundials.realtype*>(<np.ndarray>self.parameters).data)')
        # TODO: Use a separate environment for each ontology
        self.writeln('self.env = Env.ModelWrapperEnvironment(self)')
        # Initialise CVODE
        self.close_block()
        self.writeln('def SetRhsWrapper(self):')
        self.open_block()
        self.writeln('flag = Sundials.CVodeInit(self.cvode_mem, _EvaluateRhs, 0.0, self._state)')
        self.writeln('self.CheckFlag(flag, "CVodeInit")')
        self.close_block()
        # Cython-level destructor
        self.writeln('def __dealloc__(self):')
        self.open_block()
        self.writeln('if self._parameters != NULL:')
        self.writeln('    Sundials.N_VDestroy_Serial(self._parameters)')
        self.close_block()
        # Methods to match the AbstractModel class
        self.writeln('def SetOutputFolder(self, path):')
        self.open_block()
        self.writeln("if os.path.isdir(path) and path.startswith('/tmp'):")
        self.writeln('shutil.rmtree(path)', indent_offset=1)
        self.writeln('os.mkdir(path)')
        self.writeln('self.outputPath = path')
        self.close_block()
        self.writeln('def SetIndentLevel(self, indentLevel):')
        self.open_block()
        self.writeln('self.indentLevel = indentLevel')
        self.close_block()
        # Methods to match the AbstractOdeModel class
        self.writeln('def SetSolver(self, solver):')
        self.open_block()
        self.writeln('print >>sys.stderr, "  " * self.indentLevel, "SetSolver: Models implemented using Cython contain a built-in ODE solver, so ignoring setting."')
        self.close_block()
        self.writeln('def GetEnvironmentMap(self):')
        self.open_block()
        self.writeln('return {', nl=False)
        # TODO: Use a separate env for each ontology
        for i, prefix in enumerate(self.model._cml_protocol_namespaces.iterkeys()):
            if i > 0:
                self.write(', ')
            self.write("'%s': self.env" % prefix)
        self.writeln('}', indent=False)
        self.close_block()
        self.writeln('cpdef SetFreeVariable(self, double t):')
        self.open_block()
        self.writeln('self.freeVariable = t')
        self.writeln(base_class, '.SetFreeVariable(self, t)')
        self.close_block()
        self.writeln('def SaveState(self, name):')
        self.open_block()
        self.writeln('self.savedStates[name] = self.state.copy()')
        self.close_block()
        self.writeln('cpdef ResetState(self, name=None):')
        self.open_block()
        self.writeln('if name is None:')
        self.writeln(base_class, '.ResetSolver(self, self.initialState)', indent_offset=1)
        self.writeln('else:')
        self.writeln(base_class, '.ResetSolver(self, self.savedStates[name])', indent_offset=1)
        self.close_block()
        self.writeln('cpdef GetOutputs(self):')
        self.open_block()
        self.writeln('cdef np.ndarray[Sundials.realtype, ndim=1] parameters = self.parameters')
        self.param_vector_name = 'parameters'
        self.TYPE_CONST_UNSIGNED = 'cdef unsigned '
        self.output_get_outputs_content()
        del self.TYPE_CONST_UNSIGNED
        self.param_vector_name = 'self.parameters'
        self.close_block()
    
    def output_mathematics(self):
        """Output the mathematics in this model.

        This generates the ODE right-hand side function, "EvaluateRhs(self, t, y)", but as a C-style callback for CVODE.
        """
        self.writeln('cdef int _EvaluateRhs(Sundials.realtype ', self.code_name(self.free_vars[0]),
                     ', Sundials.N_Vector y, Sundials.N_Vector ydot, void* user_data):')
        self.open_block()
        self.writeln('model = <object>user_data')
        self.writeln('cdef np.ndarray[Sundials.realtype, ndim=1] parameters = <np.ndarray>model.parameters')
        self.param_vector_name = 'parameters'
        # Work out what equations are needed to compute the derivatives
        derivs = set(map(lambda v: (v, self.free_vars[0]), self.state_vars))
        nodeset = self.calculate_extended_dependencies(derivs)
        # Code to do the computation
        self.output_comment('State variables')
        for i, var in enumerate(self.state_vars):
            if var in nodeset:
                self.writeln(self.TYPE_DOUBLE, self.code_name(var), ' = (<Sundials.N_VectorContent_Serial>y.content).data[', i, ']')
        self.writeln()
        self.TYPE_CONST_UNSIGNED = 'cdef unsigned '
        nodes_used = self.output_data_table_lookups(nodeset)
        del self.TYPE_CONST_UNSIGNED
        self.writeln()
        self.output_comment('Mathematics')
        self.output_equations(nodeset - nodes_used)
        self.writeln()
        # Assign to derivatives vector
        for i, var in enumerate(self.state_vars):
            self.writeln('(<Sundials.N_VectorContent_Serial>ydot.content).data[', i, '] = ', self.code_name(var, True))
        self.param_vector_name = 'self.parameters'
        self.close_block()

    def output_array_definition(self, array_name, array_data):
        """Output code to create and fill a fixed-size 1d array."""
        self.writeln('cdef Sundials.realtype[', len(array_data), ']', array_name)
        self.writeln(array_name, '[:] = [', ', '.join(map(lambda f: "%.17g" % f, array_data)), ']')

    def fast_floor(self, arg):
        """Return code to compute the floor of an argument as an integer quickly, typically by casting."""
        return "int(%s)" % arg

    def write_setup_py(self):
        """Write our subsidiary setup.py file for building the extension."""
        self.out2.write("""
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

ext_modules=[
    Extension("%(filebase)s",
              ["%(filebase)s.pyx"],
              include_dirs=[numpy.get_include(), '%(fcpath)s'],
              #library_dirs=['%(fcpath)s/fc/sundials'],
              libraries=['sundials_cvode', 'sundials_nvecserial', 'm'])
              # users can set CFLAGS and LDFLAGS in their env if needed
]

setup(
  name = "%(filebase)s",
  cmdclass = {"build_ext": build_ext},
  ext_modules = ext_modules
)
""" % {'filebase': os.path.splitext(os.path.basename(self.output_filename))[0],
       'fcpath': os.path.join(os.path.dirname(__file__), '../../projects/FunctionalCuration/src/python')})



###############################################
# Register translation classes in this module #
###############################################

CellMLTranslator.register(CellMLTranslator, 'C++')
CellMLTranslator.register(CellMLToChasteTranslator, 'Chaste')
CellMLTranslator.register(CellMLToCvodeTranslator, 'CVODE')
CellMLTranslator.register(CellMLToMapleTranslator, 'Maple')
CellMLTranslator.register(CellMLToMatlabTranslator, 'Matlab')
CellMLTranslator.register(CellMLToHaskellTranslator, 'Haskell')
CellMLTranslator.register(CellMLToPythonTranslator, 'Python')
CellMLTranslator.register(CellMLToCythonTranslator, 'Cython')




class SolverInfo(object):
    """Add information for specialised translator classes into a model."""
    def __init__(self, model, force=False):
        """Add information for the solvers as XML.

        The Jacobian and linearity analyses store their results in
        Python data structures as attributes of this object.
        Transcribe these into XML in a child <solver_info> element.

        If any of these elements exist in the model they will be left
        unaltered, unless force is set to True.
        
        This constructor just sets up the container element; call one
        of the add_* methods to actually add the information to it.
        """
        self._model = model
        if force and hasattr(model, u'solver_info'):
            model.xml_remove_child(model.solver_info)
        if hasattr(model, u'solver_info'):
            solver_info = model.solver_info
        else:
            solver_info = model.xml_create_element(u'solver_info', NSS[u'solver'])
            model.xml_append(solver_info)
        self._solver_info = solver_info
        self._component = None
        self._dt = None
    
    def add_all_info(self):
        """Actually add the info."""
        self.add_transmembrane_potential_name()
        self.add_membrane_ionic_current()
        self.add_linearised_odes()
        self.add_jacobian_matrix()
        self.add_dt_reference()
    
    def add_dt_reference(self):
        """Add a reference to the variable representing dt."""
        solver_info = self._solver_info
        model = self._model
        if not hasattr(solver_info, u'dt'):
            dt = self.get_dt()
            elt = model.xml_create_element(u'dt', NSS[u'solver'], content=dt.fullname(cellml=True))
            solver_info.xml_append(elt)
            self._model._add_sorted_assignment(dt)
    
    def add_transmembrane_potential_name(self):
        """The name of the transmembrane potential."""
        solver_info = self._solver_info
        model = self._model
        if not hasattr(solver_info, u'transmembrane_potential'):
            v_elt = model.xml_create_element(
                u'transmembrane_potential', NSS[u'solver'],
                content=model._cml_transmembrane_potential.fullname())
            solver_info.xml_append(v_elt)
    
    def add_linearised_odes(self):
        """Linearised ODEs - where du/dt = g + hu (and g, h are not functions of u).
        
        Structure looks like:
        <linear_odes>
            <math>
                <apply><eq/>
                    <apply><diff/>
                        <bvar><ci>t</ci></bvar>
                        <ci>u</ci>
                    </apply>
                    <apply><plus/>
                        g
                        <apply><times/>
                            h
                            <ci>u</ci>
                        </apply>
                    </apply>
                </apply>
                .
                .
                .
            </math>
        </linear_odes>
        """
        solver_info = self._solver_info
        model = self._model
        if not hasattr(solver_info, u'linear_odes') and model._cml_linear_update_exprs:
            odes_elt = model.xml_create_element(u'linear_odes', NSS[u'solver'])
            solver_info.xml_append(odes_elt)
            odes_math = model.xml_create_element(u'math', NSS[u'm'])
            odes_elt.xml_append(odes_math)
            linear_vars = model._cml_linear_update_exprs.keys()
            linear_vars.sort(key=lambda v: v.fullname())
            free_var = model._cml_free_var
            for var in linear_vars:
                g, h = model._cml_linear_update_exprs[var]
                hu = mathml_apply.create_new(model, u'times', [h, var.fullname()])
                rhs = mathml_apply.create_new(model, u'plus', [g, hu])
                odes_math.xml_append(mathml_diff.create_new(
                    model, free_var.fullname(), var.fullname(), rhs))
            # Ensure that the model has a special component
            self._get_special_component()
    
    def _fix_jac_var_name(self, vname):
        """
        If PE will be performed on a model with a single component, then we'll need full names in
        the variable attributes.
        """
        if vname[:4] == 'var_' and len(self._model.component) == 1 and not self._model.component.ignore_component_name:
            name = unicode('var_' + self._model.component.name + '__' + vname[4:])
        else:
            name = unicode(vname)
        return name
    
    def add_jacobian_matrix(self):
        """Jacobian matrix elements.
        
        Structure looks like:
        <jacobian>
            [<math> assignments of common sub-terms </math>]
            <entry var_i='varname' var_j='varname'>
                <math> apply|cn|ci ...</math>
            </entry>
        </jacobian>
        """
        solver_info = self._solver_info
        model = self._model
        if model._cml_jacobian and model._cml_jacobian_full:
            jac = model._cml_jacobian[1]
        else:
            # Old-style partial jacobian, or no jacobian
            jac = model._cml_jacobian
        if not hasattr(solver_info, u'jacobian') and jac:
            jac_elt = model.xml_create_element(u'jacobian', NSS[u'solver'])
            solver_info.xml_append(jac_elt)
            
            if model._cml_jacobian_full:
                # There may be temporaries
                temporaries = model._cml_jacobian[0]
                if temporaries:
                    jac_elt.xml_append(amara_parse_cellml(temporaries).math)

            jac_vars = jac.keys()
            jac_vars.sort() # Will sort by variable name
            for v_i, v_j in jac_vars:
                # Add (i,j)-th entry
                attrs = {u'var_i': self._fix_jac_var_name(v_i),
                         u'var_j': self._fix_jac_var_name(v_j)}
                entry = model.xml_create_element(u'entry', NSS[u'solver'], attributes=attrs)
                jac_elt.xml_append(entry)
                entry_doc = amara_parse_cellml(jac[(v_i, v_j)].xml())
                entry.xml_append(entry_doc.math)
            # Ensure that the model has a special component
            self._get_special_component()
        return
    
    def use_canonical_variable_names(self):
        """
        PE has just been performed, so we need to update variable names occurring outside
        the modifiable mathematics sections.
        """
        jac_elt = getattr(self._solver_info, u'jacobian', None)
        for entry in getattr(jac_elt, u'entry', []):
            for vlabel in ['var_i', 'var_j']:
                vname = getattr(entry, vlabel)
                var = self._get_variable(vname)
                new_name = var.get_source_variable(recurse=True).fullname()
                setattr(entry, vlabel, new_name)
        dt_elt = getattr(self._solver_info, u'dt', None)
        if dt_elt:
            var = self._get_variable(unicode(dt_elt))
            new_name = var.get_source_variable(recurse=True).fullname()
            dt_elt.xml_remove_child(unicode(dt_elt))
            dt_elt.xml_append(unicode(new_name))

    def add_membrane_ionic_current(self):
        """Add ionic current information as XML for solvers to use."""
        solver_info = self._solver_info
        model = self._model
        # The total ionic current.  This relies on having a configuration store.
        if hasattr(model.xml_parent, '_cml_config') and not hasattr(solver_info, u'ionic_current'):
            conf = model.xml_parent._cml_config
            if conf.i_ionic_vars:
                ionic_elt = model.xml_create_element(u'ionic_current', NSS[u'solver'])
                # Adds each ionic var to the xml doc from the config store
                for var in conf.i_ionic_vars:
                    varelt = model.xml_create_element(u'var', NSS[u'solver'],
                                                      content=var.fullname())
                    ionic_elt.xml_append(varelt)
                solver_info.xml_append(ionic_elt)
        return
    
    def add_linear_ode_update_equations(self):
        """Add the update equations for the linear ODEs.
        
        A linear ODE has the form du/dt = g+h.u where g & h are not functions of u.  The
        update expression then looks like u = (u + g.dt)/(1 - h.dt).
        
        This replaces the linear_odes block with the structure:
        <linear_odes>
            <math>
                <ci>u</ci>
                <ci>t</ci>
                <apply> <!-- (u + g.dt)/(1 - h.dt) --> </apply>
            </math>
            .
            .
            .
        </linear_odes>
        """
        block = getattr(self._solver_info, u'linear_odes', None)
        dt = self._model.get_config().dt_variable.fullname() # was dt = u'delta_t'
        # Add the new equations
        for u, t, gh in self.get_linearised_odes():
            g, h = gh
            g.safe_remove_child(g, g.xml_parent)
            g_dt = mathml_apply.create_new(block, u'times', [g, dt])
            numer = mathml_apply.create_new(block, u'plus', [u.fullname(), g_dt])
            h.safe_remove_child(h, h.xml_parent)
            h_dt = mathml_apply.create_new(block, u'times', [h, dt])
            denom = mathml_apply.create_new(block, u'minus', [(u'1', u'dimensionless'), h_dt])
            eqn = mathml_apply.create_new(block, u'divide', [numer, denom])
            math = block.xml_create_element(u'math', NSS[u'm'])
            math.xml_append(mathml_ci.create_new(block, u.fullname()))
            math.xml_append(mathml_ci.create_new(block, t.fullname()))
            math.xml_append(eqn)
            block.xml_append(math)
            self._add_variable_links(math)
        # Remove the old equations (first math element)
        block.xml_remove_child(block.math)
    
    def add_variable_links(self):
        """Link ci elements in the added XML to cellml_variable objects.
        
        This analyses the names in the ci elements to determine which variable in
        the model they refer to.
        """
        self._process_mathematics(self._add_variable_links)
        #1795 - classify temporary variables for the Jacobian matrix, and append
        # to the main list of assignments in the model
        solver_info = self._solver_info
        if hasattr(solver_info, u'jacobian') and hasattr(solver_info.jacobian, u'math'):
            for elt in solver_info.jacobian.math.apply:
                elt.classify_variables(root=True)
            for elt in solver_info.jacobian.math.apply:
                self._model.topological_sort(elt)
        #2418 - check if any state variables have been units-converted
        self._check_state_var_units_conversions()
    
    def _check_state_var_units_conversions(self):
        """Check if any Jacobian entries need to be altered because the units of state variables have changed.
        
        If any variable considered a state variable by the Jacobian is now of type Computed then it has been
        converted.  We figure out the conversion factor, update the Jacobian to reference the new state variable,
        and units-convert the derivative.
        """
        if not hasattr(self._solver_info, u'jacobian'):
            return
        # Helper methods
        def set_var_values(elt, vars=None):
            """Fake all variables appearing in the given expression being set to 1.0, and return them."""
            if vars is None:
                vars = []
            if isinstance(elt, mathml_ci):
                elt.variable.set_value(1.0)
                vars.append(elt.variable)
            else:
                for child in getattr(elt, 'xml_children', []):
                    set_var_values(child, vars)
            return vars
        # Find any converted state variables
        converted_state_vars = set()
        for entry in getattr(self._solver_info.jacobian, u'entry', []):
            var = self._get_variable(entry.var_i)
            if var.get_type() == VarTypes.Computed:
                converted_state_vars.add(var)
        if not converted_state_vars:
            return
        # Figure out the conversion factor in each case
        state_var_map = {}
        for var in converted_state_vars:
            defn = var.get_dependencies()[0]
            defn_vars = set_var_values(defn.eq.rhs)
            assert len(defn_vars) == 1, "Unexpected form of units conversion expression found"
            factor = defn.eq.rhs.evaluate()
            state_var_map[var] = (defn_vars[0], factor)
            defn_vars[0].unset_values()
        # Apply the conversion to relevant Jacobian entries
        for entry in getattr(self._solver_info.jacobian, u'entry', []):
            factor = 1
            var_i = self._get_variable(entry.var_i)
            if var_i in converted_state_vars:
                var_i, factor_i = state_var_map[var_i]
                var_i = var_i.get_source_variable(recurse=True)
                entry.var_i = unicode(var_i.fullname())
                factor /= factor_i
            var_j = self._get_variable(entry.var_j)
            if var_j in converted_state_vars:
                var_j, factor_j = state_var_map[var_j]
                var_j = var_j.get_source_variable(recurse=True)
                entry.var_j = unicode(var_j.fullname())
                factor *= factor_j
            if factor != 1:
                # Replace rhs with rhs * factor
                rhs = list(entry.math.xml_element_children())[0]
                entry.math.safe_remove_child(rhs)
                new_rhs = mathml_apply.create_new(entry, 'times', [(factor, 'dimensionless'), rhs])
                entry.math.xml_append(new_rhs)
    
    def do_binding_time_analysis(self):
        """Do a binding time analysis on the additional mathematics.
        
        This requires self.add_variable_links to have been called already.
        """
        self._process_mathematics(lambda elt: elt._get_binding_time())
        
    def _process_mathematics(self, func):
        """Apply func to each top-level mathematical construct in the solver info blocks.
        
        func must be able to accept mathml_piecewise, mathml_apply, mathml_ci and mathml_cn elements.
        """
        solver_info = self._solver_info
        # Jacobian
        if hasattr(solver_info, u'jacobian'):
            if hasattr(solver_info.jacobian, u'math'):
                for elt in solver_info.jacobian.math.apply:
                    func(elt)
            for entry in solver_info.jacobian.entry:
                for elt in entry.math.xml_element_children():
                    func(elt)
        # Linearised ODEs
        if hasattr(solver_info, u'linear_odes'):
            for math in solver_info.linear_odes.math:
                for elt in math.xml_element_children():
                    func(elt)
    
    def has_modifiable_mathematics(self):
        """Check if the solver info blocks contain any modifiable mathematics."""
        try:
            self.get_modifiable_mathematics().next()
            return True
        except StopIteration:
            return False
    
    def get_modifiable_mathematics(self):
        """Get an iterable over mathematical constructs in the solver info blocks that can be changed.
        
        Returned elements will be mathml_piecewise, mathml_apply, mathml_ci or mathml_cn instances.
        """
        solver_info = self._solver_info
        # Jacobian - entry definitions and temporaries can be changed
        if hasattr(solver_info, u'jacobian'):
            if hasattr(solver_info.jacobian, u'math'):
                for elt in solver_info.jacobian.math.apply:
                    yield elt
            for entry in solver_info.jacobian.entry:
                for elt in entry.math.xml_element_children():
                    yield elt
        # Linearised ODEs - only g & h can be changed
        if hasattr(solver_info, u'linear_odes'):
            for _, _, eqns in self.get_linearised_odes():
                for eqn in eqns:
                    yield eqn

    def get_linearised_odes(self):
        """Return an iterable over the linearised ODEs, i.e. ODEs of the form
        du/dt = g + hu (with g, h not functions of u).
        
        Yields tuples (u, t, eqns) where the form of eqns depends on whether
        add_linear_ode_update_equations has been called.  If so, it is a 1-tuple
        containing the update equation; if not, it is (g,h).
        """
        if hasattr(self._solver_info, u'linear_odes'):
            if hasattr(self._solver_info.linear_odes.math, u'ci'):
                for math in self._solver_info.linear_odes.math:
                    u, t, eqn = list(math.xml_element_children())
                    u = u.variable
                    t = t.variable
                    yield (u, t, (eqn,))
            else:
                for ode in self._solver_info.linear_odes.math.apply:
                    u = ode.apply.ci.variable
                    t = ode.apply.bvar.ci.variable
                    opers = ode.apply[1].operands()
                    g = opers.next()
                    h = opers.next().operands().next()
                    yield (u, t, (g,h))
    
    def _add_variable_links(self, elt):
        """Recursively link ci elements in the given XML tree to cellml_variable objects.
        
        Also sets component links: for ci elements, to the component containing the linked
        variable, and for cn elements, to the first component in the model.
        """
        if isinstance(elt, mathml_ci):
            var = self._get_variable(unicode(elt))
            elt._cml_variable = var
            elt._cml_component = var.component
        elif isinstance(elt, mathml_cn):
            # Fake a component, since it doesn't really have one
            elt._cml_component = elt.model.component
        elif hasattr(elt, 'xml_children'):
            for child in elt.xml_children:
                self._add_variable_links(child)

    _jac_temp_re = re.compile(r't[0-9]+')
    def _get_variable(self, varname):
        """Return the variable in the model with name varname."""
        try:
            if varname == 'delta_t':
                # Special case for the timestep in ComputeJacobian and elsewhere
                var = self.get_dt()
            elif self._jac_temp_re.match(varname):
                var = self._get_special_variable(varname, VarTypes.Unknown)
            else:
                var = cellml_variable.get_variable_object(self._model, varname)
        except KeyError:
            raise ValueError("Cannot find variable '%s' referenced in SolverInfo" % varname)
        return var
    
    def create_dt(self, modifier, comp, units):
        """Create the special 'dt' variable in the given component."""
        self._dt = modifier.add_variable(comp, modifier._uniquify_var_name(u'dt', comp), units)
        self._dt._set_type(VarTypes.Free)
        return self._dt
    
    def get_dt(self):
        """Get or create a special 'dt' variable."""
        if not self._dt:
            self._dt = self._get_special_variable(u'dt', VarTypes.Free)
        return self._dt

    def _get_special_variable(self, varname, ptype=VarTypes.Unknown):
        """Get or create a special variable object that doesn't really exist in the model."""
        comp = self._get_special_component()
        try:
            var = comp.get_variable_by_name(varname)
        except KeyError:
            var = cellml_variable.create_new(self._model, varname, u'dimensionless')
            comp._add_variable(var)
            var._set_type(ptype)
        return var

    def _get_special_component(self):
        """Get or create a special component for containing special variables."""
        if not self._component:
            self._component = cellml_component.create_new(self._model, u'')
            self._model._add_component(self._component, special=True)
        return self._component



class ConfigurationStore(object):
    """
    A container for configuration information, read in from XML
    configuration files.  The file structure is described in the
    read_configuration_file method.
    """
    def __init__(self, doc, options=None):
        """Create a new store.

        doc specifies a CellML document, the processing of which this configuration store will affect.

        If given, options should be an optparse.Values instance containing command-line options.
        """
        self.doc = doc
        doc._cml_config = self
        self.options = options
        self.unit_definitions = cellml_component.create_new(doc.model, '*lookup_table_units*')
        self.unit_definitions.xml_parent = doc.model # Needed for looking up standard units
        # Transmembrane potential
        self.V_definitions = [u'membrane,V']
        self.V_variable = None
        # Membrane capacitance
        self.Cm_definitions = []
        self.Cm_variable = None
        # Cytosolic calcium concentration
        self.cytosolic_calcium_concentration = None
        # Lookup table configuration
        self.lut_config = {}
        # Ionic currents configuration
        self.i_stim_definitions = [u'membrane,i_Stim']
        self.i_stim_var = None
        self.i_ionic_definitions = [u'membrane,i_.*']
        self.i_ionic_vars = []
        # Whether GetIIonic will need to negate the sum of i_ionic_vars
        self.i_ionic_negated = False
        # Whether the stimulus magnitude is positive, rather than negative
        self.i_stim_negated = False
        # Other variables that may be set by other code, for example an InterfaceGenerator
        self.dt_variable = None
        self.i_data_clamp_current = None
        self.i_data_clamp_conductance = None
        return

    def read_configuration_file(self, config_file):
        """Read configuration stored in config_file.

        The configuration file is expected to be XML, conforming to
        the following structure.  Currently little checking is done on
        the file format; incorrectly formatted files are unlikely to
        give particularly helpful error messages.

        The root element may contain a 'global' element, giving global
        configuration options.  These include:

         * 'lookup_tables'
           Contains one or more 'lookup_table' elements, one for each
           type of lookup table available.  These contain (a selection of)
           the elements:
           * 'var' - the variable to key on.  The component name given
             should be that from which the variable is exported.  Must be
             present.
           * 'min', 'max', 'step' - table bounds parameters.  Optional.
           Default values are used for parameters that are not present.
         * 'currents'
           Defines which variables hold the ionic and stimulus currents,
           if any.  It contains 2 elements:
           * 'stimulus' - the full name of the stimulus current variable
           * 'ionic_match' - a regular expression matching full names of
             ionic current variables.  It may also match the stimulus
             current, but the stimulus will never be considered an ionic
             current.  The value is split on ','; the first part is then
             matched against component names, and the second against
             variables in matching components.
             
             This is mostly redundant now, because the equation for dV/dt
             is used first to determine the ionic currents (see documentation
             for _find_transmembrane_currents_from_voltage_ode), and only
             if this fails to find suitable currents will the ionic_match
             definition be used.
         * 'transmembrane_potential'
           Defines which variable holds the transmembrane potential.
           Defaults to 'membrane,V' if not present.
         * 'membrane_capacitance'
           Defines which variable holds the cell membrane capacitance.
           
        The root element also contains 0 or more 'for_model' elements,
        which contain settings for individual models.  These must have
        at least one of an 'id' or 'name' attribute, which specify the
        model in question.  They can also contain anything allowable as
        global configuration options.  Options given here override those
        specified globally.

        Configuration which is identical for groups of models may be given
        using the 'for_models' element.  This has an 'ids' element as its
        first child, which contains 'id' elements holding either the name
        or id of a model.  The remaining contents of the 'for_models'
        element are as for 'for_model'.

        There are 3 ways of specifying variables:
        1. By name (var type='name')
           Variable names are given in full form, i.e. component_name,variable_name
        2. By standardised name (var type='oxmeta')
           Use the name from the oxmeta annotations
        3. By reference to a section of the config file (when defining lookup table keys)
           e.g. <var type='config-name'>transmembrane_potential</var>

        Within any element that specifies a variable, a list of <var> elements can be
        provided.  Each will be tried in turn to see if a match can be found in the model,
        and the first match wins.

        Some items are overridden if oxmeta annotations are present in the model, with
        the annotated variable taking precedence over the config file specification.
        """
        DEBUG('config', "Reading configuration from ", config_file)
        binder = amara.bindery.binder()
        binder.set_binding_class(None, "units", cellml_units)
        binder.set_binding_class(None, "unit", cellml_unit)
        rules = [bt.ws_strip_element_rule(u'*')]
        config_doc = amara_parse(config_file, rules=rules, binderobj=binder)
        # Store extra units definitions
        for defn in config_doc.xml_xpath(u'/*/units'):
            defn.xml_parent = self.unit_definitions # Needed for looking up the units this definition is derived from
            self.unit_definitions.add_units(defn.name, defn)
        # Overrides for command-line options
        if self.options and hasattr(config_doc.pycml_config, 'command_line_args'):
            args = map(str, config_doc.pycml_config.command_line_args.arg)
            args.append('dummy-file')
            get_options(args, self.options)
        # Sections to use in configuration; later sections take precedence
        sections = []
        # Use global configuration?
        glo = config_doc.xml_xpath(u'/*/global')
        if glo:
            sections.append(glo[0])
        # Get the config section(s) for our model.  Sections
        # specifically for this model come after sections covering
        # multiple models, so they take precedence.
        model_id = getattr(self.doc.model, u'id', self.doc.model.name)
        sections.extend(config_doc.xml_xpath(
            u'/*/for_models[ids/id="%s" or ids/id="%s"]'
            % (self.doc.model.name, model_id)))
        sections.extend(config_doc.xml_xpath(
            u'/*/for_model[@name="%s" or @id="%s"]'
            % (self.doc.model.name, model_id)))
        # Main items of configuration
        for section in sections:
            if hasattr(section, u'lookup_tables'):
                self._parse_lookup_tables(section.lookup_tables)
            if hasattr(section, u'currents'):
                self._parse_currents(section.currents)
            if hasattr(section, u'transmembrane_potential'):
                self._parse_Vm(section.transmembrane_potential)
            if hasattr(section, u'membrane_capacitance'):
                self._parse_Cm(section.membrane_capacitance)
    
    def finalize_config(self):
        """Having read all the configuration files, apply to the model."""
        # If no LT options given, add defaults
        if not self.lut_config:
            config_key = ('config-name', 'transmembrane_potential')
            self.lut_config[config_key] = {}
            self._set_lut_defaults(self.lut_config[config_key])
        # Identify the variables in the model
        self.find_transmembrane_potential()
        self.find_membrane_capacitance()
        self.find_cytosolic_calcium_concentration()
        if not self.options.protocol:
            self.find_current_vars()

    def _create_var_def(self, content, defn_type):
        """Create a variable definition object."""
        xml_fragment = '<var type="%s">%s</var>' % (defn_type, content)
        return amara.parse(str(xml_fragment)).var

    def _check_var_def(self, var_elt, var_desc):
        """Check a variable definition is syntactically valid.
        
        If type == 'name', it must have text content of the form "component_name,variable_name".
        If type == 'oxmeta', it must have text content drawn from METADATA_NAMES.
        If type == 'config-name', it must have text content either 'stimulus' or 'transmembrane_potential'.
        """
        defn_type = getattr(var_elt, u'type', u'name')
        if defn_type == u'name':
            name_parts = unicode(var_elt).strip().split(',')
            if len(name_parts) != 2:
                raise ConfigurationError('Invalid definition of ' + var_desc + ': '
                                         + unicode(var_elt))
        elif defn_type == u'oxmeta':
            if unicode(var_elt) not in cellml_metadata.METADATA_NAMES:
                raise ConfigurationError('"' + unicode(var_elt) + '" is not a valid oxmeta name')
        elif defn_type == u'config-name':
            if unicode(var_elt) not in [u'stimulus', u'transmembrane_potential', u'membrane_capacitance']:
                raise ConfigurationError('"' + unicode(var_elt) + '" is not a name known to the config file')
        else:
            raise ConfigurationError('"' + defn_type + '" is not a valid variable definition type')
        return

    def _parse_var(self, elt, name):
        """Parse definition of a special variable."""
        if hasattr(elt, 'var'):
            # List of possibilities
            defs = []
            for vardef in elt.var:
                self._check_var_def(vardef, name)
                defs.append(vardef)
        else:
            # Old style - single variable given by text content
            self._check_var_def(elt, name)
            defs = [elt]
        return defs

    def _parse_Vm(self, vm_elt):
        """Parse definition of variable holding the transmembrane potential."""
        self.V_definitions = self._parse_var(vm_elt, 'transmembrane_potential')
    
    def _parse_Cm(self, cm_elt):
        """Parse definition of variable holding the cell membrane capacitance."""
        self.Cm_definitions = self._parse_var(cm_elt, 'membrane_capacitance')

    def _parse_currents(self, currents):
        """Parse definitions of ionic and stimulus currents."""
        if hasattr(currents, u'stimulus'):
            self.i_stim_definitions = self._parse_var(currents.stimulus, 'stimulus current')
        if hasattr(currents, u'ionic_match'):
            self.i_ionic_definitions = self._parse_var(currents.ionic_match, 'ionic currents')
        return
    
    def _find_variable(self, defn, pe_done=False):
        """Find a variable matching the given definition.

        If pe_done is True, then partial evaluation has been performed
        on the model, so looking for variables by name needs to look for
        variables called compname__varname in the single component.
        """
        defn_type = getattr(defn, u'type', u'name')
        if defn_type == u'name':
            name_parts = unicode(defn).strip().split(',')
            if pe_done:
                try:
                    var = self.doc.model.component.get_variable_by_name(u'__'.join(name_parts))
                except KeyError:
                    var = None
            else:
                var = self.doc.model.xml_xpath(u'cml:component[@name="%s"]/cml:variable[@name="%s"]'
                                               % tuple(name_parts))
                if var:
                    var = var[0]
        elif defn_type == u'oxmeta':
            var = self.doc.model.get_variable_by_oxmeta_name(str(defn), throw=False)
        elif defn_type == u'config-name':
            if unicode(defn) == u'stimulus':
                var = self.i_stim_var
            elif unicode(defn) == u'transmembrane_potential':
                var = self.V_variable
            elif unicode(defn) == u'membrane_capacitance':
                var = self.Cm_variable
            else:
                raise ConfigurationError('"' + str(defn) + '" is not a valid configuration file variable name')
        else:
            raise ConfigurationError('"' + defn_type + '" is not a valid variable definition type')
        return var
    
    def _process_ci_elts(self, elt, func, **kwargs):
        """Recursively apply func to any ci elements in the tree rooted at elt."""
        if isinstance(elt, mathml_ci):
            func(elt, **kwargs)
        else:
            for child in getattr(elt, 'xml_children', []):
                self._process_ci_elts(child, func, **kwargs)
    
    def _find_transmembrane_currents_from_voltage_ode(self):
        """Analyse the expression for dV/dt to determine the transmembrane currents.
        
        Looks for an equation defining dV/d(something) and assumes the something is
        time; this will be checked during code generation for Chaste.  It then finds
        all variables on the RHS of this equation which have the same units as the
        stimulus current (self.i_stim_var) and identifies these as transmembrane
        currents.  Will automatically exclude the stimulus current.
        
        If self.V_variable is not set, returns the empty list.
        """
        if not self.V_variable:
            DEBUG('config', "Transmembrane potential not configured, so can't determine currents from its ODE")
            return []
        if self.i_stim_var:
            current_units = [self.i_stim_var.component.get_units_by_name(self.i_stim_var.units)]
        else:
            current_units = CellMLToChasteTranslator.get_current_units_options(self.doc.model)
        ionic_vars = []
        
        def find_units_match(test_units, units_list, remove_match=False, keep_only_match=False):
            """Look for a units definition dimensionally equivalent to test_units within units_list.
            
            If remove_match is True, remove the first match from the list.
            If keep_only_match is True, remove all entries except the first match from the list.
            Return the matching units, or None if there are no matches.
            """
            for units in units_list:
                if test_units.dimensionally_equivalent(units):
                    match = units
                    break
            else:
                match = None
            if match and remove_match:
                units_list.remove(match)
            if match and keep_only_match:
                units_list[:] = [match]
            return match

        def clear_values(expr, process_definitions=False):
            """Recursively clear saved values for variables in this expression.
            
            If process_definitions is True, recursively treat expressions defining variables
            used in this expression, too.
            """
            def process_var(var):
                var.unset_values()
                var._unset_binding_time(only_temporary=True)
                if process_definitions:
                    defn = var.get_dependencies()
                    if defn:
                        if isinstance(defn[0], mathml_apply):
                            clear_values(defn[0].eq.rhs, process_definitions=True)
                        elif isinstance(defn[0], cellml_variable):
                            process_var(defn[0])
            def process_ci(ci_elt):
                process_var(ci_elt.variable)
            self._process_ci_elts(expr, process_ci)
        
        def check_if_current(ci_elt, vars_found):
            """Check if this is a transmembrane current."""
            v = ci_elt.variable
            if v.get_source_variable(recurse=True) is not self.i_stim_var:
                vars_found.append(v)
                # Check units
                u = v.component.get_units_by_name(v.units)
                if find_units_match(u, current_units, keep_only_match=True):
                    ionic_vars.append(v.get_source_variable(recurse=True))
                    ionic_vars[-1]._cml_ref_in_dvdt = ci_elt # Hack for data clamp support (#2708)
            # Fake this variable being 1 so we can check the sign of GetIIonic
            if not v.is_statically_const(ignore_annotations=True):
                v.set_value(1.0)
        
        def bfs(func, vars, *args, **kwargs):
            """Do a breadth first search of the definitions of variables in vars.
            
            func is the recursive function to call.  It will be given the list of defining expressions
            as its first argument, and args and kwargs as remaining arguments.
            """
            def get_defn(var):
                defn = var.get_dependencies()
                if defn:
                    var._set_binding_time(BINDING_TIMES.static, temporary=True)
                    if isinstance(defn[0], cellml_variable):
                        defn = get_defn(defn[0])
                    else:
                        assert isinstance(defn[0], mathml_apply)
                        var.unset_values()
                        defn = defn[0].eq.rhs
                return defn
            defns = []
            for var in vars:
                defn = get_defn(var)
                if defn:
                    defns.append(defn)
            if defns:
                func(defns, *args, **kwargs)

        def find_currents(exprs, depth=0, maxdepth=2):
            """Find ionic currents by searching the given expressions.
            
            On the initial call, exprs should contain just the definition of dV/dt (i.e. the RHS).
            
            Uses breadth-first search of the equation dependency tree to find variables that
            have units dimensionally equivalent to one of the current formulations that Chaste
            can handle, or equivalent to the stimulus current's units if one is defined.
            
            Initially, A_per_F is removed from the list, since the RHS of dV/dt should always
            have equivalent dimensions.  If another option can't be found within maxdepth levels,
            we restart the search with A_per_F included.  The depth limit is intended to guard against
            unexpectedly finding something that isn't a current; it's somewhat dodgy, but won't
            break on any model I know, and I haven't thought of a better approach yet.
            
            When one variable with suitable units is found, further ionic currents must have units
            equivalent to its to be found.  Also once one ionic current is found, only the remaining
            expressions at its depth will be processed.
            """
            if depth == 0 and maxdepth > 0:
                dvdt_units = exprs[0].xml_parent.eq.lhs.get_units()
                A_per_F = find_units_match(dvdt_units, current_units, remove_match=True)
#                 # We could do this check, but actually it doesn't catch much and later checks will pick up on the problem
#                 if A_per_F is None and not self.i_stim_var:
#                     raise ConfigurationError('Units ' + dvdt_units.description() + ' of dV/dt are not equivalent to V/s - unable to continue.')
            # Process all expressions at this depth
            vars_found = []
            for expr in exprs:
                self._process_ci_elts(expr, check_if_current, vars_found=vars_found)
            if not ionic_vars and depth != maxdepth:
                # Process the definitions of expressions at this depth
                bfs(find_currents, vars_found, depth+1, maxdepth)
            # If we reached maxdepth unsuccessfully, try again with A_per_F included (if it was an option)
            if not ionic_vars and depth == 0 and maxdepth > 0 and A_per_F:
                current_units.append(A_per_F)
                find_currents(exprs, depth, maxdepth=-1)

        def assign_values_for_stimulus_check(exprs, found_stim=Sentinel()):
            """Assign temporary values to variables in order to check the stimulus sign.
            
            This will process defining expressions in a breadth first search until the stimulus
            current is found.  Each variable that doesn't have its definitions processed will
            be given a value as follows:
             - stimulus current = 1
             - other currents = 0
             - other variables = 1
            The stimulus current is then negated from the sign expected by Chaste if evaluating
            dV/dt gives a positive value.
            """
            assert len(current_units) == 1 # We are using the stimulus units
            vars = []
            def f(ci_elt):
                v = ci_elt.variable
                if v.get_source_variable(recurse=True) is self.i_stim_var:
                    v.set_value(1.0)
                    found_stim.set()
                else:
                    u = v.component.get_units_by_name(v.units)
                    if u.dimensionally_equivalent(current_units[0]):
                        v.set_value(0.0)
                    elif not v.is_statically_const(ignore_annotations=True):
                        v.set_value(1.0)
                    vars.append(v)
            for expr in exprs:
                self._process_ci_elts(expr, f)
            if not found_stim:
                bfs(assign_values_for_stimulus_check, vars, found_stim=found_stim)

        # Iterate over all expressions in the model, to find the one for dV/d(something)
        for expr in (e for e in self.doc.model.get_assignments() if isinstance(e, mathml_apply) and e.is_ode()):
            # Assume the independent variable is time; if it isn't, we'll catch this later
            (dep_var, time_var) = expr.assigned_variable()
            if dep_var.get_source_variable(recurse=True) is self.V_variable:
                # Recursively search for ionic currents
                find_currents([expr.eq.rhs])
                if not ionic_vars:
                    # The sign checks below will be nonsense in this case. An error will be raised later.
                    break
                # Check the sign of the RHS
                self.i_ionic_negated = expr.eq.rhs.evaluate() > 0.0
                clear_values(expr.eq.rhs, process_definitions=True)
                if self.i_stim_var:
                    # Check the sign of the stimulus current
                    assign_values_for_stimulus_check([expr.eq.rhs])
                    self.i_stim_negated = expr.eq.rhs.evaluate() > 0.0
                    clear_values(expr.eq.rhs, process_definitions=True)
                # Found dV/d(something); don't check any more expressions
                break
        DEBUG('config', "Found ionic currents from dV/dt: ", ionic_vars)
        call_if(self.i_ionic_negated, DEBUG, 'config', "Ionic current is negated")
        call_if(self.i_stim_negated, DEBUG, 'config', "Stimulus current is negated")
        return ionic_vars
    
    def _find_var(self, oxmeta_name, definitions):
        """Find the variable object in the model for a particular concept.
        
        Will look for a variable annotated with the given oxmeta_name first, then
        try the list of definitions from the configuration file in turn.
        """
        var = None
        # Prepend an oxmeta definition
        oxmeta_defn = self._create_var_def(oxmeta_name, 'oxmeta')
        for defn in [oxmeta_defn] + definitions:
            var = self._find_variable(defn)
            if var:
                break
        return var

    def find_current_vars(self):
        """Find the variables representing currents."""
        # Find the stimulus current, if it exists for this kind of model (some are self-excitatory)
        if not self.doc.model.is_self_excitatory():
            self.i_stim_var = self._find_var('membrane_stimulus_current', self.i_stim_definitions)
            if self.i_stim_var.get_type() == VarTypes.Mapped:
                print >>sys.stderr, "Mapped variable specified as stimulus; using its source instead."
                self.i_stim_var = self.i_stim_var.get_source_variable(recurse=True)
            DEBUG('config', 'Found stimulus', self.i_stim_var)
            if not self.i_stim_var:
                # No match :(
                msg = "No stimulus current found; you'll have trouble generating Chaste code"
                if self.options.fully_automatic:
                    raise ConfigurationError(msg)
                else:
                    print >>sys.stderr, msg
                    self.i_stim_var = None
        # For other ionic currents, try using the equation for dV/dt unless told otherwise
        if not self.options.use_i_ionic_regexp:
            self.i_ionic_vars = self._find_transmembrane_currents_from_voltage_ode()
        else:
            for defn in self.i_ionic_definitions:
                if getattr(defn, u'type', u'name') != u'name':
                    raise ConfigurationError('Ionic current definitions have to have type "name"')
                regexps = unicode(defn).strip().split(',')
                comp_re = re.compile(regexps[0] + '$')
                var_re = re.compile(regexps[1] + '$')
                for component in getattr(self.doc.model, u'component', []):
                    if comp_re.match(unicode(component.name).strip()):
                        for var in getattr(component, u'variable', []):
                            if (var is not self.i_stim_var and
                                var_re.match(unicode(var.name).strip())):
                                self.i_ionic_vars.append(var)
        if not self.i_ionic_vars:
            msg = "No ionic currents found; you'll have trouble generating Chaste code"
            if self.options.fully_automatic:
                raise ConfigurationError(msg)
            else:
                print >>sys.stderr, msg
        return

    def _parse_lookup_tables(self, lookup_tables):
        """Parse a lookup_tables element."""
        for lt in lookup_tables.lookup_table:
            var_type = getattr(lt.var, u'type', u'name')
            var_name = unicode(lt.var).strip()
            config_key = (var_type, var_name)
            if not config_key in self.lut_config:
                self.lut_config[config_key] = {}
                self._set_lut_defaults(self.lut_config[config_key])
            for elt in lt.xml_element_children():
                if elt.localName != u'var':
                    self.lut_config[config_key]['table_' + elt.localName] = unicode(elt).strip()
            if hasattr(lt, u'units'):
                try:
                    units = self.unit_definitions.get_units_by_name(lt.units)
                except KeyError:
                    raise ConfigurationError('The units "%s" referenced by the lookup table for "%s" do not exist'
                                             % (lt.units, var_name))
                self.lut_config[config_key]['table_units'] = units
        return

    def _set_lut_defaults(self, lut_dict):
        """Set default configuration for a lookup table."""
        def_dict = optimize.LookupTableAnalyser._LT_DEFAULTS
        for k, v in def_dict.iteritems():
            if k != 'table_var':
                lut_dict[k] = v
        lut_dict['table_units'] = None
        return

    def annotate_currents_for_pe(self):
        """Annotate ionic & stimulus current vars so PE doesn't remove them.
        Also annotate the membrane capacitance, if defined."""
        if self.i_stim_var:
            self.i_stim_var.set_pe_keep(True)
        for var in self.i_ionic_vars:
            var.set_pe_keep(True)
        if self.Cm_variable:
            self.Cm_variable.set_pe_keep(True)
        return
    
    def expose_variables(self):
        """Expose variables for access with GetAnyVariable if desired."""
        def annotate(var):
            t = var.get_type()
            if t == VarTypes.Constant:
                var.set_is_modifiable_parameter(True)
            elif t in [VarTypes.Computed, VarTypes.Free, VarTypes.Mapped]:
                var.set_is_derived_quantity(True)
        if self.options.expose_annotated_variables:
            for var in self.metadata_vars:
                if (not self.options.use_chaste_stimulus or
                    not var.oxmeta_name in cellml_metadata.STIMULUS_NAMES):
                    annotate(var)
            DEBUG('translate', "+++ Exposed annotated variables")
        if self.options.expose_all_variables:
            for var in self.doc.model.get_all_variables():
                annotate(var)
            DEBUG('translate', "+++ Exposed all variables")
    
    def annotate_metadata_for_pe(self):
        "Annotate all vars tagged with metadata so PE doesn't remove them."
        for var in self.metadata_vars:
            var.set_pe_keep(True)
        return

    def find_transmembrane_potential(self):
        """Find and store the variable object representing V.

        Tries metadata annotation first.  If that fails, uses the name given in
        the command line options, if present.  If that fails, uses the config file.
        """
        if not self.options:
            raise ValueError('No command line options given')
        # Check command line option before config file
        if self.options.transmembrane_potential:
            self.V_definitions[0:0] = [self.options.transmembrane_potential.strip().split(',')]
            if len(self.V_definitions[0]) != 2:
                raise ConfigurationError('The name of V must contain both component and variable name')
        self.V_variable = self._find_var('membrane_voltage', self.V_definitions)
        DEBUG('config', 'Found V', self.V_variable)
        if not self.V_variable and not self.options.protocol:
            raise ConfigurationError('No transmembrane potential found; check your configuration')
        return self.V_variable
    
    def find_membrane_capacitance(self):
        """Find and store the variable object representing the cell membrane capacitance.
        
        Uses first metadata, if present, then the configuration file."""
        self.Cm_variable = self._find_var('membrane_capacitance', self.Cm_definitions)
        DEBUG('config', 'Found capacitance', self.Cm_variable)

    def find_cytosolic_calcium_concentration(self):
        """Find and store the variable object representing the cytosolic_calcium_concentration.
        
        Uses metadata only, if does not store."""
        self.cytosolic_calcium_variable = None
        self.cytosolic_calcium_variable = self.doc.model.get_variable_by_oxmeta_name('cytosolic_calcium_concentration', throw=False)
        if(self.cytosolic_calcium_variable):
            DEBUG('config', 'Found capaccytosolic_calcium_variable', self.cytosolic_calcium_variable)
        else:
            DEBUG('config', 'capaccytosolic_calcium_variable NOT found', self.cytosolic_calcium_variable)

    def find_lookup_variables(self):
        """Find the variable objects used as lookup table keys.

        This method translates the variable names given in the configuration file into objects
        in the document, and then uses those objects as keys in our lut_config dictionary.
        The ultimate source variable for the variable specified is used, in order to avoid
        complications caused by intermediaries being removed (e.g. by PE).

        The table settings are also units-converted to match the units of the key variable.
        """
        new_config = {}
        for key in self.lut_config:
            defn_type, content = key
            defn = self._create_var_def(content, defn_type)
            var = self._find_variable(defn)
            if not var:
                # Variable doesn't exist, so we can't index on it
                LOG('lookup-tables', logging.WARNING, 'Variable', content, 'not found, so not using as table index.')
            else:
                var = var.get_source_variable(recurse=True)
                if not var in new_config:
                    new_config[var] = {}
                new_config[var].update(self.lut_config[key])
                # Apply units conversions to the table settings if required
                table_units = new_config[var]['table_units']
                if table_units:
                    var_units = var.get_units()
                    if not table_units.dimensionally_equivalent(var_units):
                        LOG('lookup-tables', logging.WARNING, 'Variable', content, 'is in units', var_units.description(),
                            'which are incompatible with', table_units.description(), 'so not using as table index.')
                    elif not table_units.equals(var_units):
                        # New setting[var_units] = m[var_units/table_units]*(setting-o1[table_units]) + o2[var_units]
                        # c.f. mathml_units_mixin._add_units_conversion
                        print 'LT conversion:', table_units.description(), 'to', var_units.description(), 'equal?', table_units.equals(var_units)
                        m = table_units.get_multiplicative_factor() / var_units.get_multiplicative_factor()
                        for setting in new_config[var]:
                            try:
                                old_value = float(new_config[var][setting])
                                new_value = m * (old_value - table_units.get_offset()) + var_units.get_offset()
                                new_config[var][setting] = unicode(new_value)
                                print 'LT conversion', setting, old_value, new_value
                            except (ValueError, TypeError):
                                pass
        self.lut_config = new_config
        DEBUG('config', 'Lookup tables configuration:', new_config)
        return

    # TODO - move into seperate metadata class?
    def validate_metadata(self, assume_valid=False):
        """Check that the metadata annotations are 'sensible'.
        
        Ensures that only names we know are used, and that the same name isn't used for multiple variables.
        """
        vars = cellml_metadata.find_variables(self.doc.model, ('bqbiol:is', NSS['bqbiol']))
        self.metadata_vars = filter(lambda v: v.oxmeta_name, vars)
        if assume_valid:
            return
        names_used = [var.oxmeta_name for var in self.metadata_vars]
        DEBUG('metadata', 'Names found: ', names_used)
        # Check all metadata is allowed
        unknown_names = frozenset(names_used) - cellml_metadata.METADATA_NAMES
        if unknown_names:
            msg = ['Unrecognised oxmeta variable names found (run with --assume-valid to ignore):']
            msg.extend(sorted(unknown_names))
            raise ConfigurationError('\n  '.join(msg))
        # Check for duplicates
        d = {}
        for name in names_used:
            if name in d:
                raise ConfigurationError(name + ' metadata attribute is duplicated in the cellml file.')
            else:
                d[name] = name


######################################################################
#                    For running as an executable                    #
######################################################################

def get_options(args, default_options=None):
    """get_options(args):
    Process our command-line options.

    args is a list of options & positional arguments.
    
    default_options, if given, is an instance of optparse.Values created by a
    previous call to this function.
    """
    usage = 'usage: %prog [options] <cellml file or URI>'
    parser = optparse.OptionParser(version="%%prog %s" % __version__,
                                   usage=usage)
    parser.add_option('-q', '--quiet', action='store_true', default=False,
                      help="don't show warning messages, only errors")
    # What type of translation is being performed
    parser.add_option('-T', '--translate',
                      dest='translate', action='store_true',
                      default=True,
                      help="output computer code [default]")
    parser.add_option('-C', '--output-cellml',
                      dest='translate', action='store_false',
                      help="output an annotated CellML file instead of translating, on stdout unless -o specified")
    translators = sorted(CellMLTranslator.translators)
    parser.add_option('-t', '--translate-type',
                      type='choice', choices=translators,
                      default='Chaste', metavar='TYPE',
                      help="the type of code to output [default: %default].  "
                      "Choices: " + str(translators))
    parser.add_option('-o', dest='outfilename', metavar='OUTFILE',
                      help="write program code to OUTFILE [default action is to use the input filename with a different extension]")
    # Global adjustment settings
    parser.add_option('--config-file',
                      action='append', default=[],
                      help="pathname of configuration file")
    parser.add_option('-A', '--fully-automatic',
                      action='store_true', default=False,
                      help="if human intervention is required, fail noisily")
    parser.add_option('--assume-valid',
                      action='store_true', default=False,
                      help="skip some of the model validation checks")
    parser.add_option('--warn-on-unit-conversions',
                      action='store_true', default=False,
                      help="generate a warning if unit conversions are required")
    parser.add_option('--Wu', '--warn-on-units-errors',
                      action='store_true', default=False,
                      dest='warn_on_units_errors',
                      help="give a warning instead of an error for dimensional inconsistencies")
    parser.add_option('-V', '--transmembrane-potential', default=None, metavar='POT_VAR',
                      help="POT_VAR is the full name of the variable representing the transmembrane potential."
                      "  If not specified here, the configuration file will be used, which is the prefered method."
                      "  Defaults to 'membrane,V'.")
    parser.add_option('-d', '--debug', action='store_true', default=False,
                      help="output debug info to stderr")
    parser.add_option('-D', '--debug-source', action='append',
                      help="only show debug info from the specified part of the code."
                      "  This option may appear more than once to select multiple sources.  Implies -d.")
    parser.add_option('--profile', action='store_true', default=False,
                      help="turn on profiling of PyCml")
    # To examine the profile do something like:
    #    import os,pstats
    #    os.chdir('/tmp')
    #    files = filter(lambda f: f.startswith('pycml'), os.listdir('.'))
    #    p = pstats.Stats(*files)
    #    p.strip_dirs().sort_stats('cum').print_stats(15)
    # What optimisations/transformations to do
    group = optparse.OptionGroup(parser, 'Transformations',
                                 "These options control which transformations (typically optimisations) are applied in the generated code")
    group.add_option('-l', '--lookup-tables',
                     dest='lut', action='store_true', default=False,
                     help="perform a lookup table analysis")
    group.add_option('-p', '--pe', '--partial-evaluation',
                     dest='pe', action='store_true', default=False,
                     help="partially evaluate the model")
    group.add_option('-u', '--units-conversions',
                     action='store_true', default=False,
                     help="add explicit units conversion mathematics")
    group.add_option('-j', '--maple-output',
                     metavar='FILENAME', default=None,
                     help="file containing output from a Maple script generated using -J.  The generated"
                     " code/CellML will then contain a symbolic Jacobian as computed by Maple.")
    group.add_option('-J', '--do-jacobian-analysis',
                     action='store_true', default=False,
                     help="generate code to perform Jacobian analysis for backward Euler & CVODE; implies -t Maple")
    group.add_option('--backward-euler',
                     action='store_true', default=False,
                     help="generate a specialised cell model that solves itself using a decoupled"
                     " backward Euler method.  Not compatible with --rush-larsen.  Implies -t Chaste."
                     "  Requires -j.")
    group.add_option('--rush-larsen',
                     action='store_true', default=False,
                     help="use the Rush-Larsen method to solve Hodgkin-Huxley style gating variable"
                     " equations.  Not compatible with --backward-euler.  Implies -t Chaste.")
    group.add_option('--grl1',
                     action='store_true', default=False,
                     help="use the GRL1 method to solve Hodgkin-Huxley style gating variable"
                     " equations.  Not compatible with the backward Euler transformation."
                     " Implies -t Chaste.")
    group.add_option('--grl2',
                     action='store_true', default=False,
                     help="use the GRL2 method to solve Hodgkin-Huxley style gating variable"
                     " equations.  Not compatible with the backward Euler transformation."
                     " Implies -t Chaste.")
    parser.add_option_group(group)
    # Settings tweaking the generated code
    group = optparse.OptionGroup(parser, 'Generated code options')
    group.add_option('-c', '--class-name', default=None,
                     help="explicitly set the name of the generated class")
    group.add_option('-a', '--augment-class-name',
                     dest='augment_class_name', action='store_true',
                     default=False,
                     help="alter the class name to show what transformations are used")
    group.add_option('--no-timestamp',
                     action='store_true', default=False,
                     help="don't add a timestamp comment to generated files")
    parser.add_option_group(group)
    # Options specific to Maple output
    group = optparse.OptionGroup(parser, 'Maple options', "Options specific to Maple code output")
    group.add_option('--dont-omit-constants',
                     dest='omit_constants', action='store_false', default=True,
                     help="when generating Maple code, include assignments of constants")
    group.add_option('--compute-partial-jacobian', dest='compute_full_jacobian',
                     action='store_false', default=True,
                     help="make generated Maple code compute a Jacobian specific to a Newton solve"
                     " of the nonlinear portion of the ODE system, rather than the full system Jacobian")
    parser.add_option_group(group)
    # Options specific to Python output
    group = optparse.OptionGroup(parser, 'Python options', "Options specific to Python code output")
    group.add_option('--no-numba', dest='numba', default=True, action='store_false',
                     help="turn off using Numba to optimise code on-the-fly")
    parser.add_option_group(group)
    # Options specific to Chaste output
    group = optparse.OptionGroup(parser, 'Chaste options', "Options specific to Chaste code output")
    group.add_option('-y', '--dll', '--dynamically-loadable',
                     dest='dynamically_loadable',
                     action='store_true', default=False,
                     help="add code to allow the model to be compiled to a shared library and dynamically loaded"
                     " (only works if -t Chaste is used)")
    group.add_option('--use-chaste-stimulus',
                     action='store_true', default=False,
                     help="when generating Chaste code, use Chaste's stimulus rather than that defined in the model")
    group.add_option('--no-use-chaste-stimulus', dest='use_chaste_stimulus',
                     action='store_false',
                     help="when generating Chaste code, use the model's stimulus, not Chaste's")
    group.add_option('-i', '--convert-interfaces',
                     action='store_true', default=False,
                     help="perform units conversions at interfaces to Chaste (only works if -t Chaste is used)")
    group.add_option('--use-i-ionic-regexp', dest='use_i_ionic_regexp',
                     action='store_true', default=False,
                     help="determine ionic currents from the regexp specified in the config file"
                     " rather than analysing the voltage derivative equation")
    group.add_option('--include-dt-in-tables',
                     action='store_true', default=False,
                     help="[experimental] allow timestep to be included in lookup tables.  By default"
                     " uses the timestep of the first cell created.  Requires support from external"
                     " code if timestep changes.  Only really useful for backward Euler cells.")
    group.add_option('-m', '--use-modifiers',
                     action='store_true', default=False,
                     help="[experimental] add modifier functions for certain"
                     " metadata-annotated variables for use in sensitivity analysis (only works if -t Chaste is used)")
    group.add_option('--use-data-clamp',
                     action='store_true', default=False,
                     help="[experimental] generate a data clamp subclass of CVODE cells"
                     " which contains data clamp currents for fitting experimental data (only works if -t CVODE is used)")
    group.add_option('--expose-annotated-variables',
                     action='store_true', default=False,
                     help="expose all oxmeta-annotated variables for access via the GetAnyVariable functionality")
    group.add_option('--expose-all-variables',
                     action='store_true', default=False,
                     help="expose all variables for access via the GetAnyVariable functionality")
    parser.add_option_group(group)
    # Options specific to Functional Curation
    group = optparse.OptionGroup(parser, 'Functional Curation options', "Options specific to use by Functional Curation")
    def protocol_callback(option, opt_str, value, parser):
        """
        Protocols don't always produce normal cardiac cell models.
        However, we want to allow a later option to override these changes.
        """
        parser.values.protocol = value
        parser.values.convert_interfaces = False
        parser.values.use_chaste_stimulus = False
    group.add_option('--protocol',
                     action='callback', callback=protocol_callback, type='string', nargs=1,
                     help="specify a simulation protocol to apply to the model prior to translation")
    group.add_option('--protocol-options', action='store', type='string',
                     help="extra options for the protocol")
    group.add_option('--expose-named-parameters',
                     action='store_true', default=False,
                     help="expose all constant variables with 'name' annotations for access as model parameters")
    parser.add_option_group(group)
    # Settings for lookup tables
    group = optparse.OptionGroup(parser, 'Lookup tables options', "Options specific to the lookup tables optimisation")
    lookup_type_choices = ['entry-below', 'nearest-neighbour', 'linear-interpolation']
    group.add_option('--lookup-type', choices=lookup_type_choices,
                     default='linear-interpolation',
                     help="the type of table lookup to perform [default: %default]."
                     " Choices: " + str(lookup_type_choices))
    group.add_option('--no-separate-lut-class', dest='separate_lut_class',
                     action='store_false', default=True,
                     help="don't put lookup tables in a separate class")
    group.add_option('--row-lookup-method',
                     action='store_true', default=True,
                     help="add and use a method to look up a whole row of a table")
    group.add_option('--no-row-lookup-method', dest='row_lookup_method',
                     action='store_false',
                     help="don't add and use a method to look up a whole row of a table")
    group.add_option('--combine-commutative-tables',
                     action='store_true', default=False,
                     help="optimise a special corner case to reduce the number of tables."
                     " See documentation for details.")
    group.add_option('--lt-index-uses-floor',
                     action='store_true', default=False,
                     help="use floor() to calculate LT indices, instead of just casting")
    group.add_option('--constrain-table-indices',
                     action='store_true', default=False,
                     help="constrain lookup table index variables to remain within the bounds specified,"
                     " rather than throwing an exception if they go outside the bounds")
    group.add_option('--no-check-lt-bounds', dest='check_lt_bounds',
                     action='store_false', default=True,
                     help="[unsafe] don't check for LT indexes going outside the table bounds")
    parser.add_option_group(group)
    # Settings for partial evaluation
    group = optparse.OptionGroup(parser, 'Partial evaluation options', "Options specific to the partial evaluation optimisation")
    group.add_option('--pe-convert-power',
                     action='store_true', default=False,
                     help="convert pow(x,3) to x*x*x; similarly for powers 2 & 4.")
    group.add_option('--no-partial-pe-commutative', dest='partial_pe_commutative',
                     action='store_false', default=True,
                     help="don't combine static operands of dynamic commutative associative applys")
    group.add_option('--no-pe-instantiate-tables', dest='pe_instantiate_tables',
                     action='store_false', default=True,
                     help="don't instantiate definitions that will be tables regardless of usage")
    parser.add_option_group(group)

    options, args = parser.parse_args(args, values=default_options)
    if len(args) != 1:
        parser.error("exactly one input CellML file must be specified")

    # Some options imply others
    if options.debug_source:
        options.debug = True
    if options.do_jacobian_analysis:
        options.translate_type = 'Maple'
        options.maple_output = False
        options.rush_larsen = False
        options.backward_euler = False
    if options.backward_euler:
        if not options.maple_output:
            parser.error("Backward Euler code generation requires maple output (-j)")
        options.rush_larsen = False
        options.grl1 = False
        options.grl2 = False
    if options.rush_larsen or options.backward_euler or options.grl1 or options.grl2:
        options.translate_type = 'Chaste'
    if options.use_data_clamp and not options.translate_type=='CVODE':
        parser.error("Data clamp option '--use-data-clamp' also requires CVODE ('-t CVODE'). If you are calling this via ConvertCellModel use '--cvode-data-clamp'.")
    # Numba may not be available
    if options.numba:
        try:
            import numba
        except:
            options.numba = False 
    return options, args[0]


def load_model(model_file, options):
    """Load and validate a CellML model."""
    # Setup logging
    logging.thread = None # Hack: we're not multi-threaded, so be slightly quicker...
    if options.debug:
        formatter = logging.Formatter(fmt="%(name)s: %(message)s")
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(formatter)
        handler.addFilter(OnlyDebugFilter())
        if options.debug_source:
            handler.addFilter(OnlyTheseSourcesFilter(options.debug_source))
        logging.getLogger().addHandler(handler)
        logging.getLogger().setLevel(logging.DEBUG)

    # We can't translate if some warnings occur, as well as if the model is invalid
    notifier = NotifyHandler(level=logging.WARNING_TRANSLATE_ERROR)
    logging.getLogger('validator').addHandler(notifier)
    v = validator.CellMLValidator(create_relaxng_validator=not options.assume_valid)
    valid, doc = v.validate(model_file, return_doc=True, show_warnings=not options.quiet,
                            check_for_units_conversions=options.warn_on_unit_conversions,
                            warn_on_units_errors=options.warn_on_units_errors,
                            assume_valid=options.assume_valid)
    v.quit()
    del v

    if not valid or notifier.messages:
        print >>sys.stderr, model_file,
        if not valid:
            print >>sys.stderr, "is not a valid CellML file"
        else:
            print >>sys.stderr, "contains untranslatable constructs (see warnings above for details)"
        sys.exit(1)
    
    return doc

def run():
    """Translate the file given on the command line."""
    options, model_file = get_options(sys.argv[1:])
    doc = load_model(model_file, options)
    DEBUG('translate', "+++ Loaded model")

    config = ConfigurationStore(doc, options=options)
    for config_file in options.config_file:
        config.read_configuration_file(config_file)
    DEBUG('translate', "+++ Read config")
    
    # Apply protocol, if given
    if options.protocol:
        import protocol
        protocol.apply_protocol_file(doc, options.protocol)
        if options.debug:
            post_proto_cellml = options.outfilename or model_file
            post_proto_cellml = os.path.splitext(post_proto_cellml)[0] + '-proto.cellml.ppp'
            stream = open_output_stream(post_proto_cellml)
            doc.xml(indent=u'yes', stream=stream)
            close_output_stream(stream)
        DEBUG('translate', "+++ Applied protocol")

    config.finalize_config()
    DEBUG('translate', "+++ Processed config")

    solver_info = SolverInfo(doc.model)

    # Generate an interface component, if desired
    translator_klass = CellMLTranslator.translators[options.translate_type]
    if not options.protocol:
        translator_klass.generate_interface(doc, solver_info)
    config.validate_metadata(options.assume_valid)
    DEBUG('translate', "+++ Generated interface")
    
    if options.lut:
        config.find_lookup_variables()
        DEBUG('translate', "+++ Found LT keys")

    # These bits could do with improving, as they annotate more than is really needed!
    if options.pe:
        # We need to ensure PE doesn't remove ionic currents needed for GetIIonic
        config.annotate_currents_for_pe()
        # "Need" to ensure pe doesn't remove metadata-annotated variables (when using modifiers or default stimulus?)
        config.annotate_metadata_for_pe()
        DEBUG('translate', "+++ Annotated variables")
    # Deal with the 'expose' options
    config.expose_variables()

    class_name = options.class_name
    if not class_name:
        class_name = doc.model.name.replace('-', '_')
        if options.augment_class_name:
            class_name = u'CML_' + class_name
            if options.pe:
                class_name += '_pe'
            if options.lut:
                class_name += '_lut'
            if options.backward_euler:
                class_name += '_be'
            if options.use_modifiers:
                class_name += '_sens'
    if options.protocol:
        # Try to avoid OdeSystemInformation conflicts
        class_name += "_Proto_" + os.path.splitext(os.path.basename(options.protocol))[0]

    output_filename = getattr(options, 'outfilename', None)
    if not options.translate and not output_filename:
        output_filename = 'stdout'

    if options.units_conversions:
        doc.model.add_units_conversions()
        DEBUG('translate', "+++ Added units conversions")

    if options.do_jacobian_analysis:
        lin = optimize.LinearityAnalyser()
        lin.analyse_for_jacobian(doc, V=config.V_variable)
        DEBUG('translate', "+++ Analysed model for Jacobian")

    if options.maple_output:
        # Parse Jacobian matrix
        from maple_parser import MapleParser
        mp = MapleParser()
        jacobian_file = file(options.maple_output) # TODO: Error checking
        doc.model._cml_jacobian = mp.parse(jacobian_file)
        doc.model._cml_jacobian_full = mp.JacobianWasFullSize
        jacobian_file.close()
        if not options.backward_euler and doc.model._cml_jacobian_full:
            # Add full jacobian to XML
            solver_info.add_jacobian_matrix()
            solver_info.add_variable_links()

    if options.backward_euler:
        # Rearrange linear ODEs
        lin = optimize.LinearityAnalyser()
        lin.analyse_for_jacobian(doc, V=config.V_variable)
        lin.rearrange_linear_odes(doc)
        # Remove jacobian entries that don't correspond to nonlinear state variables
        jacobian = doc.model._cml_jacobian
        if isinstance(jacobian, tuple):
            assert doc.model._cml_jacobian_full
            jacobian = jacobian[1]
        nonlinear_vars = set([v.get_source_variable(recurse=True) for v in doc.model._cml_nonlinear_system_variables])
        def gv(vname):
            return cellml_variable.get_variable_object(doc.model, vname).get_source_variable(recurse=True)
        for var_i, var_j in jacobian.keys():
            #hack to not remove converted variables.
            if gv(var_i).get_type()==VarTypes.Computed:
                source_var_i = gv(var_i).get_dependencies()[0].get_dependencies()[0].get_source_variable()
            else:
                source_var_i = gv(var_i)
            if gv(var_j).get_type()==VarTypes.Computed:
                source_var_j = gv(var_j).get_dependencies()[0].get_dependencies()[0].get_source_variable()
            else:
                source_var_j = gv(var_j)

            if source_var_i not in nonlinear_vars or source_var_j not in nonlinear_vars:
                del jacobian[(var_i, var_j)]
        if doc.model._cml_jacobian_full:
            # Transform the Jacobian into the form needed by the Backward Euler code
            import maple_parser
            for key, expr in jacobian.iteritems():
                new_expr = None
                if key[0] == key[1]:
                    # 1 on the diagonal
                    new_expr = maple_parser.MNumber(['1'])
                if not (isinstance(expr, maple_parser.MNumber) and str(expr) == '0'):                
                    # subtract delta_t * expr
                    args = []
                    if new_expr:
                        args.append(new_expr)
                    args.append(maple_parser.MOperator([maple_parser.MVariable(['delta_t']), expr], 'prod', 'times'))
                    new_expr = maple_parser.MOperator(args, '', 'minus')
#                    jacobian[key] = new_expr
                if new_expr:
                    jacobian[key] = new_expr                    

        # Add info as XML
        solver_info.add_all_info()
        # Analyse the XML, adding cellml_variable references, etc.
        solver_info.add_variable_links()
        solver_info.add_linear_ode_update_equations()
        DEBUG('translate', "+++ Parsed and incorporated Maple output")
    else:
        options.include_dt_in_tables = False

    if options.lut:
        # Create the analyser so PE knows which variables are table keys
        lut = optimize.LookupTableAnalyser()
    else:
        lut = None

    if options.pe:
        # Do partial evaluation
        pe = optimize.PartialEvaluator()
        pe.parteval(doc, solver_info, lut)
        DEBUG('translate', "+++ Done PE")

    if options.lut:
        # Do the lookup table analysis
        lut.analyse_model(doc, solver_info)
        DEBUG('translate', "+++ Done LT analysis")
    
    if options.rush_larsen:
        rl = optimize.RushLarsenAnalyser()
        rl.analyse_model(doc)
        DEBUG('translate', "+++ Done Rush-Larsen analysis")

    if options.translate:
        # Translate to code
        initargs = {'add_timestamp': not options.no_timestamp,
                    'options': options}
        transargs = {'v_variable': config.V_variable}
        transargs['row_lookup_method'] = options.row_lookup_method
        transargs['lt_index_uses_floor'] = options.lt_index_uses_floor
        transargs['constrain_table_indices'] = options.constrain_table_indices
        if issubclass(translator_klass, CellMLToMapleTranslator):
            initargs['omit_constants'] = options.omit_constants
            initargs['compute_full_jacobian'] = options.compute_full_jacobian
        elif issubclass(translator_klass, CellMLToChasteTranslator):
            solver_info.add_membrane_ionic_current()
            transargs['use_chaste_stimulus'] = options.use_chaste_stimulus
            transargs['separate_lut_class'] = options.separate_lut_class
            transargs['convert_interfaces'] = options.convert_interfaces
            transargs['use_modifiers'] = options.use_modifiers
            transargs['use_data_clamp'] = options.use_data_clamp
            transargs['dynamically_loadable'] = options.dynamically_loadable
            transargs['use_protocol'] = bool(options.protocol)
        t = translator_klass(**initargs)
        t.translate(doc, model_file, output_filename, class_name=class_name, **transargs)
        cellml_metadata.remove_model(doc.model)
    else:
        # Add a comment element
        comment = pycml.comment_base(
            body=u'\n' + version_comment(not options.no_timestamp) + u'\n')
        doc.xml_insert_before(doc.model, comment)
        # Output annotated model
        stream = open_output_stream(output_filename)
        doc.xml(indent=u'yes', stream=stream)
        close_output_stream(stream)
        
    DEBUG('translate', "+++ Done translation")

