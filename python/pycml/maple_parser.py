
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
Parser for Maple output.

We use PyParsing to parse the output from Maple, and get a Python
datatype representing the elements of the Jacobian matrix.

The output from Maple consists of a series of lines.  The initial
lines are the input expressions, and so should be ignored.  The
interesting part comes once we see a line of the form:
"--<variable_i>/<variable_j>--"
This line is then followed by the (i,j)-th entry of the Jacobian,
possibly occuring on multiple lines.  There will be one such pair
of 'lines' for each entry of the matrix.

Identifying these pairs is done by iterating over the lines, and
checking the first character to identify the interesting
portions---only the marker lines will start with a double quote.
PyParsing can then be used to process the Jacobian entries.
"""

import copy
import logging
import re
import sys

from pyparsing import *


__version__ = "$Revision$"
__all__ = ['MapleParser']




# Class hierarchy for expressions we can use
class MExpression(object):
    """A mathematical expression."""
    _cache = {}
    _cache_uses = {}
    _temporaries = {}
    _op_and_func_names = set()
    def uniquify(self):
        """Ensure only one copy of this expression exists.

        If an equal expression is in the cache, return that.
        Otherwise, return self and add us to the cache.

        This must only be called after the expression has been normalized.

        Operates recursively to uniquify children before checking the cache for this expression.
        """
        if hasattr(self, '_children'):
            # Uniquify children
            orig_children = self._children[:]
            self._children[:] = []
            for child in orig_children:
                self._children.append(child.uniquify())
        if self in self._cache:
            # This expression is used in multiple locations.
            self._cache_uses[self] += 1
            # If it is complex, generated code should assign it to a temporary variable
            if self.is_complex():
                self._temporaries[self] = self._cache[self]
            return self._cache[self]
        else:
            self._cache[self] = self
            self._cache_uses[self] = 1
            return self
    
    @classmethod
    def filter_temporaries(cls, exprs, _temps={}, _threshold=1):
        for expr in exprs:
            if expr in cls._temporaries:
                if not expr in _temps and cls._cache_uses[expr] >= _threshold:
                    _temps[expr] = cls._cache_uses[expr]
                new_threshold = cls._cache_uses[expr] + 1
            else:
                new_threshold = cls._cache_uses[expr]
            if hasattr(expr, '_children'):
                cls.filter_temporaries(expr._children, _temps, new_threshold)
        return _temps

    @classmethod
    def clear_cache(cls):
        """Clear the expression caches."""
        cls._cache.clear()
        cls._cache_uses.clear()
        cls._temporaries.clear()
    
    def __repr__(self):
        return str(self)
    def __ne__(self, other):
        """Test for non-equality of expressions.

        Equality is only defined in subclasses, but non-equality
        is always the negation of equality.
        """
        return not self.__eq__(other)
    
    def normalize(self):
        """Normalize the expression."""
        if hasattr(self, '_children'):
            orig_children = self._children[:]
            self._children[:] = []
            for child in orig_children:
                self._children.append(child.normalize())
        return self

    def xml(self):
        """Return a serialised XML representation of this expression.

        Encapulates the expression in a math element with namespace
        declarations for CellML and MathML.
        """
        xml = ['<math xmlns="http://www.w3.org/1998/Math/MathML"',
               '      xmlns:cellml="http://www.cellml.org/cellml/1.0#">',
               self.mathml(), '</math>']
        return ''.join(xml)

class MNumber(MExpression):
    """A number."""
    __slots__ = ['_value']
    def __init__(self, toklist):
        self._value = toklist[0]
        if self._value[0] == '.':
            self._value = '0' + self._value
    def __str__(self):
        return self._value
    def __eq__(self, other):
        return isinstance(other, MNumber) and (self._value == other._value)
    def __hash__(self):
        return hash((self.__class__.__name__, self._value))
    
    def is_complex(self):
        return False

    def mathml(self):
        """Return a serialised MathML representation of this expression.

        Units have to be fudged, unfortunately; they will be given as dimensionless.
        """
        return '<cn cellml:units="dimensionless">' + self._value + '</cn>'

class MVariable(MExpression):
    """A variable."""
    __slots__ = ['_name']
    def __init__(self, toklist):
        self._name = toklist[0]
    def __str__(self):
        return self._name
    def __eq__(self, other):
        return isinstance(other, MVariable) and (self._name == other._name)
    def __hash__(self):
        return hash((self.__class__.__name__, self._name))

    def is_complex(self):
        return False

    def mathml(self):
        """Return a serialised MathML representation of this expression."""
        if self._name == 'Pi':
            return '<pi/>'
        else:
            return '<ci>' + self._name + '</ci>'

class MPiecewise(MExpression):
    """A piecewise expression."""
    def __init__(self, s, loc, toklist):
        self._children = []
        for piece in toklist.piecewise.pieces:
            self._children.extend(piece)
    def __str__(self):
        return 'piecewise<' + ', '.join(map(str, self._children)) + '>'
    def __eq__(self, other):
        return isinstance(other, MPiecewise) and (self._children == other._children)
    def __hash__(self):
        k = [self.__class__.__name__] + self._children
        return hash(tuple(k))
    
    def is_complex(self):
        return True
    
    _otherwise = MVariable(['otherwise'])
    def _is_otherwise(self, cond):
        """Determine if the given condition is an 'otherwise'."""
        return cond == self._otherwise
    
    def normalize(self):
        """Make this expression nicely tree-structured.
        
        Maple often doesn't supply a derivative at the discontinuity of a piecewise expression.
        Code generation can't cope with this, so we arbitrarily pick the value from one side.
        If one side has a result that's just a variable or number, we use that, otherwise we
        use the first side that can have the condition changed to include the boundary point.
        """
        # Normalize our children
        super(MPiecewise, self).normalize()
        # Get a list of pieces
        pieces = []
        for i in range(len(self._children)/2):
            result = self._children[2*i]
            cond = self._children[2*i+1]
            pieces.append([result, cond])
        #print "Piecewise normalize:", len(self._children), len(pieces)
        #print "Before", self
        # Helper functions
        undef = MVariable(['undefined'])
        def is_undef(result):
            """Test if a result is 'Float(undefined)' or 'undefined'."""
            return ((isinstance(result, MFunction) and result._name == 'Float' and result._args[0] == undef)
                    or result == undef)
        def update_piece(piece):
            """Change a piece so the condition includes equality."""
            cond = piece[1]
            if self._is_otherwise(cond):
                return
            if isinstance(cond, MVariable):
                cond = copy.copy(MSequence.definitions[cond._name])
                piece[1] = cond
            assert isinstance(cond, MOperator)
            if cond._operator == 'lt':
                cond._operator = 'leq'
            elif cond._operator == 'gt':
                cond._operator = 'geq'
            elif cond._operator in ['leq', 'geq']:
                pass
            elif is_undef(piece[0]):
                pass # Two undefineds next to each other - don't bother changing the second!
            else:
                raise ValueError('Not a suitable operator!\t' + str(cond))
        def is_nice(piece):
            """Check if this is a piece that should definitely be modified,."""
            result, cond = piece
            return isinstance(result, (MNumber, MVariable)) or self._is_otherwise(cond)
        def is_updateable(piece):
            """Check if we can alter the condition in this piece."""
            _, cond = piece
            return isinstance(cond, MOperator) and cond._operator in ['lt', 'gt']
        # Check for boundary pieces
        num_pieces = len(pieces)
        indices_to_delete = []
        for i, piece in enumerate(pieces):
            result, cond = piece
            if is_undef(result):
                # Delete this piece, and adjust one of its neighbours
                indices_to_delete.append(i)
                neighbours = []
                if i > 0:
                    neighbours.append(pieces[i-1])
                if i < num_pieces - 1:
                    neighbours.append(pieces[i+1])
                for n in neighbours:
                    if is_nice(n):
                        update_piece(n)
                        break
                else:
                    for n in neighbours:
                        if is_updateable(n):
                            update_piece(n)
                            break
                    else:
                        update_piece(neighbours[0])
        for i in reversed(indices_to_delete):
            del pieces[i]
        # Re-create our children from the remaining pieces
        self._children = []
        for piece in pieces:
            self._children.extend(piece)
        #print "After", self
        return self
        
    def mathml(self):
        """Return a serialised MathML representation of this expression."""
        xml = ['<piecewise>']
        for i in range(len(self._children)/2):
            result = self._children[2*i]
            cond = self._children[2*i+1]
            if self._is_otherwise(cond):
                xml.append('<otherwise>%s</otherwise>' % result.mathml())
            else:
                xml.append('<piece>%s%s</piece>' % (result.mathml(), cond.mathml()))
        xml.append('</piecewise>')
        return ''.join(xml)

class MFunction(MExpression):
    """A function call."""
    NAME_MAP = {'log10': 'log'}
    def __init__(self, toklist):
        name = toklist.function.fn_name
        self._name = self.NAME_MAP.get(name, name)
        self._op_and_func_names.add(self._name)
        self._args = list(toklist.function.fn_args)
        self._children = self._args
    def __str__(self):
        return self._name + '<' + ','.join(map(str, self._args)) + '>'
    def __eq__(self, other):
        return isinstance(other, MFunction) and (self._name == other._name) and \
               (self._args == other._args)
    def __hash__(self):
        k = [self.__class__.__name__, self._name]
        k.extend(self._args)
        return hash(tuple(k))

    def is_complex(self):
        return True

    def mathml(self):
        """Return a serialised MathML representation of this expression."""
        if self._name in ['`if`', 'piecewise']:
            # Generate a piecewise expression
            xml = ['<piecewise><piece>']
            xml.append(self._args[1].mathml()) # 'then' case
            xml.append(self._args[0].mathml()) # condition
            xml.append('</piece><otherwise>')
            xml.append(self._args[2].mathml()) # 'else' case
            xml.append('</otherwise></piecewise>')
        else:
            xml = ['<apply><', self._name, '/>']
            for arg in self._args:
                xml.append(arg.mathml())
            xml.append('</apply>')
        return ''.join(xml)

class MOperator(MExpression):
    """An application of an operator."""
    OPERATOR_MAP = {'-': 'minus', '+': 'plus', '*': 'times', '/': 'divide',
                    '!': 'factorial', '^': 'power',
                    '=': 'eq', '<>': 'neq', '<=': 'leq', '<': 'lt',
                    '>=': 'geq', '>': 'gt'}
    def __init__(self, toklist, op_type, operator=None):
        self._op_type = op_type
        if isinstance(toklist, ParseResults):
            self._toks = toklist[0].asList()
        if operator:
            self._operator = operator
            self._op_and_func_names.add(operator)
            self._operands = self._children = toklist
        else:
            self._operator = '#'
            self._operands = self._children = []
    def __str__(self):
        return self._operator + '<' + ','.join(map(str, self._operands)) + '>'
    def __eq__(self, other):
        return isinstance(other, MOperator) and \
               (self._operator == other._operator) and \
               (self._op_type == other._op_type) and \
               (self._operands == other._operands) and \
               (not (hasattr(self, '_toks') and hasattr(other, '_toks')) or
                (self._toks == other._toks))
    def __hash__(self):
        k = [self.__class__.__name__, self._operator]
        k.extend(self._operands)
        return hash(tuple(k))

    def is_complex(self):
        result = False
        nested_children = False
        for child in self._children:
            if child.is_complex():
                result = True
                break
            if hasattr(child, '_children'):
                nested_children = True
        if not result and (self._operator not in ['minus', 'plus', 'times', 'eq', 'neq', 'lt', 'leq', 'gt', 'geq']
                           or len(self._children) > 2 or nested_children):
            result = True
        return result

    def normalize(self):
        """Make this expression properly tree-structured."""
        if self._operands:
            # We've been normalized already
            return self
        def op_name(op):
            try:
                return self.OPERATOR_MAP[op]
            except KeyError:
                return op
        if self._op_type == 'sign':
            if self._toks[0] == '+':
                # Leading plus is a no-op
                return self._toks[1].normalize()
            else:
                # Unary minus
                self._operator = 'minus'
                self._op_and_func_names.add(self._operator)
                self._children = self._operands = [self._toks[1].normalize()]
                return self
        elif self._op_type in ['prod', 'plus', 'rel']:
            # Convert to a left-factored binary tree
            # e.g. a-b-c+d+e = (((a-b)-c)+d)+e
            res = self._toks[0].normalize()
            for i in range((len(self._toks)-1)/2):
                operator = op_name(self._toks[1+i*2])
                operand = self._toks[2+i*2].normalize()
                res = MOperator([res, operand], self._op_type, operator)
            return res
        else:
            # Just normalize operands and set operator name
            self._operator = op_name(self._toks[1])
            self._op_and_func_names.add(self._operator)
            self._children = self._operands = map(lambda op: op.normalize(), self._toks[::2])
            return self

    def mathml(self):
        """Return a serialised MathML representation of this expression."""
        xml = ['<apply><', self._operator, '/>']
        for arg in self._operands:
            xml.append(arg.mathml())
        xml.extend(['</apply>'])
        return ''.join(xml)

class MDerivative(MExpression):
    """A derivative expression."""
    __slots__ = ['_dependent_var', '_independent_var']
    def __init__(self, dependent_var, independent_var):
        self._dependent_var = dependent_var
        self._independent_var = independent_var
    def __str__(self):
        return 'diff<' + self._dependent_var + ',' + self._independent_var + '>'
    def __eq__(self, other):
        return isinstance(other, MDerivative) and \
               (self._independent_var == other._independent_var) and \
               (self._dependent_var == other._dependent_var)
    def __hash__(self):
        return hash((self.__class__.__name__,
                     self._independent_var, self._dependent_var))

    def is_complex(self):
        return False

    def mathml(self):
        """Return a serialised MathML representation of this expression."""
        xml = ['<apply><diff/><bvar>']
        xml.append(self._independent_variable.mathml())
        xml.append('</bvar>')
        xml.append(self._dependent_variable.mathml())
        xml.append('</apply>')
        return ''.join(xml)

class MSequence(MExpression):
    """A computation sequence."""
    def __init__(self, toklist):
        names, children = [], []
        for assignment in toklist[0]:
            names.append(assignment[0])
            children.append(assignment[1])
        self._names = names
        self._children = children

    def __eq__(self, other):
        return isinstance(other, MSequence) and (self._names == other._names) and (self._children == other._children)
    def __hash__(self):
        return hash(tuple(self._names) + tuple(self._children))
    def __str__(self):
        return '\n'.join([str(n) + ' = ' + str(self._children[i]) for i, n in enumerate(self._names)])
    
    def is_complex(self):
        return True
    
    definitions = {}
    def normalize(self):
        """Normalize the expression, and record assignments for easy looking up."""
        new_children = []
        for i, child in enumerate(self._children):
            new_child = child.normalize()
            new_children.append(new_child)
            self.definitions[self._names[i]] = new_child
        self._children = new_children
        return self
    
    def mathml(self):
        """Return a serialised MathML representation of the temporary variable definitions."""
        xml = []
        for i, name in enumerate(self._names):
            if isinstance(name, basestring):
                value = self._children[i]
                xml.append('<apply><eq/><ci>%s</ci>%s</apply>' % (name, value.mathml()))
        return ''.join(xml)
    
    def get_jacobian_entries(self):
        """Return a dictionary mapping (i,j) indices to Jacobian matrix entries."""
        J = {}
        for i, name in enumerate(self._names):
            if not isinstance(name, basestring):
                value = self._children[i]
                index = (int(name[1]), int(name[2]))
                J[index] = value
        return J


def make_moperator(op_type):
    def pa(s, l, t):
        #print "Creating an", op_type
        return MOperator(t, op_type)
    return pa


# Necessary for reasonable speed
ParserElement.enablePackrat()


# Punctuation
oparen = Literal('(').suppress()
cparen = Literal(')').suppress()
comma = Literal(',').suppress()

# Expressions must be constructed recursively
expr = Forward().setName('Expression')
#wrapped_expr = oparen + expr + cparen

# Identifiers start with a letter or underscore.
# Subsequent chars may also be digits.
initchars = alphas + '_'
ident = Word(initchars, initchars+nums)

# Variables
var = ident.copy().setName('Variable').setParseAction(MVariable)
# Function calls
fname = (ident | "`if`").setResultsName('fn_name')
fname.setName('Func_name')
fargs = Group(expr + ZeroOrMore(comma + expr)).setResultsName('fn_args')
fargs.setName('Func_args')
func = Group(fname + oparen + fargs + cparen).setResultsName('function')
func.setName('Function')
func.setParseAction(MFunction)

# Piecewise
osquare = Literal('[').suppress()
csquare = Literal(']').suppress()
piece = Group(osquare + expr + comma + expr + csquare).setResultsName('piece')
piece.setName('Piece')
pieces = Group(piece + ZeroOrMore(comma + piece)).setResultsName('pieces')
pieces.setName('Pieces')
piecewise = Group('PIECEWISE' + oparen + pieces + cparen).setResultsName('piecewise')
piecewise.setName('Piecewise')
piecewise.setParseAction(MPiecewise)

# Numbers can be given in scientific notation, with an optional sign.
real_re = Regex(r'([0-9]+\.?([0-9]+(e[-+]?[0-9]+)?)?)|(\.[0-9]+(e[-+]?[0-9]+)?)')
real_re.setName('Number')
real_re.setParseAction(MNumber)

# Operators.
# All are left-associative, apart from:
#  not - right-assoc
#  ^   - non-assoc
#  *, +, and, or - really n-ary, but binary is ok
op_fact = Literal('!') # Check if use 2! or !2 (expects 2! at present)
op_expt = Literal('^')
#op_sign = oneOf('+ -')
op_sign = Literal('-')
op_prod = oneOf('* /')
op_plus = oneOf('+ -')
op_rel = oneOf('= <> <= >= < >')
op_not = Literal('not')
op_and = Literal('and')
op_or = Literal('or')
op_xor = Literal('xor')

for op in ['fact', 'expt', 'sign', 'prod', 'plus', 'rel', 'not', 'and', 'or', 'xor']:
    d = {'op': 'op_' + op}
    exec "%(op)s = %(op)s.setName('%(op)s').setResultsName('%(op)s')" % d

# Base terms, in order (hopefully) of speed of determining non-match
atom = real_re | piecewise | func | var

# The main expression grammar
expr << operatorPrecedence(atom, [
    (op_fact, 1, opAssoc.LEFT, make_moperator('fact')),
    (op_expt, 2, opAssoc.LEFT, make_moperator('expt')), # Fudge associativity
    (op_sign, 1, opAssoc.RIGHT, make_moperator('sign')),
    (op_prod, 2, opAssoc.LEFT, make_moperator('prod')),
    (op_plus, 2, opAssoc.LEFT, make_moperator('plus')),
    (op_rel, 2, opAssoc.LEFT, make_moperator('rel')),
    (op_not, 1, opAssoc.RIGHT, make_moperator('not')),
    (op_and, 2, opAssoc.LEFT, make_moperator('and')),
    (op_or, 2, opAssoc.LEFT, make_moperator('or')),
    (op_xor, 2, opAssoc.LEFT, make_moperator('xor'))])

# We can also see "computation sequences" of name=expr
equals = Literal('=').suppress()
#extended_ident = Regex(r'[a-zA-Z_][a-zA-Z0-9_]+( ?\[ ?[0-9]+ ?, ?[0-9]+ ?\])?')
extended_ident = Group(Literal('jacobian') + osquare + Word(nums) + comma + Word(nums) + csquare) | ident
assign = Group(extended_ident + equals + expr).setResultsName('assign')
computation_seq = Group(assign + ZeroOrMore(comma + assign)).setParseAction(MSequence)




class MapleParser(object):
    """
    A parser for results output by Maple.
    
    The results consist of a series of lines.  The initial lines are the
    input expressions, and so should be ignored.  The interesting part comes
    once we see a line of the form:
    
    "--<variable_i>/<variable_j>--"

    This line is then followed by the (i,j)-th entry of the Jacobian,
    possibly occuring on multiple lines.
    There will be one such pair of 'lines' for each entry of the matrix.
    
    Lines may be explicitly continued, if they terminate with a backslash,
    or implicitly, if context shows that the expression is incomplete.
    
    The class has one method, parse, which parses a file-like object.
    """
    
    def __init__(self):
        """Initialise the parser."""
        # We just store the original stack limit here, so we can increase
        # it for the lifetime of this object if needed for parsing, on the
        # basis that if one expression needs to, several are likely to.
        self._original_stack_limit = sys.getrecursionlimit()
        self._stack_limit_multiplier = 1
        
    def __del__(self):
        """Reset the stack limit if it changed."""
        sys.setrecursionlimit(self._original_stack_limit)
    
    def parse_full_jacobian(self, stream):
        """Parse some Maple output containing a Jacobian for the full system.
        
        Input is a file-like object supporting the iterator protocol.
        
        Output is a pair (temporaries, J), where
          temporaries is a MathML string containing a sequence of apply elements,
          which represent definitions of common temporary variables
          J is a dictionary, keyed by variable name ordered pairs; the entries
          are the parsed Jacobian matrix entries.
        """
        state_var_names = {}
        var_name_re = re.compile(r'"--([0-9]+)--(.+)--')
        J = '' # The Maple code for the Jacobian
        in_J = False
        for line in stream:
            if line[0] == '"':
                # Tells us the index of a state variable: --%d--%d--
                while line.rstrip()[-1] == '\\':
                    # Explicit continuation
                    line = line.rstrip() + stream.next()
                m = var_name_re.match(line)
                state_var_names[int(m.group(1))] = m.group(2)
                continue
            if line.startswith('bytes used') or line.startswith('memory used'):
                # Maple footer, which sometimes appears in the middle...
                continue
            if line.startswith('J := '):
                assert not in_J
                in_J = True
                line = line[5:] # Strip 'J := '
            if in_J:
                J += line.strip()
                if J[-1] == '\\':
                    # Explicit continuation: remove the backslash
                    J = J[:-1]
                else:
                    # Potentially implicit continuation: add whitespace to replace the '\n'
                    J += ' '
        #print J
        result = self._parse_sequence(J)
        #print result
        temporaries = result.xml()
        jacobian = result.get_jacobian_entries()
        # Change the indices from numbers to variable names
        for i, j in jacobian.keys():
            jacobian[(state_var_names[i], state_var_names[j])] = jacobian[(i,j)]
            del jacobian[(i,j)]
        #print jacobian
        return (temporaries, jacobian)
    
    def parse(self, stream, debug=False):
        """Parse some Maple output.
        
        Input is a file-like object supporting the iterator protocol.
        
        Output is a dictionary, keyed by variable name ordered pairs;
        the entries are the parsed Jacobian matrix entries.
        
        If debug is given as True, then also return a similar dictionary
        containing the expressions as strings.
        """
        results, debug_res = {}, {}
        curr_key = None
        s = "" # Current 'line' contents
        in_header = False
        self.JacobianWasFullSize = False
        for line in stream:
            if line.strip() == '"FULL JACOBIAN"':
                self.JacobianWasFullSize = True
                return self.parse_full_jacobian(stream)
            if line.startswith('bytes used') or line.startswith('memory used'):
                # Maple footer
                if curr_key and s:
                    self._parse_expr(curr_key, s, results, debug_res)
                    curr_key = None
                    s = ""
                # The 'footer' sometimes appears in the middle...
                continue
            if line[0] == '"':
                # Start of a header
                if s and not in_header:
                    self._parse_expr(curr_key, s, results, debug_res)
                    s = ""
                in_header = True
                curr_key = None
            if curr_key:
                # This line is part of an expression...
                s = s + line.strip()
                if s[-1] == '\\':
                    # Explicit continuation: remove the backslash
                    s = s[:-1]
                else:
                    # Potentially implicit continuation: add whitespace to replace the '\n'
                    s += ' '
            if in_header:
                # This line is part of an entry header
                s = s + line.strip()
                # Is it continued? (Must be explicit continuation in this case)
                if s[-1] == '\\':
                    # remove backslash
                    s = s[:-1]
                else:
                    # Extract the variable names
                    pos = s.find('/')
                    if pos == -1:
                        raise ValueError("Bad header in Maple output: " + s)
                    var_i = s[3:pos]
                    var_j = s[pos+1:-3]
                    curr_key = (var_i, var_j)
                    in_header = False
                    s = ""
#        if 1:
#            temps = MExpression.filter_temporaries(results.values())
#            print "Filtered temporaries", len(temps)
#            for t in temps:
#                print t, MExpression._cache_uses[t]
#        print "DEBUG: Operators and functions in Maple:", sorted(MExpression._op_and_func_names)
        if debug:
            return results, debug_res
        else:
            return results
    
    def _parse_sequence(self, seq_str):
        """Parse a computation sequence, and return a list of (name, value) pairs."""
        r = self._try_parse(seq_str, computation_seq)
        n = r[0].normalize()
        u = n.uniquify()
        return u

    def _parse_expr(self, key, expr_str, results, debug_res):
        """Parse a single expression, and store the result under key."""
        self._debug("Parsing derivative", key)
        r = self._try_parse(expr_str, expr)
        n = r[0].normalize()
        u = n.uniquify()
        debug_res[key] = (expr_str, r, n)
        results[key] = u
        return

    def _try_parse(self, string, grammar_rule):
        """Try parsing a string with the given rule, repeating with a higher recursion limit if needed."""
        r = None
        while self._stack_limit_multiplier < 3: # Magic number!
            try:
                r = grammar_rule.parseString(string, parseAll=True)
            except RuntimeError, msg:
                self._debug("Got RuntimeError:", msg)
                self._stack_limit_multiplier += 0.5
                new_limit = int(self._stack_limit_multiplier * self._original_stack_limit)
                self._debug("Increasing recursion limit to", new_limit)
                sys.setrecursionlimit(new_limit)
            else:
                break # Parsed OK
        if not r:
            raise RuntimeError("Failed to parse expression even with a recursion limit of %d; giving up!"
                               % (int(self._stack_limit_multiplier * self._original_stack_limit),))
        return r

    def set_debug(self, debug=True):
        """Turn debugging on or off."""
        for v in globals().itervalues():
            if isinstance(v, ParserElement):
                v.setDebug(debug)
        return

    def _debug(self, *args):
        """Log a debug message."""
        logger = logging.getLogger('maple-parser')
        logger.debug(' '.join(map(str, args)))

    def add_assignments(self, results):
        """Turn a results dictionary into assignment expressions.

        For each entry k, v in results, change the value into an
        assignment expression, assigning v to the appropriate derivative.

        Returns a new dictionary, with the same keys but different values.
        """
        res = {}
        for k, v in results.iteritems():
            var_i, var_j = map(lambda varname: MVariable([varname]), k)
            deriv = MDerivative(var_i, var_j)
            ass_expr = MOperator([deriv, v], 'rel', 'eq')
            res[k] = ass_expr
        return res
        


if __name__ == '__main__':
    # Parse files given on the command line
    mp = MapleParser()

    def prettyprint(results):
        ks = results.keys()
        ks.sort()
        for k in ks:
            if results[k] != MNumber(['0']):
                print k, ':', results[k]
    
    for fname in sys.argv[1:]:
        res, dbg = mp.parse(file(fname), True)
        prettyprint(res)






#def mathml_to_mexpr(elt):
#    """Convert a MathML element tree to an MExpression."""
#    result = None
#    if isinstance(elt, mathml_ci):
#        var = elt.variable
#        result = MVariable([var.fullname(cellml=True)])
#    elif isinstance(elt, mathml_cn):
#        result = MNumber([unicode(elt)])
#    elif isinstance(elt, mathml_apply):
#        # Note: will have to treat diff seperately.  Others can be MFunction?
#        pass
#    elif isinstance(elt, mathml_piecewise):
#        # Make into MFunction('`if`') ?
#        pass
#    else:
#        # Hrm.
#        pass
