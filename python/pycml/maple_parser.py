
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

import logging
from pyparsing import *


__version__ = "$Revision$"
__all__ = ['MapleParser']




# Class hierarchy for expressions we can use
class MExpression(object):
    """A mathematical expression."""
    _cache = {}
    _cache_uses = {}
    _temporaries = {}
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

class MFunction(MExpression):
    """A function call."""
    def __init__(self, toklist):
        self._name = toklist.function.fn_name
        self._args = list(toklist.function.fn_args) # Check if asList() needed
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
        if self._name == '`if`':
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
        complex = False
        nested_children = False
        for child in self._children:
            if child.is_complex():
                complex = True
                break
            if hasattr(child, '_children'):
                nested_children = True
        if not complex and (self._operator not in ['minus', 'plus', 'times', 'eq', 'neq', 'lt', 'leq', 'gt', 'geq']
                            or len(self._children) > 2 or nested_children):
            complex = True
        return complex

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

def dummy(cls):
    def pa(s, l, t):
        print "Creating a", cls.__name__
        cls(t)
    return pa

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

# Numbers can be given in scientific notation, with an optional sign.
real_re = Regex(r'([0-9]+(\.[0-9]+(e[-+][0-9]+)?)?)|(\.[0-9]+(e[-+][0-9]+)?)')
real_re.setName('Number')
real_re.setParseAction(MNumber)
uint = Word(nums)
sci_e='e' + Optional(oneOf('- +')) + uint
sci_dec='.' + uint + Optional(sci_e)
real = Combine( (uint + Optional(sci_dec)) | sci_dec )
real.setName('Number')

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

# Base terms
atom = func | var | real_re
#atom = var | real


class OpParseResults(object):
    def __init__(self, toks, op_type):
        self.toks = toks[:]
        self.op_type = op_type
    def __str__(self):
        return 'OP' + self.op_type + '<' + str(self.toks) + '>'
    def __repr__(self):
        return str(self)
def make_op_action(op_type):
    def parse_action(st, loc, toks):
        pass
        #return OpParseResults(toks, op_type)
    return parse_action


_prec = operatorPrecedence(atom, [
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

expr << _prec






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
        for line in stream:
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
                    # remove the backslash
                    s = s[:-1]
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
        if debug:
            return results, debug_res
        else:
            return results

    def _parse_expr(self, key, expr_str, results, debug_res):
        """Parse a single expression, and store the result under key."""
        self._debug("Parsing derivative", key)
        old_limit = sys.getrecursionlimit()
        try:
            r = expr.parseString(expr_str)
        except RuntimeError, msg:
            self._debug("Got RuntimeError:", msg)
            new_limit = int(old_limit * 1.5)
            self._debug("Failed to parse with recursion limit of %d; increasing to %d" % (old_limit, new_limit))
            sys.setrecursionlimit(new_limit)
            r = expr.parseString(expr_str)
            sys.setrecursionlimit(old_limit)
        n = r[0].normalize()
        u = n.uniquify()
        debug_res[key] = (expr_str, r, n)
        results[key] = u
        return

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
