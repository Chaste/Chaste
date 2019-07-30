#!/usr/bin/env python

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
This part of PyCml applies various optimising transformations to CellML
models, in particular partial evaluation and the use of lookup tables.
"""

import operator

# Common CellML processing stuff
import pycml
from pycml import *  # Put contents in the local namespace as well

__version__ = "$Revision$"[11:-2]



######################################################################
#                         Partial Evaluation                         #
######################################################################

class PartialEvaluator(object):
    """Perform partial evaluation of a CellML model."""
    def _debug(self, *args):
        """Output debug info from the PE process."""
        logger = logging.getLogger('partial-evaluator')
        logger.debug(' '.join(map(str, args)))

    def _expr_lhs(self, expr):
        """Display the LHS of this expression."""
        lhs = expr.assigned_variable()
        if isinstance(lhs, cellml_variable):
            return lhs.fullname()
        else:
            return lhs[0].fullname() + u'/' + lhs[1].fullname()
        
    def _describe_expr(self, expr):
        """Describe this expression for debug info."""
        if isinstance(expr, mathml_apply):
            if expr.is_assignment() or expr.is_ode():
                return self._expr_lhs(expr)
            else:
                return '[nested apply]'
        elif isinstance(expr, mathml_ci):
            return expr.variable.fullname()
        elif isinstance(expr, mathml_cn):
            return u'cn[' + unicode(expr) + u']'
        else:
            return '[unknown]'
        
    def _process_ci_elts(self, elt, func):
        """Apply func to all ci elements in the tree rooted at elt."""
        if isinstance(elt, mathml_ci):
            func(elt)
        else:
            for e in elt.xml_element_children():
                self._process_ci_elts(e, func)
    
    def _rename_var(self, elt):
        """Change this ci element to use a canonical name."""
        if elt.xml_parent.localName == u'bvar':
            # The free variable in a derivative must refer directly to the ultimate source,
            # since this is assumed in later stages and in code generation.
            elt._set_variable_obj(elt.variable.get_source_variable(recurse=True))
        elt._rename()
        self._debug("Using canonical name", unicode(elt))

    def _do_reduce_eval_loop(self, expr_source):
        """Do the reduce/evaluate loop.
        
        expr_source is a callable that returns an iterable over expressions.
        """
        while True:
            self.doc.model._pe_repeat = u'no'
            for expr in list(expr_source()):
                self._reduce_evaluate_expression(expr)
            if self.doc.model._pe_repeat == u'no':
                break
            self._debug("----- looping -----")
        del self.doc.model._pe_repeat
    
    def _reduce_evaluate_expression(self, expr):
        """Reduce or evaluate a single expression."""
        if hasattr(expr, '_pe_process'):
            # This expression has been reduced or evaluated already, but needs further
            # processing later so hasn't been removed yet.
            return
        if expr._get_binding_time() is BINDING_TIMES.static:
            # Evaluate
            try:
                value = expr.evaluate() # Needed here for the is_assignment case
            except:
                print "Error evaluating", self._describe_expr(expr)
                raise
            self._debug("Evaluated", self._describe_expr(expr), "to", value)
            if isinstance(expr, mathml_apply):
                if expr.is_ode():
                    # Replace the RHS with a <cn> element giving the value
                    rhs = expr.eq.rhs
                    new_elt = expr._eval_self()
                    expr.xml_insert_after(rhs, new_elt)
                    expr.xml_remove_child(rhs)
                elif expr.is_assignment():
                    # The variable assigned to will have its initial_value set,
                    # so we don't need the expression any more.  Flag it for removal.
                    expr._pe_process = u'remove'
                else:
                    # Replace the expression with a <cn> element giving the value
                    new_elt = expr._eval_self()
                    expr.xml_parent.xml_insert_after(expr, new_elt)
                    expr.xml_parent.xml_remove_child(expr)
            else:
                # Replace the expression with a <cn> element giving the value
                expr._reduce()
            # Update variable usage counts for the top-level apply case
            if isinstance(expr, mathml_apply):
                if expr.is_ode() or expr.is_assignment():
                    expr._update_usage_counts(expr.eq.rhs, remove=True)
                else:
                    expr._update_usage_counts(expr, remove=True)
        else:
            # Reduce
            expr._reduce()
    
    def _get_assignment_exprs(self, skip_solver_info=True):
        """Get an iterable over all assignments in the model that are mathml_apply instances."""
        if not skip_solver_info:
            skip = set()
        else:
            skip = set(self.solver_info.get_modifiable_mathematics())
        for e in self.doc.model.get_assignments():
            if isinstance(e, mathml_apply) and e not in skip:
                assert e.is_ode() or e.is_assignment()
                yield e

    def is_instantiable(self, expr):
        """Determine whether special conditions mean that this assignment can be instantiated.
        
        Normally an assignment can only be instantiated if the assigned-to variable is used only
        once, in order to avoid code duplication.
        However, if the definition under consideration for instantiation is a function only of a
        single LT keying variable (and we will do LT) then code duplication doesn't really matter,
        since the whole expression will be converted to a table anyway.  So we should instantiate
        regardless of multiple uses in this case.
        
        Note: this does check that only a single keying variable appears, but doesn't check for
        the presence of expensive functions.  Of course, if there aren't any expensive functions,
        the code duplication isn't that worrying.
        """
        instantiate = False
        if self.lookup_tables_analyser and self.doc.model.get_option('pe_instantiate_tables'):
            keying_vars = set()
            all_keying = [True]
            def func(ci_elt):
                if self.lookup_tables_analyser.is_keying_var(ci_elt.variable):
                    keying_vars.add(ci_elt.variable)
                else:
                    all_keying[0] = False
            self._process_ci_elts(expr.eq.rhs, func)
            instantiate = len(keying_vars) == 1 and all_keying[0]
        return instantiate
    
    def _is_source_of(self, v1, v2):
        """Test if v1 is a source of the mapped variable v2."""
        if v1 is v2:
            return True
        elif v2.get_type() == VarTypes.Mapped:
            return self._is_source_of(v1, v2.get_source_variable())
        else:
            return False
    
    def _check_retargetting(self, ci_elt):
        """Check if this new variable reference means a retarget needs to change.
        
        A retarget occurs when a kept dynamic mapped variable is changed to computed
        because its source variable(s) are only used once and are not kept.  But new
        mathematics may use one or more of those sources, in which case we need to
        revert the retarget.
        """
        # Is this a retarget?
        var = ci_elt.variable
        root_defn = var.get_source_variable(recurse=True).get_dependencies()
        if root_defn:
            root_defn = root_defn[0]
        else:
            return
        if (isinstance(root_defn, mathml_apply)
             and hasattr(root_defn, '_pe_process')
             and root_defn._pe_process == 'retarget'):
            # Are we a source of the 'new' assignee?
            assignee = root_defn._cml_assigns_to
            if not self._is_source_of(assignee, var):
                if var.get_type() == VarTypes.Computed:
                    # We were the original source variable; stop retargetting
                    self._debug('Ceasing re-target of', root_defn, 'to', assignee)
                    del root_defn._pe_process
                    root_defn._cml_assigns_to = var
                else:
                    # Re-target to var instead
                    self._debug('Changing re-target of', root_defn, 'from', assignee, 'to', var)
                    var._cml_var_type = VarTypes.Computed
                    root_defn._cml_assigns_to = var
                    var._cml_depends_on = [root_defn]
                    var._cml_source_var = None
                # assignee should now map to var
                assignee._cml_var_type = VarTypes.Mapped
                assignee._cml_source_var = var
                assignee._cml_depends_on = [var]
    
    def parteval(self, doc, solver_info, lookup_tables_analyser=None):
        """Do the partial evaluation."""
        self.doc = doc
        self.solver_info = solver_info
        self.lookup_tables_analyser = lookup_tables_analyser
        if lookup_tables_analyser:
            lookup_tables_analyser.doc = doc
        doc.partial_evaluator = self
        # Do BTA and reduce/eval of main model
        doc.model.do_binding_time_analysis()
        self._do_reduce_eval_loop(self._get_assignment_exprs)
        
        if solver_info.has_modifiable_mathematics():
            # Do BTA and reduce/eval of solver info section
            for expr in solver_info.get_modifiable_mathematics():
                if not (isinstance(expr, mathml_apply) and expr.is_top_level()):
                    self._process_ci_elts(expr, lambda ci: ci.variable._used())
                self._process_ci_elts(expr, self._check_retargetting)
            solver_info.do_binding_time_analysis()
            self._do_reduce_eval_loop(solver_info.get_modifiable_mathematics)
        
        # Process flagged expressions
        for expr in list(self._get_assignment_exprs()):
            if hasattr(expr, '_pe_process'):
                if expr._pe_process == u'remove':
                    if (expr._get_binding_time() == BINDING_TIMES.dynamic and
                        isinstance(expr._cml_assigns_to, cellml_variable) and
                        expr._cml_assigns_to.get_usage_count() > 1):
                        self._debug("Keeping", repr(expr), "due to SolverInfo")
                        continue
                    expr.xml_parent.xml_remove_child(expr)
                    self.doc.model._remove_assignment(expr)
                elif expr._pe_process == u'retarget':
                    lhs = expr.eq.lhs
                    var = expr._cml_assigns_to
                    ci = mathml_ci.create_new(lhs, var.fullname(cellml=True))
                    self._debug('Re-targetting', lhs, var, ci)
                    ci._set_variable_obj(var)
                    lhs.xml_parent.xml_insert_after(lhs, ci)
                    lhs.xml_parent.xml_remove_child(lhs)
        
        # Use canonical variable names in all ci elements
        for expr in list(self._get_assignment_exprs(False)):
            # If the assigned-to variable isn't used or kept, remove the assignment
            if isinstance(expr.eq.lhs, mathml_ci):
                var = expr.eq.lhs.variable
                if not (var.get_usage_count() or var.pe_keep):
                    expr.xml_parent.xml_remove_child(expr)
                    doc.model._remove_assignment(expr)
                    continue
            self._process_ci_elts(expr, self._rename_var)
        for expr in solver_info.get_modifiable_mathematics():
            self._process_ci_elts(expr, self._rename_var)
        solver_info.use_canonical_variable_names()

        # Tidy up kept variables, in case they aren't referenced in an eq'n.
        for var in doc.model.get_all_variables():
            if var.pe_keep:
                var._reduce()
        
        # Remove unused variables
        for var in list(doc.model.get_all_variables()):
            assert var.get_usage_count() >= 0
            if var.get_usage_count() == 0 and not var.pe_keep:
                var.xml_parent._del_variable(var)

        # Collapse into a single component
        new_comp = cellml_component.create_new(doc, u'c')
        new_comp._cml_created_by_pe = True
        old_comps = list(getattr(doc.model, u'component', []))
        doc.model._add_component(new_comp)
        # We iterate over a copy of the component list so we can delete components
        # from the model in this loop, and so the new component exists in the model
        # so we can add content to it.
        for comp in old_comps:
            # Move relevant contents into new_comp
            for units in list(getattr(comp, u'units', [])):
                # Copy all <units> elements
                # TODO: Just generate the ones we need, using _ensure_units_exist
                comp.xml_remove_child(units)
                new_comp.xml_append(units)
            for var in list(getattr(comp, u'variable', [])):
                # Only move used source variables
                self._debug('Variable', var.fullname(), 'usage', var.get_usage_count(),
                           'type', var.get_type(), 'kept', var.pe_keep)
                if (var.get_usage_count() and var.get_type() != VarTypes.Mapped) or var.pe_keep:
                    self._debug('Moving variable', var.fullname(cellml=True))
                    # Remove from where it was
                    comp._del_variable(var, keep_annotations=True)
                    # Set name to canonical version
                    var.name = var.fullname(cellml=True)
                    # Place in new component
                    new_comp._add_variable(var)
            # Don't copy reactions
            for math in list(getattr(comp, u'math', [])):
                # Copy all <math> elements with content
                if math.xml_children:
                    comp.xml_remove_child(math)
                    new_comp.xml_append(math)
                    # Invalidate cached links
                    math._unset_cached_links()
            doc.model._del_component(comp)
        # Remove groups & connections
        for group in list(getattr(doc.model, u'group', [])):
            doc.model.xml_remove_child(group)
        for conn in list(getattr(doc.model, u'connection', [])):
            doc.model.xml_remove_child(conn)

        # Remove unused variable assignments from the list
        vs = [v for v in doc.model.get_assignments() if isinstance(v, cellml_variable)]
        for v in vs:
            if not v.xml_parent is new_comp:
                doc.model._remove_assignment(v)

        # Remove interface attributes from variables
        for v in new_comp.variable:
            for iface in [u'public', u'private']:
                try:
                    delattr(v, iface+u'_interface')
                except AttributeError:
                    pass

        # Refresh expression dependency lists
        for expr in self._get_assignment_exprs(False):
            expr._cml_depends_on = list(expr.vars_in(expr.eq.rhs))
            if expr.is_ode():
                # Add dependency on the independent variable
                indep_var = expr.eq.lhs.diff.independent_variable
                if not indep_var in expr._cml_depends_on:
                    expr._cml_depends_on.append(indep_var)
                # Update ODE definition dependency if needed
                expr.eq.lhs.diff.dependent_variable._update_ode_dependency(indep_var, expr)
        return


######################################################################
#                        Lookup table analysis                       #
######################################################################

class LookupTableAnalyser(object):
    """
    Analyses & annotates a CellML model to indicate where lookup
    tables can be used.
    """

    def __init__(self):
        """Create an analyser."""
        # No model to analyse yet
        self.doc = None
        # Set default parameter values
        self.set_params()

    @property
    def config(self):
        """Get the current document's configuration store."""
        return getattr(self.doc, '_cml_config', None)

    def var_is_membrane_potential(self, var):
        """Determine if the given variable represents the transmembrane potential.

        This method takes an instance of cellml_variable and returns a boolean.
        """
        return (var.name in [u'V', u'membrane__V'] and
                var.get_type(follow_maps=True) == VarTypes.State)

    def is_allowed_variable(self, var):
        """Return True iff the given variable is allowed in a lookup table.

        This method uses the config store in the document to check the variable object.
        """
        var = var.get_source_variable(recurse=True)
        allowed = (var in self.config.lut_config or
                   (self.config.options.include_dt_in_tables and
                    var is self.solver_info.get_dt().get_source_variable(recurse=True)))
        return allowed

    def is_keying_var(self, var):
        """Return True iff the given variable can be used as a table key.

        Will check the config store if it exists.  If not, the variable name must match self.table_var.
        """
        if self.config:
            return var.get_source_variable(recurse=True) in self.config.lut_config
        else:
            return var.name == self.table_var

    _LT_DEFAULTS = {'table_min': u'-100.0001',
                    'table_max': u'49.9999',
                    'table_step': u'0.01',
                    'table_var': u'V'}
    def set_params(self, **kw):
        """Set parameters controlling lookup table generation.

        Keyword parameters control the lookup table settings, which are
        stored as attributes on suitable expressions.
        table_min - minimum table entry (unicode) -> lut:min
        table_max - maximum table entry (unicode) -> lut:max
        table_step - table step size (unicode) -> lut:step
        table_var - the name of the variable indexing the table (unicode) -> lut:var
        """
        defaults = self._LT_DEFAULTS
        for attr in defaults:
            if attr in kw:
                setattr(self, attr, kw[attr])
            else:
                setattr(self, attr, getattr(self, attr, defaults[attr]))
        return

    def get_param(self, param_name, table_var):
        """Get the value of the lookup table parameter.

        table_var is the variable object being used to key this table.

        If the document has a config store, lookup the value there.
        If that doesn't give us a value, use that given using set_params.
        """
        try:
            val = self.config.lut_config[
                table_var.get_source_variable(recurse=True)][param_name]
        except AttributeError, KeyError:
            val = getattr(self, param_name)
        return val

    # One of these functions is required for a lookup table to be worthwhile
    lut_expensive_funcs = frozenset(('exp', 'log', 'ln', 'root',
                                     'sin', 'cos', 'tan',
                                     'sec', 'csc', 'cot',
                                     'sinh', 'cosh', 'tanh',
                                     'sech', 'csch', 'coth',
                                     'arcsin', 'arccos', 'arctan',
                                     'arcsinh', 'arccosh', 'arctanh',
                                     'arcsec', 'arccsc', 'arccot',
                                     'arcsech', 'arccsch', 'arccoth'))

    class LUTState(object):
        """Represents the state for lookup table analysis."""
        def __init__(self):
            """Set the initial state.

            We assume at first a lookup table would not be suitable.
            """
            self.has_var = False
            self.bad_vars = set()
            self.has_func = False
            self.table_var = None

        def update(self, res):
            """Update the state with the results of a recursive call.

            res is the result of checking a sub-expression for suitability,
            and should be another instance of this class.
            """
            self.has_var = (self.has_var or res.has_var) and \
                           (not (self.table_var and res.table_var) or
                            self.table_var.get_source_variable(recurse=True) is
                             res.table_var.get_source_variable(recurse=True))
            # The second condition above specifies that the keying variables must be the same if they both exist
            self.bad_vars.update(res.bad_vars)
            self.has_func = self.has_func or res.has_func
            self.table_var = self.table_var or res.table_var
            if not self.has_var and self.table_var:
                # Two sub-expressions have different keying variables, so consider them as bad variables
                self.bad_vars.add(self.table_var.name)
                self.bad_vars.add(res.table_var.name)
                self.table_var = None

        def suitable(self):
            """Return True iff this state indicates a suitable expression for replacement with a lookup table."""
            return (self.has_var and
                    not self.bad_vars and
                    self.has_func)

        def reason(self):
            """
            Return a unicode string describing why this state indicates the
            expression is not suitable for replacement with a lookup table.

            This can be:
             'no_var' - doesn't contain the table variable
             'bad_var <vname>' - contains a variable which isn't permitted
             'no_func' - doesn't contain an expensive function
            or a comma separated combination of the above.
            """
            r = []
            if not self.has_var:
                r.append(u'no_var')
            if not self.has_func:
                r.append(u'no_func')
            for vname in self.bad_vars:
                r.append(u'bad_var ' + vname)
            return u','.join(r)

    def create_state_from_annotations(self, expr):
        """Create a LUTState instance from an already annotated expression."""
        state = self.LUTState()
        possible = expr.getAttributeNS(NSS['lut'], u'possible', '')
        if possible == u'yes':
            varname = expr.getAttributeNS(NSS['lut'], u'var')
            state.table_var = expr.component.get_variable_by_name(varname)
        elif possible == u'no':
            reason = expr.getAttributeNS(NSS['lut'], u'reason', '')
            reasons = reason.split(u',')
            for reason in reasons:
                if reason == u'no_var':
                    state.has_var = False
                elif reason == u'no_func':
                    state.has_func = False
                elif reason.startswith(u'bad_var '):
                    state.bad_vars.add(reason[8:])
        return state

    def analyse_for_lut(self, expr, var_checker_fn):
        """Check if the given expression can be replaced by a lookup table.

        The first argument is the expression to check; the second is a
        function which takes a variable object and returns True iff this
        variable is permitted within a lookup table expression.

        If self.annotate_failures is True then annotate <apply> and
        <piecewise> expressions which don't qualify with the reason
        why they do not.
        This can be:
         'no_var' - doesn't contain the table variable
         'bad_var <vname>' - contains a variable which isn't permitted
         'no_func' - doesn't contain an expensive function
        or a comma separated combination of the above.
        The annotation is stored as the lut:reason attribute.

        If self.annotate_outermost_only is True then only annotate the
        outermost qualifying expression, rather than also annotating
        qualifying subexpressions.
        """
        # If this is a cloned expression, then just copy any annotations on the original.
        if isinstance(expr, mathml):
            source_expr = expr.get_original_of_clone()
        else:
            source_expr = None
        if source_expr and source_expr.getAttributeNS(NSS['lut'], u'possible', '') != '':
            LookupTableAnalyser.copy_lut_annotations(source_expr, expr)
            state = self.create_state_from_annotations(source_expr)
            DEBUG('lookup-tables', "No need to analyse clone", expr.xml(), state.suitable(), state.reason())
        else:
            # Initialise the indicators
            state = self.LUTState()
            # Process current node
            if isinstance(expr, mathml_ci):
                # Variable reference
                if var_checker_fn(expr.variable):
                    # Could be a permitted var that isn't a keying var
                    if self.is_keying_var(expr.variable):
                        assert state.table_var is None # Sanity check
                        state.has_var = True
                        state.table_var = expr.variable
                else:
                    state.bad_vars.add(expr.variable.name)
            elif isinstance(expr, mathml_piecewise):
                # Recurse into pieces & otherwise options
                if hasattr(expr, u'otherwise'):
                    r = self.analyse_for_lut(child_i(expr.otherwise, 1),
                                             var_checker_fn)
                    state.update(r)
                for piece in getattr(expr, u'piece', []):
                    r = self.analyse_for_lut(child_i(piece, 1), var_checker_fn)
                    state.update(r)
                    r = self.analyse_for_lut(child_i(piece, 2), var_checker_fn)
                    state.update(r)
            elif isinstance(expr, mathml_apply):
                # Check function
                if (not state.has_func and
                    expr.operator().localName in self.lut_expensive_funcs):
                    state.has_func = True
                # Check operands
                operand_states = {}
                for operand in expr.operands():
                    r = self.analyse_for_lut(operand, var_checker_fn)
                    state.update(r)
                    operand_states[id(operand)] = r
                # Check qualifiers
                for qual in expr.qualifiers():
                    r = self.analyse_for_lut(qual, var_checker_fn)
                    state.update(r)
                # Special case additional optimisations
                if self.config and self.config.options.combine_commutative_tables and not state.suitable():
                    if isinstance(expr.operator(), reduce_commutative_nary):
                        self.check_commutative_tables(expr, operand_states)
                    elif isinstance(expr.operator(), mathml_divide):
                        self.check_divide_by_table(expr, operand_states)
            else:
                # Just recurse into children
                for e in expr.xml_children:
                    if getattr(e, 'nodeType', None) == Node.ELEMENT_NODE:
                        r = self.analyse_for_lut(e, var_checker_fn)
                        state.update(r)
        # Annotate the expression if appropriate
        if isinstance(expr, (mathml_apply, mathml_piecewise)):
            if state.suitable():
                self.annotate_as_suitable(expr, state.table_var)
            else:
                if self.annotate_failures:
                    expr.xml_set_attribute((u'lut:reason', NSS['lut']), state.reason())
        return state
    
    def check_divide_by_table(self, expr, operand_states):
        """Convert division by a table into multiplication.
        
        This is called if expr, a division, cannot be replaced as a whole by a lookup table.
        If the denominator can be replaced by a table, then convert expr into a multiplication
        by the reciprocal, moving the division into the table.
        """
        numer, denom = list(expr.operands())
        state = operand_states[id(denom)]
        if state.suitable():
            expr.safe_remove_child(numer)
            expr.safe_remove_child(denom)
            recip = mathml_apply.create_new(expr, u'divide', [(u'1', u'dimensionless'), denom])
            times = mathml_apply.create_new(expr, u'times', [numer, recip])
            expr.replace_child(expr, times, expr.xml_parent)
            self.annotate_as_suitable(recip, state.table_var)

    def check_commutative_tables(self, expr, operand_states):
        """Check whether we can combine suitable operands into a new expression.
        
        If expr has a commutative (and associative) n-ary operator, but is not suitable as a
        whole to become a lookup table (checked by caller) then we might still be able to
        do slightly better than just analysing its operands.  If multiple operands can be
        replaced by tables keyed on the same variable, these can be combined into a new
        application of the same operator as expr, which can then be replaced as a whole
        by a single lookup table, and made an operand of expr.
        
        Alternatively, if at least one operand can be replaced by a table, and a subset of
        other operands do not contain other variables, then they can be included in the single
        table.
        """
        # Operands that can be replaced by tables
        table_operands = filter(lambda op: operand_states[id(op)].suitable(), expr.operands())
        if not table_operands:
            return
        # Sort these suitable operands by table_var (map var id to var & operand list, respectively)
        table_vars, table_var_operands = {}, {}
        for oper in table_operands:
            table_var = operand_states[id(oper)].table_var
            table_var_id = id(table_var)
            if not table_var_id in table_vars:
                table_vars[table_var_id] = table_var
                table_var_operands[table_var_id] = []
            table_var_operands[table_var_id].append(oper)
        # Figure out if any operands aren't suitable by themselves but could be included in a table
        potential_operands = {id(None): []}
        for table_var_id in table_vars.keys():
            potential_operands[table_var_id] = []
        for op in expr.operands():
            state = operand_states[id(op)]
            if not state.suitable() and not state.bad_vars:
                if not state.table_var in potential_operands:
                    potential_operands[id(state.table_var)] = []
                potential_operands[id(state.table_var)].append(op)
        # Do any combining
        for table_var_id in table_vars.keys():
            suitable_opers = table_var_operands[table_var_id] + potential_operands[table_var_id] + potential_operands[id(None)]
            if len(suitable_opers) > 1:
                # Create new sub-expression with the suitable operands
                for oper in suitable_opers:
                    expr.safe_remove_child(oper)
                new_expr = mathml_apply.create_new(expr, expr.operator().localName, suitable_opers)
                expr.xml_append(new_expr)
                self.annotate_as_suitable(new_expr, table_vars[table_var_id])
                # Remove the operands with no table_var from consideration with other keying vars
                potential_operands[id(None)] = []
    
    def annotate_as_suitable(self, expr, table_var):
        """Annotate the given expression as being suitable for a lookup table."""
        if self.annotate_outermost_only:
            # Remove annotations from (expr and) child expressions
            self.remove_lut_annotations(expr)
        for param in ['min', 'max', 'step']:
            expr.xml_set_attribute((u'lut:' + param, NSS['lut']),
                                   self.get_param('table_' + param, table_var))
        expr.xml_set_attribute((u'lut:var', NSS['lut']), table_var.name)
        expr.xml_set_attribute((u'lut:possible', NSS['lut']), u'yes')
        self.doc.lookup_tables[expr] = True
    
    @staticmethod
    def copy_lut_annotations(from_expr, to_expr):
        """Copy any lookup table annotations from one expression to another."""
        for pyname, fullname in from_expr.xml_attributes.iteritems():
            if fullname[1] == NSS['lut']:
                to_expr.xml_set_attribute(fullname, getattr(from_expr, pyname))

    def remove_lut_annotations(self, expr, remove_reason=False):
        """Remove lookup table annotations from the given expression.

        By default this will only remove annotations from expressions
        (and sub-expressions) that can be converted to use lookup tables.
        If remove_reason is True, then the lut:reason attributes on
        non-qualifying expressions will also be removed.
        """
        # Remove from this expression
        delete_table = False
        for pyname in getattr(expr, 'xml_attributes', {}).keys():
            fullname = expr.xml_attributes[pyname]
            if fullname[1] == NSS['lut']:
                if remove_reason or fullname[0] != u'lut:reason':
                    expr.__delattr__(pyname)
                    if fullname[0] != u'lut:reason':
                        delete_table = True
        # Delete expr from list of lookup tables?
        if delete_table:
            del self.doc.lookup_tables[expr]
        # Recurse into children
        for e in expr.xml_children:
            if getattr(e, 'nodeType', None) == Node.ELEMENT_NODE:
                self.remove_lut_annotations(e, remove_reason)

    def analyse_model(self, doc, solver_info,
                      annotate_failures=True,
                      annotate_outermost_only=True):
        """Analyse the given document.

        This method checks all expressions (and subexpressions)
        in the given document for whether they could be converted to
        use a lookup table, and annotates them appropriately.

        By default expressions which don't qualify will be annotated
        to indicate why; set annotate_failures to False to suppress
        this.

        Also by default only the outermost suitable expression in any
        given tree will be annotated; if you want to annotate suitable
        subexpressions of a suitable expression then pass
        annotate_outermost_only as False.
        """
        self.doc = doc
        self.solver_info = solver_info
        self.annotate_failures = annotate_failures
        self.annotate_outermost_only = annotate_outermost_only
        doc.lookup_tables = {}
        # How to check for allowed variables
        if hasattr(doc, '_cml_config'):
            checker_fn = self.is_allowed_variable
        else:
            checker_fn = self.var_is_membrane_potential

        # Check all expressions
        for expr in (e for e in doc.model.get_assignments()
                     if isinstance(e, mathml_apply)):
            ops = expr.operands()
            ops.next()
            e = ops.next()
            self.analyse_for_lut(e, checker_fn)
        for expr in solver_info.get_modifiable_mathematics():
            self.analyse_for_lut(expr, checker_fn)

        if solver_info.has_modifiable_mathematics():
            self._determine_unneeded_tables()

        # Assign names (numbers) to the lookup tables found.
        # Also work out which ones can share index variables into the table.
        doc.lookup_tables = doc.lookup_tables.keys()
        doc.lookup_tables.sort(cmp=element_path_cmp)
        doc.lookup_table_indexes, n = {}, 0
        for i, expr in enumerate(doc.lookup_tables):
            expr.xml_set_attribute((u'lut:table_name', NSS['lut']), unicode(i))
            comp = expr.get_component()
            var = comp.get_variable_by_name(expr.var).get_source_variable(recurse=True)
            key = (expr.min, expr.max, expr.step, var)
            if not key in doc.lookup_table_indexes:
                doc.lookup_table_indexes[key] = unicode(n)
                n += 1
            expr.xml_set_attribute((u'lut:table_index', NSS['lut']), doc.lookup_table_indexes[key])
        
        if solver_info.has_modifiable_mathematics():
            self._determine_duplicate_tables()
            
        # Re-do dependency analysis so that an expression using lookup
        # tables only depends on the keying variable.
        for expr in (e for e in doc.model.get_assignments()
                     if isinstance(e, mathml_apply)):
            expr.classify_variables(root=True,
                                    dependencies_only=True,
                                    needs_special_treatment=self.calculate_dependencies)

    def _find_tables(self, expr, table_dict):
        """Helper method for _determine_unneeded_tables."""
        if expr.getAttributeNS(NSS['lut'], u'possible', '') == u'yes':
            table_dict[id(expr)] = expr
        else:
            for e in self.doc.model.xml_element_children(expr):
                self._find_tables(e, table_dict)

    def _determine_unneeded_tables(self):
        """Determine whether some expressions identified as lookup tables aren't actually used.
        
        This occurs if some ODEs have been linearised, in which case the original definitions
        will have been analysed for lookup tables, but aren't actually used.
        
        TODO: The original definitions might be used for computing derived quantities...
        """
        original_tables = {}
        new_tables = {}
        def f(exprs, table_dict):
            exprs = filter(lambda n: isinstance(n, (mathml_ci, mathml_apply, mathml_piecewise)), exprs)
            for node in self.doc.model.calculate_extended_dependencies(exprs):
                if isinstance(node, mathml_apply):
                    self._find_tables(node, table_dict)
        for u, t, eqns in self.solver_info.get_linearised_odes():
            original_defn = u.get_ode_dependency(t)
            f([original_defn], original_tables)
            f(eqns, new_tables)
        for id_ in set(original_tables.keys()) - set(new_tables.keys()):
            expr = original_tables[id_]
            self.remove_lut_annotations(expr)
            expr.xml_set_attribute((u'lut:reason', NSS['lut']),
                                   u'Expression will not be used in generated code.')
            DEBUG('lookup-tables', 'Not annotating probably unused expression', expr)
    
    def _determine_duplicate_tables(self):
        """Determine whether we have multiple tables for the same expression.
        
        Any expression that is identical to a previous table will be re-annotated to refer to the
        previous table, instead of declaring a new one.
        
        This is a temporary measure until we have proper sub-expression elimination for the Jacobian
        and residual calculations.
        """
        uniq_tables = []
        for expr in self.doc.lookup_tables:
            for table in uniq_tables:
                if expr.same_tree(table):
                    lt_name = table.getAttributeNS(NSS['lut'], u'table_name', u'')
                    # Need to remove old name before we can set a new one (grr amara)
                    del expr.table_name
                    expr.xml_set_attribute((u'lut:table_name', NSS['lut']), lt_name)
                    break
            else:
                uniq_tables.append(expr)

    def calculate_dependencies(self, expr):
        """Determine the dependencies of an expression that might use a lookup table.

        This method is suitable for use as the needs_special_treatment function in
        mathml_apply.classify_variables.  It is used to override the default recursion
        into sub-trees.  It takes a single sub-tree as argument, and returns either
        the dependency set for that sub-tree, or None to use the default recursion.
        
        Expressions that can use a lookup table only depend on the keying variable.
        """
        if expr.getAttributeNS(NSS['lut'], u'possible', '') == u'yes':
            key_var_name = expr.getAttributeNS(NSS['lut'], u'var')
            key_var = expr.component.get_variable_by_name(key_var_name).get_source_variable(recurse=True)
            return set([key_var])
        # If not a table, use default behaviour
        return None


######################################################################
#                          Jacobian analysis                         #
######################################################################

class LinearityAnalyser(object):
    """Analyse linearity aspects of a model.

    This class performs analyses to determine which ODEs have a linear
    dependence on their state variable, discounting the presence of
    the transmembrane potential.

    This can be used to decouple the ODEs for more efficient solution,
    especially in a multi-cellular context.

    analyse_for_jacobian(doc) must be called before
    rearrange_linear_odes(doc).
    """
    LINEAR_KINDS = Enum('None', 'Linear', 'Nonlinear')

    def analyse_for_jacobian(self, doc, V=None):
        """Analyse the model for computing a symbolic Jacobian.

        Determines automatically which variables will need to be solved
        for using Newton's method, and stores their names in
        doc.model._cml_nonlinear_system_variables, as a list of variable
        objects.

        Also stores doc.model._cml_linear_vars, doc.model._cml_free_var,
        doc.model._cml_transmembrane_potential.
        """
        # TODO: Add error checking and tidy
        stvs = doc.model.find_state_vars()
        if V is None:
            Vcname, Vvname = 'membrane', 'V'
            V = doc.model.get_variable_by_name(Vcname, Vvname)
        V = V.get_source_variable(recurse=True)
        doc.model._cml_transmembrane_potential = V
        free_var = doc.model.find_free_vars()[0]
        lvs = self.find_linear_odes(stvs, V, free_var)
        # Next 3 lines for benefit of rearrange_linear_odes(doc)
        lvs.sort(key=lambda v: v.fullname())
        doc.model._cml_linear_vars = lvs
        doc.model._cml_free_var = free_var
        # Store nonlinear vars in canonical order
        nonlinear_vars = list(set(stvs) - set([V]) - set(lvs))
        nonlinear_vars.sort(key=lambda v: v.fullname())
        doc.model._cml_nonlinear_system_variables = nonlinear_vars
        # Debugging
        f = lambda var: var.fullname()
        DEBUG('linearity-analyser', 'V=', V.fullname(), '; free var=',
              free_var.fullname(), '; linear vars=', map(f, lvs),
              '; nonlinear vars=', map(f, nonlinear_vars))
        return

    def _get_rhs(self, expr):
        """Return the RHS of an assignment expression."""
        ops = expr.operands()
        ops.next()
        return ops.next()

    def _check_expr(self, expr, state_var, bad_vars):
        """The actual linearity checking function.

        Recursively determine the type of dependence expr has on
        state_var.  The presence of any members of bad_vars indicates
        a non-linear dependency.

        Return a member of the self.LINEAR_KINDS enum.
        """
        kind = self.LINEAR_KINDS
        result = None
        if isinstance(expr, mathml_ci):
            var = expr.variable.get_source_variable(recurse=True)
            if var is state_var:
                result = kind.Linear
            elif var in bad_vars:
                result = kind.Nonlinear
            elif var.get_type(follow_maps=True) == VarTypes.Computed:
                # Recurse into defining expression
                src_var = var.get_source_variable(recurse=True)
                src_expr = self._get_rhs(src_var.get_dependencies()[0])
                DEBUG('find-linear-deps', "--recurse for", src_var.name,
                      "to", src_expr)
                result = self._check_expr(src_expr, state_var, bad_vars)
            else:
                result = kind.None
            # Record the kind of this variable, for later use when
            # rearranging linear ODEs
            var._cml_linear_kind = result
        elif isinstance(expr, mathml_cn):
            result = kind.None
        elif isinstance(expr, mathml_piecewise):
            # If any conditions have a dependence, then we're
            # nonlinear.  Otherwise, all the pieces must be the same
            # (and that's what we are) or we're nonlinear.
            pieces = getattr(expr, u'piece', [])
            conds = map(lambda p: child_i(p, 2), pieces)
            chld_exprs = map(lambda p: child_i(p, 1), pieces)
            if hasattr(expr, u'otherwise'):
                chld_exprs.append(child_i(expr.otherwise, 1))
            for cond in conds:
                if self._check_expr(cond, state_var, bad_vars) != kind.None:
                    result = kind.Nonlinear
                    break
            else:
                # Conditions all OK
                for e in chld_exprs:
                    res = self._check_expr(e, state_var, bad_vars)
                    if result is not None and res != result:
                        # We have a difference
                        result = kind.Nonlinear
                        break
                    result = res
        elif isinstance(expr, mathml_apply):
            # Behaviour depends on the operator
            operator = expr.operator().localName
            operands = expr.operands()
            if operator in ['plus', 'minus']:
                # Linear if any operand linear, and none non-linear
                op_kinds = map(lambda op: self._check_expr(op, state_var,
                                                           bad_vars),
                               operands)
                result = max(op_kinds)
            elif operator == 'divide':
                # Linear iff only numerator linear
                numer = operands.next()
                denom = operands.next()
                if self._check_expr(denom, state_var, bad_vars) != kind.None:
                    result = kind.Nonlinear
                else:
                    result = self._check_expr(numer, state_var, bad_vars)
            elif operator == 'times':
                # Linear iff only 1 linear operand
                op_kinds = map(lambda op: self._check_expr(op, state_var,
                                                           bad_vars),
                               operands)
                lin, nonlin = 0, 0
                for res in op_kinds:
                    if res == kind.Linear: lin += 1
                    elif res == kind.Nonlinear: nonlin += 1
                if nonlin > 0 or lin > 1:
                    result = kind.Nonlinear
                elif lin == 1:
                    result = kind.Linear
                else:
                    result = kind.None
            else:
                # Cannot be linear; may be no dependence at all
                result = max(map(lambda op: self._check_expr(op, state_var,
                                                             bad_vars),
                                 operands))
                if result == kind.Linear:
                    result = kind.Nonlinear
        else:
            # Either a straightforward container element
            try:
                child = child_i(expr, 1)
            except ValueError:
                # Assume it's just a constant
                result = kind.None
            else:
                result = self._check_expr(child, state_var, bad_vars)
        DEBUG('find-linear-deps', "Expression", expr, "gives result", result)
        return result

    def find_linear_odes(self, state_vars, V, free_var):
        """Identify linear ODEs.

        For each ODE (except that for V), determine whether it has a
        linear dependence on the dependent variable, and thus can be
        updated directly, without using Newton's method.

        We also require it to not depend on any other state variable,
        except for V.
        """
        kind = self.LINEAR_KINDS
        candidates = set(state_vars) - set([V])
        linear_vars = []
        for var in candidates:
            ode_expr = var.get_ode_dependency(free_var)
            if self._check_expr(self._get_rhs(ode_expr), var,
                                candidates - set([var])) == kind.Linear:
                linear_vars.append(var)
        return linear_vars

    def _clone(self, expr):
        """Properly clone a MathML sub-expression."""
        if isinstance(expr, mathml):
            clone = expr.clone_self(register=True)
        else:
            clone = mathml.clone(expr)
        return clone

    def _make_apply(self, operator, ghs, i, filter_none=True,
                    preserve=False):
        """Construct a new apply expression for g or h.

        ghs is an iterable of (g,h) pairs for operands.
        
        i indicates whether to construct g (0) or h (1).
        
        filter_none indicates the behaviour of 0 under this operator.
        If True, it's an additive zero, otherwise it's a
        multiplicative zero.
        """
        # Find g or h operands
        ghs_i = map(lambda gh: gh[i], ghs)
        if not filter_none and None in ghs_i:
            # Whole expr is None
            new_expr = None
        else:
            # Filter out None subexprs
            if operator == u'minus':
                # Do we need to retain a unary minus?
                if len(ghs_i) == 1 or ghs_i[0] is None:
                    # Original was -a or 0-a
                    retain_unary_minus = True
                else:
                    # Original was a-0 or a-b
                    retain_unary_minus = False
            else:
                # Only retain if we're told to preserve as-is
                retain_unary_minus = preserve
            ghs_i = filter(None, ghs_i)
            if ghs_i:
                if len(ghs_i) > 1 or retain_unary_minus:
                    new_expr = mathml_apply.create_new(
                        self.__expr, operator, ghs_i)
                else:
                    new_expr = self._clone(ghs_i[0]) # Clone may be unneeded
            else:
                new_expr = None
        return new_expr

    def _transfer_lut(self, expr, gh, var):
        """Transfer lookup table annotations from expr to gh.

        gh is a pair (g, h) s.t. expr = g + h*var.

        If expr can be represented by a lookup table, then the lookup
        variable cannot be var, since if it were, then expr would have
        a non-linear dependence on var.  Hence h must be 0, since
        otherwise expr would contain a (state) variable other than the
        lookup variable, and hence not be a candidate for a table.
        Thus expr=g, so we transfer the annotations to g.
        """
        if expr.getAttributeNS(NSS['lut'], u'possible', '') != u'yes':
            return
        # Paranoia check that our reasoning is correct
        g, h = gh
        assert h is None
        # Transfer the annotations into g
        LookupTableAnalyser.copy_lut_annotations(expr, g)
        # Make sure g has a reference to its component, for use by code generation.
        g._cml_component = expr.component
        return

    def _rearrange_expr(self, expr, var):
        """Rearrange an expression into the form g + h*var.

        Performs a post-order traversal of this expression's tree,
        and returns a pair (g, h)
        """
#        import inspect
#        depth = len(inspect.stack())
#        print ' '*depth, "_rearrange_expr", prid(expr, True), var.name, expr
        gh = None
        if isinstance(expr, mathml_ci):
            # Variable
            ci_var = expr.variable.get_source_variable(recurse=True)
            if var is ci_var:
                gh = (None, mathml_cn.create_new(expr,
                                                 u'1', u'dimensionless'))
            else:
                if ci_var._cml_linear_kind == self.LINEAR_KINDS.None:
                    # Just put the <ci> in g, but with full name
                    gh = (mathml_ci.create_new(expr, ci_var.fullname()), None)
                    gh[0]._set_variable_obj(ci_var)
                else:
                    # ci_var is a linear function of var, so rearrange
                    # its definition
                    if not hasattr(ci_var, '_cml_linear_split'):
                        ci_defn = ci_var.get_dependencies()[0]
                        ci_var._cml_linear_split = self._rearrange_expr(
                            self._get_rhs(ci_defn), var)
                    gh = ci_var._cml_linear_split
        elif isinstance(expr, mathml_piecewise):
            # The tests have to move into both components of gh:
            # "if C1 then (a1,b1) elif C2 then (a2,b2) else (a0,b0)"
            # maps to "(if C1 then a1 elif C2 then a2 else a0,
            #           if C1 then b1 elif C2 then b2 else b0)"
            # Note that no test is a function of var.
            # First rearrange child expressions
            pieces = getattr(expr, u'piece', [])
            cases = map(lambda p: child_i(p, 1), pieces)
            cases_ghs = map(lambda c: self._rearrange_expr(c, var), cases)
            if hasattr(expr, u'otherwise'):
                ow_gh = self._rearrange_expr(child_i(expr.otherwise, 1), var)
            else:
                ow_gh = (None, None)
            # Now construct the new expression
            conds = map(lambda p: self._clone(child_i(p, 2)), pieces)
            def piecewise_branch(i):
                pieces_i = zip(map(lambda gh: gh[i], cases_ghs),
                               conds)
                pieces_i = filter(lambda p: p[0] is not None,
                                  pieces_i) # Remove cases that are None
                ow = ow_gh[i]
                if pieces_i:
                    new_expr = mathml_piecewise.create_new(
                        expr, pieces_i, ow)
                elif ow:
                    new_expr = ow
                else:
                    new_expr = None
                return new_expr
            gh = (piecewise_branch(0), piecewise_branch(1))
            self._transfer_lut(expr, gh, var)
        elif isinstance(expr, mathml_apply):
            # Behaviour depends on the operator
            operator = expr.operator().localName
            operands = expr.operands()
            self.__expr = expr # For self._make_apply
            if operator in ['plus', 'minus']:
                # Just split the operation into each component
                operand_ghs = map(lambda op: self._rearrange_expr(op, var),
                                  operands)
                g = self._make_apply(operator, operand_ghs, 0)
                h = self._make_apply(operator, operand_ghs, 1)
                gh = (g, h)
            elif operator == 'divide':
                # (a, b) / (c, 0) = (a/c, b/c)
                numer = self._rearrange_expr(operands.next(), var)
                denom = self._rearrange_expr(operands.next(), var)
                assert denom[1] is None
                denom_g = denom[0]
                g = h = None
                if numer[0]:
                    g = mathml_apply.create_new(expr, operator,
                                                [numer[0], denom_g])
                if numer[1]:
                    if g:
                        denom_g = self._clone(denom_g)
                    h = mathml_apply.create_new(expr, operator,
                                                [numer[1], denom_g])
                gh = (g, h)
            elif operator == 'times':
                # (a1,b1)*(a2,b2) = (a1*a2, b1*a2 or a1*b2 or None)
                # Similarly for the nary case - at most one b_i is not None
                operand_ghs = map(lambda op: self._rearrange_expr(op, var),
                                  operands)
                g = self._make_apply(operator, operand_ghs, 0,
                                     filter_none=False)
                # Find non-None b_i, if any
                for i, ab in enumerate(operand_ghs):
                    if ab[1] is not None:
                        operand_ghs[i] = (ab[1], None)
                        # Clone the a_i to avoid objects having 2 parents
                        for j, ab in enumerate(operand_ghs):
                            if j != i:
                                operand_ghs[j] = (self._clone(operand_ghs[j][0]), None)
                        h = self._make_apply(operator, operand_ghs, 0, filter_none=False)
                        break
                else:
                    h = None
                gh = (g, h)
            else:
                # (a, None) op (b, None) = (a op b, None)
                operand_ghs = map(lambda op: self._rearrange_expr(op, var),
                                  operands)
                g = self._make_apply(operator, operand_ghs, 0, preserve=True)
                gh = (g, None)
            self._transfer_lut(expr, gh, var)
        else:
            # Since this expression is linear, there can't be any
            # occurrence of var in it; all possible such cases are covered
            # above.  So just clone it into g.
            gh = (self._clone(expr), None)
#        print ' '*depth, "Re-arranged", prid(expr, True), "to", prid(gh[0], True), ",", prid(gh[1], True)
        return gh
    
    def rearrange_linear_odes(self, doc):
        """Rearrange the linear ODEs so they can be updated directly
        on solving.

        Each ODE du/dt = f(u, t) can be written in the form
            du/dt = g(t) + h(t)u.
        A backward Euler update step is then as simple as
            u_n = (u_{n-1} + g(t)dt) / (1 - h(t)dt)
        (assuming that the transmembrane potential has already been
        solved for at t_n.

        Stores the results in doc.model._cml_linear_update_exprs, a
        mapping from variable object u to pair (g, h).
        """

        odes = map(lambda v: v.get_ode_dependency(doc.model._cml_free_var),
                   doc.model._cml_linear_vars)
        result = {}
        for var, ode in itertools.izip(doc.model._cml_linear_vars, odes):
            # Do the re-arrangement for this variable.
            # We do this by a post-order traversal of its defining ODE,
            # constructing a pair (g, h) for each subexpression recursively.
            # First, get the RHS of the ODE
            rhs = self._get_rhs(ode)
            # And traverse
            result[var] = self._rearrange_expr(rhs, var)
        # Store result in model
        doc.model._cml_linear_update_exprs = result
        return result

    def show(self, d):
        """Print out a more readable report on a rearrangement,
        as given by self.rearrange_linear_odes."""
        for var, expr in d.iteritems():
            print var.fullname()
            print "G:", expr[0].xml()
            print "H:", expr[1].xml()
            print "ODE:", var.get_ode_dependency(
                var.model._cml_free_var).xml()
            print


######################################################################
#                        Rush-Larsen analysis                        #
######################################################################

class ExpressionMatcher(object):
    """Test whether a MathML expression matches a given tree pattern.
    
    Patterns are instances of the nested Pattern class, or more specifically
    one of its subclasses.  The static method match on this class checks an
    expression against a pattern, returning True iff there is a match.
    """
    
    class Pattern(object):
        """Abstract base class for tree patterns."""
        def match(self, expr):
            """
            Method implemented by concrete subclasses to test a given expression.
            Returns True iff there is a match.
            """
            raise NotImplementedError
    
    class A(Pattern):
        """An apply expression."""
        def __init__(self, operator, operands):
            self.operator = operator
            self.operands = operands
        
        def match(self, expr):
            matched = False
            if isinstance(expr, mathml_apply):
                if expr.operator().localName == self.operator:
                    expr_operands = list(expr.operands())
                    if len(expr_operands) == len(self.operands):
                        matched = reduce(operator.and_,
                                         map(lambda (pat, op): pat.match(op),
                                             zip(self.operands, expr_operands)))
            return matched
    
    class V(Pattern):
        """A variable reference."""
        def __init__(self, var=None):
            self.set_variable(var)
        
        def set_variable(self, var):
            if var:
                self.var = var.get_source_variable(recurse=True)
            else:
                self.var = var
        
        def match(self, expr):
            matched = False
            if isinstance(expr, mathml_ci) and self.var is expr.variable.get_source_variable(recurse=True):
                matched = True
            return matched
    
    class N(Pattern):
        """A constant number, optionally with the value specified."""
        def __init__(self, value=None):
            self.value = value
        
        def match(self, expr):
            matched = False
            if isinstance(expr, mathml_cn):
                value = expr.evaluate()
                if self.value is None:
                    self.value = value
                if self.value == value:
                    matched = True
            return matched
    
    class X(Pattern):
        """A placeholder, matching anything (and noting what was matched)."""
        def __init__(self):
            self.matched = None
        
        def match(self, expr):
            self.matched = expr
            return True
    
    class R(Pattern):
        """A container that matches any number of levels of indirection/recursion.
        
        This can be used to wrap a pattern where we wish to allow for variable mappings
        or equations such as "var1 = var2" before we reach the 'interesting' equation.
        If the expression we're matching is a ci element we recursively find the
        ultimate non-ci defining expression and match our sub-pattern against that.  If
        the expression isn't a ci, or the ultimate definition isn't an expression, we
        match our sub-pattern against it directly.
        """
        def __init__(self, pat):
            self.sub_pattern = pat
        
        def match(self, expr):
            while isinstance(expr, mathml_ci):
                # Find this variable's defining expression, if it is an equation
                var = expr.variable.get_source_variable(recurse=True)
                defn = var.get_dependencies()
                if defn and isinstance(defn[0], mathml_apply):
                    expr = defn[0].eq.rhs
            return self.sub_pattern.match(expr)
    
    @staticmethod
    def match(pattern, expression):
        """Test for a match."""
        return pattern.match(expression)

class RushLarsenAnalyser(object):
    """Analyse a model to identify Hodgkin-Huxley style gating variable equations.
    
    We look for ODEs whose definition matches "alpha*(1-x) - beta*x" (where x is
    the state variable, and alpha & beta are any expression).  Alternatively for
    models which already have tau & inf variables, we match against "(inf-x)/tau".
    
    To allow for when units conversions have been performed, we chase 'simple'
    assignments and (the semantically equivalent) variable mappings until reaching
    an 'interesting' defining equation.  We also need to allow the whole RHS to be
    multiplied by a constant.  If this occurs, the constant conversion factor is
    also stored; otherwise we store None.
    
    Stores a dictionary on the document root mapping cellml_variable instances to
    4-tuples.  The tuple is either ('ab', alpha, beta, conv) or ('ti', tau, inf, conv)
    depending on which formulation has been used.  Note that the expressions in these
    are not cloned copies - they are the original objects still embedded within the
    relevant ODE.  The units conversion factor 'conv' is stored as a Python double.
    """
    def __init__(self):
        """Create the patterns to match against."""
        em = ExpressionMatcher
        self._var = em.V()
        self._conv = em.N()
        # Alpha/beta form
        self._alpha = em.X()
        self._beta = em.X()
        self._ab_pattern = em.R(em.A('minus', [em.A('times', [self._alpha,
                                                              em.A('minus', [em.N(1), self._var])]),
                                               em.A('times', [self._beta, self._var])]))
        self._alt_ab_pattern = em.R(em.A('times', [self._conv, self._ab_pattern]))
        # Tau/inf form
        self._tau = em.X()
        self._inf = em.X()
        self._ti_pattern = em.R(em.A('divide', [em.A('minus', [self._inf, self._var]),
                                                self._tau]))
        self._alt_ti_pattern = em.R(em.A('times', [self._conv, self._ti_pattern]))
    
    def analyse_model(self, doc):
        # First, find linear ODEs that have the potential to be gating variables
        la = LinearityAnalyser()
        V = doc._cml_config.V_variable
        state_vars = doc.model.find_state_vars()
        free_var = doc.model.find_free_vars()[0]
        linear_vars = la.find_linear_odes(state_vars, V, free_var)
        # Next, check they match dn/dt = a (1-n) - b n
        doc._cml_rush_larsen = {}
        for var in linear_vars:
            ode_expr = var.get_ode_dependency(free_var)
            self._check_var(var, ode_expr, doc._cml_rush_larsen)
    
    def _check_var(self, var, ode_expr, mapping):
        rhs = ode_expr.eq.rhs
        self._var.set_variable(ode_expr.eq.lhs.diff.dependent_variable)
        if self._ab_pattern.match(rhs):
            mapping[var] = ('ab', self._alpha.matched, self._beta.matched, None)
        elif self._alt_ab_pattern.match(rhs):
            mapping[var] = ('ab', self._alpha.matched, self._beta.matched, self._conv.value)
        elif self._ti_pattern.match(rhs):
            mapping[var] = ('ti', self._tau.matched, self._inf.matched, None)
        elif self._alt_ti_pattern.match(rhs):
            mapping[var] = ('ti', self._tau.matched, self._inf.matched, self._conv.value)
