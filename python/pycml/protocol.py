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
Defines the Protocol class, which encapsulates the interface between a model
and a Web Lab protocol.
"""

import os
import sys

import pycml
from pycml import *
import processors

class ProtocolError(ValueError):
    """Error thrown if a Protocol instance is invalid."""
    pass

class Protocol(processors.ModelModifier):
    """A class representing part of a simulation protocol for the functional curation system.

    PyCml is responsible for implementing the 'model interface' section of protocols.
    See https://chaste.cs.ox.ac.uk/trac/wiki/SimulationProtocolNotes#Modelinterface
    for further details.  The main user-facing methods here implement the XML elements defined
    there:
     - specify_input_variable
     - specify_output_variable
     - set_independent_variable_units
     - declare_new_variable
     - add_or_replace_equation
     - define_units_conversion_rule
    """
    def __init__(self, model, multi_stage=True, namespaces={}):
        """Create a new protocol.
        
        The public methods listed above should be called to set up the internal data structures,
        and then self.modify_model called to apply the protocol to the model.
        """
        super(Protocol, self).__init__(model)
        self._protocol_component = None
        self._free_var_has_changed = None
        self._protocol_namespaces = model._cml_protocol_namespaces = {}
        self.add_protocol_namespaces(namespaces)
        self.inputs = set()
        self._input_specifications = []
        self.outputs = set()
        self._output_specifications = []
        self._vector_outputs = set()
        self._vector_outputs_detail = []
        self._optional_vars = set()
        warn_only = not model.get_option('fully_automatic') and model.get_option('warn_on_units_errors')
        self.set_units_converter(processors.UnitsConverter(self.model, warn_only))
        self.assignments_to_convert = set()
    
    @staticmethod
    def apply_protocol_file(doc, proto_file_path):
        """Parse a protocol XML file and apply it to the given model document."""
        proto = Protocol(doc.model)
        proto.units = doc.model.get_standard_units().copy()
        proto.parse_protocol(proto_file_path, proto.units)
        proto.modify_model()
    
    @staticmethod
    def find_required_annotations(proto_file_path):
        """Parse a protocol XML file and determine what model annotations will be required.
        
        Note: this method is obsolete, and doesn't pick up everything!
        """
        proto_docs = []
        model_nss = {}
        # First load all the XMLs
        def load_proto(path):
            proto_xml = amara_parse_cellml(path)
            assert hasattr(proto_xml, u'protocol')
            model_nss.update(proto_xml.xmlns_prefixes)
            proto_docs.append(proto_xml)
            for proto_import in getattr(proto_xml.protocol, u'import_', []):
                # Relative URIs must be resolved relative to this protocol file
                source = proto_import.source
                if not os.path.isabs(source):
                    source = os.path.join(os.path.dirname(path), source)
                load_proto(source)
        load_proto(proto_file_path)
        # Now analyse them
        def get_var_names(proto):
            var_names = set([unicode(e).strip() for e in
                             proto.xml_xpath(u'//m:ci') + proto.xml_xpath(u'//proto:name')
                             + proto.xml_xpath(u'/*/*/proto:specifyOutputVariable/@name')])
            optional_names = set()
            for input in proto.xml_xpath(u'/*/*/proto:specifyInputVariable'):
                if hasattr(input, u'initial_value'):
                    optional_names.add(unicode(input.name).strip())
                else:
                    var_names.add(unicode(input.name).strip())
            var_names -= optional_names
            return var_names, optional_names
        def check_names(names, term_set):
            for varname in names:
                if u':' in varname:
                    parts = varname.split(':')
                    if parts[-2] in model_nss:
                        term_set.add((parts[-2], parts[-1]))
        terms, opt_terms = set(), set()
        for proto in proto_docs:
            var_names, optional_names = get_var_names(proto)
            check_names(var_names, terms)
            check_names(optional_names, opt_terms)
        # Display results
        print "Terms required by protocol %s:" % proto_file_path
        for term in sorted(terms):
            print "   ", term[0] + ':' + term[1]
        if opt_terms:
            print "Optional terms:"
            for term in sorted(opt_terms):
                print "   ", term[0] + ':' + term[1]

    def parse_protocol(self, proto_file_path, proto_units, prefix='', units_only=False):
        """Parse a protocol XML file and set up our data structures accordingly."""
        proto_xml = amara_parse_cellml(proto_file_path)
        assert hasattr(proto_xml, u'protocol')
        self.add_protocol_namespaces(proto_xml.xmlns_prefixes)
        # Relative URIs must be resolved relative to this protocol file, or its xml:base if present
        base = os.path.dirname(getattr(proto_xml.protocol, 'base', proto_file_path))
        self.base = getattr(self, 'base', base) # Only set this if we're the main file
        # Any imports?
        for proto_import in getattr(proto_xml.protocol, u'import_', []):
            source = proto_import.source
            if not os.path.isabs(source):
                source = os.path.join(base, source)
            if not os.path.exists(source):
                # Try resolving relative to the library folder
                library = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir,
                                       'projects', 'FunctionalCuration', 'src', 'proto', 'library')
                source = os.path.join(library, proto_import.source)
            if getattr(proto_import, u'mergeDefinitions', u'0') in [u'true', u'1']:
                # Process this import immediately
                self.parse_protocol(source, proto_units)
            else:
                # Only apply model modifications from the import if requested, but
                # make all units definitions available with the prefix.
                import_units_only = True
                if hasattr(proto_xml.protocol, u'modelInterface'):
                    for use in getattr(proto_xml.protocol.modelInterface, u'useImports', []):
                        if use.prefix_ == proto_import.prefix_:
                            import_units_only = False
                            break
                self.parse_protocol(source, proto_units, prefix=proto_import.prefix_, units_only=import_units_only)
        # For now, we also need to apply model modifications from nested protocols
        for nested_proto in proto_xml.protocol.xml_xpath(u'.//proto:nestedProtocol'):
            source = nested_proto.source
            if not os.path.isabs(source):
                source = os.path.join(base, source)
            self.parse_protocol(source, self.model.get_standard_units().copy())
        if hasattr(proto_xml.protocol, u'units'):
            # Parse units definitions
            for defn in getattr(proto_xml.protocol.units, u'units', []):
                uname = defn.name
                if prefix:
                    uname = prefix + ':' + uname
                if uname in proto_units:
                    raise ProtocolError("Duplicate definition of units named '%s'" % uname)
                proto_units[uname] = defn
                defn.xml_parent = self.model
                self.add_units(defn)
        def get_units(elt, attr='units'):
            if hasattr(elt, attr):
                uname = getattr(elt, attr)
                try:
                    if not ':' in uname and prefix:
                        uname = prefix + ':' + uname
                    return proto_units[uname]
                except KeyError:
                    raise ProtocolError("Units '%s' have not been defined in the protocol" % uname)
            else:
                return None
        if not units_only and hasattr(proto_xml.protocol, u'modelInterface'):
            for optional in getattr(proto_xml.protocol.modelInterface, u'specifyOptionalVariable', []):
                self.specify_optional_variable(optional.name, optional.xml_children)
            for vardecl in getattr(proto_xml.protocol.modelInterface, u'declareNewVariable', []):
                self.declare_new_variable(vardecl.name, get_units(vardecl), getattr(vardecl, u'initial_value', None))
            for rule in getattr(proto_xml.protocol.modelInterface, u'unitsConversionRule', []):
                self.add_units_conversion_rule(get_units(rule, 'actualDimensions'),
                                               get_units(rule, 'desiredDimensions'),
                                               rule.xml_element_children().next())
            if hasattr(proto_xml.protocol.modelInterface, u'setIndependentVariableUnits'):
                self.set_independent_variable_units(get_units(proto_xml.protocol.modelInterface.setIndependentVariableUnits))
            for input in getattr(proto_xml.protocol.modelInterface, u'specifyInputVariable', []):
                self.specify_input_variable(input.name, get_units(input), getattr(input, u'initial_value', None))
            for output in getattr(proto_xml.protocol.modelInterface, u'specifyOutputVariable', []):
                self.specify_output_variable(output.name, get_units(output))
            for expr in getattr(proto_xml.protocol.modelInterface, u'addOrReplaceEquation', []):
                self.add_or_replace_equation(expr.xml_element_children().next())
    
    def specify_optional_variable(self, prefixed_name, children):
        """Specify the given variable as being optional, so that the interface copes if it isn't present.
        
        If a default expression is given, we also need to add this equation to the inputs set iff the LHS
        variable doesn't exist.  The variable itself will be added in a later processing phase, once we
        can determine whether its definition can actually be evaluated.  It will be given units even later,
        once variable connections have been fixed and hence the units of the RHS can be calculated.

        As a special case, however, if the default is just a number with units, we treat this as declaring
        a variable using the initial_value attribute, and create the whole thing here, to avoid confusing
        edge cases where the optional variable is a state variable with ODE specified using `define diff()`.
        """
        self._optional_vars.add(prefixed_name)
        assert len(children) <= 1, 'Malformed specifyOptionalVariable element with %d children' % len(children)
        if len(children) == 1:
            # Default given; check if we need it
            if self._lookup_ontology_term(prefixed_name, check_optional=True) is None:
                rhs = children[0]
                if self.expr_is_simple_constant(rhs) and self._get_cn_units(rhs):
                    # Special simple case
                    units = self.model.get_units_by_name(self._get_cn_units(rhs))
                    var = self._create_annotated_variable(prefixed_name, units)
                    var.initial_value = unicode(rhs.evaluate())
                else:
                    # General equation definition
                    rhs.xml_parent.safe_remove_child(rhs)
                    lhs = pycml.mathml_ci.create_new(rhs, prefixed_name)
                    defn = pycml.mathml_apply.create_new(rhs, u'eq', [lhs, rhs])
                    self.add_or_replace_equation(defn)
    
    def specify_output_variable(self, prefixed_name, units=None):
        """Specify the given variable as a protocol output, optionally in the given units.
        
        This method just notes the details given, ready for final processing after the model equations have been
        modified (if needed) and in particular after the input variable declarations have been processed.  This
        ordering is necessary to allow a variable to be both an input and an output.
        """
        self._output_specifications.append({'prefixed_name': prefixed_name, 'units': units, 'optional': prefixed_name in self._optional_vars})

    def process_output_declarations(self):
        """Finish processing output variable declarations after modifying the model equations and processing inputs.
        
        For each output, the units it was requested in (if any were specified) must be added to the model if they
        don't exist.  If they differ from the variable's original units, a conversion will be needed, and hence a
        new version of the variable will be added to the protocol component, with a suitable connection.  Otherwise,
        we can just record the existing variable as an output.
        
        When an ontology term matches multiple variables, then these should be treated as a single vector output.
        These structures are set up after considering single outputs, in case an output occurs both individually
        and within a vector.
        """
        # Annotate state variables with oxmeta:state_variable
        prop, targ = ('bqbiol:isVersionOf', NSS['bqbiol']), ('oxmeta:state_variable', NSS['oxmeta'])
        cellml_metadata.remove_statements(self.model, None, prop, targ)
        for var in self._find_state_variables():
            var.add_rdf_annotation(prop, targ, allow_dup=True)
        for output_spec in self._output_specifications:
            prefixed_name = output_spec['prefixed_name']
            units = output_spec['units']
            try:
                vars = self._lookup_ontology_term(prefixed_name, enforce_uniqueness=False, check_optional=True, transitive=True)
            except ValueError, e:
                raise ProtocolError(str(e))
            if vars is None:
                print >>sys.stderr, 'Ignoring missing optional output', prefixed_name
            elif len(vars) > 1:
                self._vector_outputs_detail.append((prefixed_name, units))
            else:
                var = vars[0]
                if units is None:
                    units = var.get_units()
                new_var = self.specify_as_output(var, units)
                if var.get_type() is VarTypes.Free:
                    self.set_independent_variable_units(units, new_var)
        self.process_output_variable_vectors()

    def process_output_variable_vectors(self):
        """Finish adding outputs that are vectors of variables.
        
        When an ontology term given to specify_output_variable matches multiple variables,
        then these should be treated as a single vector output.  This method sets up these
        structures.  It's called after all output variables have been specified, in case an
        output occurs both individually and within a vector.
        """
        # Re-lookup all the ontology terms that matched multiple variables
        for prefixed_name, units in self._vector_outputs_detail:
            vars = self._lookup_ontology_term(prefixed_name, enforce_uniqueness=False, transitive=True)
            vector_name = prefixed_name.split(':')[1]
            for var in vars:
                # Units convert if needed
                desired_units = self._get_units_object(units or var.get_units())
                if not desired_units.equals(var.get_units()):
                    if var.component is self._get_protocol_component():
                        raise ProtocolError("You can't ask for an output (%s) in two different units!" % var.fullname())
                    new_var = self._replace_variable(var, desired_units)
                    self.connect_variables(var, new_var)
                    var = new_var
                # Ensure it gets computed and annotated
                self._vector_outputs.add(var)
                var.add_rdf_annotation(('pycml:output-vector', NSS['pycml']), vector_name)
    
    def _create_annotated_variable(self, prefixed_name, units):
        """Create a new variable in the model, annotated with the given term, and in the given units."""
        #1903 TODO: Be more careful to create unique local names and ids
        prefix, local_name = prefixed_name.split(':')
        var = self.add_variable(self._get_protocol_component(), local_name, units, id=prefix + '_' + local_name)
        var.add_rdf_annotation(('bqbiol:is', NSS['bqbiol']), (prefixed_name, self._protocol_namespaces[prefix]))
        return var
    
    def _find_state_variables(self):
        """Find all (likely) state variables by examining the model mathematics for ODEs."""
        state_vars = []
        for expr in self.model.search_for_assignments():
            ode_vars = self._is_ode(expr, local_vars=True)
            if ode_vars:
                state_vars.append(ode_vars[0].get_source_variable(recurse=True))
        return state_vars
    
    def _is_ode(self, expr, local_vars=False):
        """Determine whether the given assignment expression is an ODE.
        
        If it is, return the tuple (dependent_var, independent_var).  Otherwise, return None.
        """
        assert isinstance(expr, mathml_apply)
        assert expr.operator().localName == u'eq', 'Expression is not an assignment'
        # Figure out what's on the LHS of the assignment
        lhs, rhs = list(expr.operands())
        if lhs.localName == u'ci':
            return None
        assert lhs.localName == u'apply', 'Expression is neither a straight assignment nor an ODE'
        assert lhs.operator().localName == u'diff', 'Expression is neither a straight assignment nor an ODE'
        dep_var = lhs.operands().next()
        assert dep_var.localName == u'ci', 'ODE is malformed'
        def get_var(comp, name):
            if u',' in name or not local_vars:
                var = self.model.get_variable_by_name(*self._split_name(name))
            else:
                var = comp.get_variable_by_name(name)
            return var
        comp = expr.component
        dep_var = get_var(comp, unicode(dep_var))
        indep_var = get_var(comp, unicode(lhs.bvar.ci))
        return (dep_var, indep_var)
    
    def expr_is_simple_constant(self, expr):
        """Check whether the given expression is a straight cn or -cn."""
        return (isinstance(expr, mathml_cn) or
                (isinstance(expr, mathml_apply) and hasattr(expr, u'minus') and len(expr.xml_children) == 2
                 and isinstance(expr.xml_children[1], mathml_cn)))
    
    def process_input_declarations(self):
        """Finish processing input declarations after modifying the model equations.
        
        For each variable (identified by ontology term) specified as a model input, we ensure it
        has the correct units and initial value.
        
        The variable must be a state variable, constant, or free variable, not computed by an
        equation.  As a special case, a variable with type Computed that is defined by assigning
        a single number is converted to be a Constant.
        
        If the units specified differ from the original units, they must be added to the model if
        not present, and a new version of the variable added to the protocol component, with a
        suitable assignment.  Connections and ODEs will be updated if needed.
        
        If an initial_value is given then this will be assigned as the attribute value.
        """
        for input in self._input_specifications:
            try:
                var = self._lookup_ontology_term(input['prefixed_name'])
            except ValueError:
                if input['prefixed_name'] not in self._optional_vars:
                    raise ProtocolError("There is no model variable annotated with the term " + prefixed_name)
                else:
                    print >>sys.stderr, 'Ignoring missing optional input', input['prefixed_name']
                    continue
            # Check it has the correct type
            deps = var.get_all_expr_dependencies()
            if len(deps) > 1:
                raise ProtocolError("Variable " + str(var) + " has been over-specified.")
            if len(deps) == 1:
                if isinstance(deps[0], cellml_variable):
                    raise ProtocolError("The mapped variable " + str(var) + " may not be specified as an input - annotate its source instead.")
                assert isinstance(deps[0], mathml_apply)
                lhs, rhs = list(deps[0].operands())
                if lhs.localName == u'ci':
                    if not self.expr_is_simple_constant(rhs):
                        raise ProtocolError("The computed variable " + str(var) + " may not be specified as an input.")
                    # It's the special Computed case - convert to a constant
                    initial_value = rhs.evaluate()
                    if not var.get_units().equals(rhs.get_units()):
                        converter = self.get_units_converter()
                        initial_value = converter.convert_constant(initial_value, rhs.get_units().extract(), var.get_units(), var.component)
                    var.initial_value = unicode(initial_value)
                    self.remove_definition(var, keep_initial_value=True)
                    if deps[0] in self.inputs:
                        self.inputs.remove(deps[0])
                else:
                    assert lhs.operator().localName == u'diff' # It's an ODE
            if not input['initial_value'] and not hasattr(var, u'initial_value'):
                raise ProtocolError("No value specified for protocol input " + input['prefixed_name'])
            # Convert units if needed
            units = self._get_units_object(input['units'] or var.get_units())
            if units.equals(var.get_units()):
                input_var = var
            else:
                input_var = self._replace_variable(var, units)
                if not input['initial_value']:
                    input_var.initial_value = unicode(self._convert_initial_value(var, units))
                self.del_attr(var, u'initial_value', None)
                # Set all variables connected to the original variable (including itself) to be mapped to the new one
                # TODO: Would connect_variables be sufficient here?
                self._update_connections(var, input_var)
                # TODO: Do we need to do anything special for state variables?
            # Update initial value if specified
            if input['initial_value']:
                input_var.initial_value = unicode(input['initial_value'])
            # Ensure the name annotation is correct
            prefix = input['prefixed_name'].split(':')[0]
            input_var.add_rdf_annotation(('bqbiol:is', NSS['bqbiol']), (input['prefixed_name'], self._protocol_namespaces[prefix]))
            # Add to the old self.inputs collection for statistics calculation
            self.inputs.add(input_var)
            if var is not input_var:
#                 print 'Removing old input', var, 'new is', input_var
                self.inputs.remove(var)

    def specify_input_variable(self, prefixed_name, units=None, initial_value=None):
        """Set the given variable as a protocol input, optionally in the given units.
        
        At this point we just note that it's going to become an input, and create the variable if
        it doesn't exist (so that a later add_or_replace_equation call can define it if appropriate),
        unless it is also specified as being optional.
        
        The model definition will only be replaced if specified using add_or_replace_equation.
        After those parts have been processed, we will iterate through all specified inputs and
        process the units and initial_value specifications, using the process_input_declarations
        method.
        
        If units are given and differ from its original units, they must be added to the model if they don't
        exist, and a new version of the variable added to the protocol component, with a suitable assignment.
        
        If the initial_value is given then this will overwrite the original setting if present.
        
        If the variable does not exist, this is only an error if an initial_value or units are not given,
        and the variable is not declared as optional.  Otherwise we just create the variable.
        """
        try:
            var = self._lookup_ontology_term(prefixed_name)
        except ValueError:
            if prefixed_name in self._optional_vars:
                print >>sys.stderr, 'Ignoring missing optional input', prefixed_name
                return
            else:
                # Create the variable
                if units is None:
                    raise ProtocolError("Units must be specified for input variables not appearing in the model; none are given for " + prefixed_name)
                var = self._create_annotated_variable(prefixed_name, units)
        # Flag for later processing
        self._input_specifications.append({'prefixed_name': prefixed_name, 'units': units, 'initial_value': initial_value})
        # Record for statistics output
        var._cml_ok_as_input = True
        self.inputs.add(var)
    
    def set_independent_variable_units(self, units, newVar=None):
        """Set the independent variable to occur in the given units.
        
        If newVar is given, then this is being called by self.specify_output_variable, since the
        independent variable is also a protocol output, and the new version has already been
        created.
        
        Since this may mean we're called twice, we ensure the second call is a no-op.
        """
        if self._free_var_has_changed:
            assert self._free_var_has_changed.get_units().equals(units)
            assert newVar is None or newVar is self._free_var_has_changed
        else:
            t = self.model.find_free_vars()[0]
            units = self._get_units_object(units)
            if not units.equals(t.get_units()):
                # We'll need a conversion, and to convert all ODEs too
                self._free_var_has_changed = newVar or self._replace_variable(t, units)
    
    def declare_new_variable(self, name, units, initial_value=None):
        """Declare a new variable for use in the model interface.
        
        The variable will be added to the protocol component, and the units added to the model
        if they're not already present.  An assertion will be tripped if a variable with the
        given name already exists.
        """
        var = self.add_variable(self._get_protocol_component(), name, units)
        if initial_value:
            var.initial_value = unicode(initial_value)
        return var
    
    def add_or_replace_equation(self, assignment):
        """Add the given assignment equation to the model.
        
        It will replace any existing definition of the same variable.
        """
        assert isinstance(assignment, mathml_apply)
        assert assignment.operator().localName == u'eq'
        self.inputs.add(assignment)
    
    def add_units_conversion_rule(self, from_units, to_units, conv_expr):
        """Add a new biology-aware units conversion rule.
        
        The third argument must be an instance of mathml_lambda taking a single parameter.
        It will effectively be passed the RHS of an assignment expression, which has units
        dimensionally equivalent to from_units, and must return an expression with units
        dimensionally equivalent to to_units.
        
        Note that the UnitsConverter class needs a function object that will modify the
        assignment equation in-place, so that's what we create and store.
        """
        converter = self.get_units_converter()
        children = list(conv_expr.xml_element_children())
        assert len(children) == 2, "A units conversion rule must have a single bound variable: " + conv_expr.xml()
        assert children[0].localName == u'bvar', "A units conversion rule must have a single bound variable: " + conv_expr.xml()
        bvar_name = unicode(children[0].ci).strip()
        body_expr = children[1]
        # Modify variable references within the body_expr so they use fully qualified names (compname, varname),
        # except for uses of the bound variable
        try:
            self._identify_referenced_variables(body_expr, bvar_name)
        except ValueError:
            # We can't apply this rule as some required variables are missing
            print >>sys.stderr, "Warning: unable to utilise units conversion rule below as required variables missing;",
            print >>sys.stderr, "this may lead to later units conversion errors:"
            if hasattr(conv_expr, 'loc'):
                print >>sys.stderr, "  ", conv_expr.loc
            else:
                print >>sys.stderr, "  From", from_units.description(), "to", to_units.description(), "via", conv_expr.xml()
            return
        func = lambda expr: self._apply_conversion_rule(expr, body_expr, bvar_name)
        converter.add_special_conversion(from_units, to_units, func)
    
    def add_protocol_namespaces(self, mapping):
        """Add to the prefix->URI mapping used by the protocol file."""
        self._protocol_namespaces.update(mapping)

    def modify_model(self):
        """Actually apply protocol modifications to the model.
        
        Prior to this being called, all variable references within self.inputs must
        use full names, i.e. 'component_name,variable_name'.  Variables without a
        component part will be placed into a new 'protocol' component.
        
        This method will add the items from self.inputs into the model, replacing
        variables with the same name, and equations that assign to the same variable.
        This may involve changing a variable's type from State to Computed, or vice
        versa.
        
        After the call, all names and name references will be 'local'.  Connections
        will be created between components as needed, and units definitions added, to
        ensure a valid model.  In order for this to work, if a variable has units that
        do not already exist in the model, the object *must* have an attribute
        _cml_units referring to a suitable cellml_units instance.
        
        Finally, the protocol outputs will be used to prune the model's assignments
        list so only assignments of interest are used to generate code.
        """
        # Add units before variables before maths so the order of inputs etc. doesn't matter so much.
        for input in filter(lambda i: isinstance(i, cellml_units), self.inputs):
            #self._check_input(input)
            self.add_units(input)
        for input in filter(lambda i: isinstance(i, cellml_variable), self.inputs):
            self._check_input(input)
            self._add_variable_to_model(input)
        for input in filter(lambda i: isinstance(i, mathml_apply), self.inputs):
            self._check_input(input)
            self._check_equation_lhs(input) # Create LHS var if needed
        for input in list(filter(lambda i: isinstance(i, mathml_apply), self.inputs)):
            # Note that we convert to list this time, because _add_maths_to_model might modify self.inputs
            self._check_input(input)
            self._add_maths_to_model(input)
        self.process_input_declarations()
        self.process_output_declarations()
        if self._free_var_has_changed:
            self._split_all_odes()
        self._fix_model_connections()
        self.finalize(self._error_handler, self._add_units_conversions)
        self._filter_assignments()
        self.report_stats()
    
    @property
    def magic_units(self):
        """Return a fake units object for use when creating variables whose units are not yet known.
        
        Also create a member list for recording variables with these units.
        """
        try:
            return self._magic_units
        except AttributeError:
            u = self._magic_units = pycml.cellml_units.create_new(self.model, u'**magic**', [])
            return u
    
    def _fix_magic_units(self):
        """Give real units (calculated from RHS of defining equation) to any variable with 'magic' units.
        
        We need to process variables in order sorted by the dependency graph, since some variables with magic units
        might be defined in terms of other variables with magic units.
        
        We also need to update the units of variables mapped to those we fix.
        """
        magic_vars = set()
        for expr in self.model.get_assignments():
            if isinstance(expr, mathml_apply):
                var = expr.assigned_variable()
                if isinstance(var, cellml_variable) and var.get_units() is self.magic_units:
                    defn = expr.eq.rhs
#                     print 'Magic units for', var, 'defined by', defn, defn.get_units()
                    units = self.add_units(defn.get_units().extract())
                    var.units = units.name
                    magic_vars.add(var)
                elif isinstance(var, tuple):
                    # It's an ODE; check the dependent var only
                    (dep_var, indep_var) = var
                    if dep_var.get_units() is self.magic_units:
                        defn = expr.eq.rhs
#                         print 'Magic ODE units for', dep_var, 'defined by', defn, defn.get_units()
                        rhs_units = defn.get_units().extract()
                        units = self.add_units(rhs_units.simplify(indep_var.get_units().extract()))
                        dep_var.units = units.name
                        magic_vars.add(dep_var)
        for var in self.model.get_all_variables():
            if var.get_units() is self.magic_units:
                src = var.get_source_variable(recurse=True)
                if src in magic_vars:
#                     print 'Magic connection', var, src, var.units, src.units
                    var.units = src.units
    
    def _check_equation_lhs(self, expr):
        """Check whether the variable on the LHS of a new equation exists, or whether we can create it.
        
        There are two cases in which the protocol can create an annotated variable not existing in the model:
        an optional variable with a default definition, or an output with a separate equation overriding any definition.
        If the latter there will also be units specified in the output definition, but in the former case we need to use
        'magic' units until variable connections have been sorted out, and we can calculate the real units.
        """
        assert isinstance(expr, mathml_apply)
        assert expr.operator().localName == u'eq', 'Expression is not an assignment'
        lhs, rhs = list(expr.operands())
        if lhs.localName == u'ci' and u':' in unicode(lhs):
            try:
                vname = unicode(lhs)
                var = self._lookup_ontology_term(vname)
            except ValueError:
                # Is this declared as an output with units?
                for output_spec in self._output_specifications:
                    if output_spec['prefixed_name'] == vname and output_spec['units']:
                        var = self._create_annotated_variable(vname, output_spec['units'])
#                         print 'Created output', vname, '=', var, 'with units', var.units
                        return
                # Is this an optional variable (with this expression being the default definition)?
                if vname in self._optional_vars:
#                     print 'Creating', vname, 'with magic units'
                    var = self._create_annotated_variable(vname, self.magic_units)
                else:
                    raise ProtocolError("Variable %s on the LHS of a 'define' does not exist in the model and is not an output with units or optional." % vname)
    
    def report_stats(self):
        """Output a short report on what the modified model looks like.
        
        Write an output file listing the key features of the modified model: what were requested as inputs & outputs,
        what state variables resulted, how many other equations were computed.  Counts of variables and equations are
        also written to standard error.
        """
        inputs = [input for input in self.inputs if isinstance(input, cellml_variable)]
        num_outputs = len(self.outputs | self._vector_outputs)
        num_equations = len([e for e in self.model._cml_assignments if isinstance(e, mathml_apply)])
        state_vars = self.model.find_state_vars()
        all_vars = list(self.model.get_all_variables())
        print >>sys.stderr, 'Statistics about the modified model:'
        print >>sys.stderr, '    # inputs:         ', len(inputs)#, map(str, inputs)
        print >>sys.stderr, '    # total outputs:  ', num_outputs#, map(str, self.outputs)
        print >>sys.stderr, '    # state variables:', len(state_vars)
        print >>sys.stderr, '    # equations:      ', num_equations
        print >>sys.stderr, '    # variables:      ', len(all_vars)
        missing_vars = []
        for prefixed_name in sorted(self._optional_vars):
            if self._lookup_ontology_term(prefixed_name, check_optional=True) is None:
                missing_vars.append(prefixed_name)
        if missing_vars:
            print >>sys.stderr, '    Missing optional variables:'
            for prefixed_name in missing_vars:
                print >>sys.stderr, ' '*7, prefixed_name
    
    def add_alias(self, var, alias):
        """Add an alias name for a variable.
        
        This is used by the SED-ML support to map XPath expressions to the names in the generated code.
        """
        var.add_rdf_annotation(('pycml:alias', NSS['pycml']), alias, allow_dup=True)
    
    def specify_as_output(self, var, units):
        """Specify the given variable within the model as a protocol output.
        
        The output is wanted in the given units, which must be added to the model if they don't exist.
        If they differ from its current units, a conversion will be needed, and hence a new version
        of this variable will be added to the protocol component, with a suitable connection.
        Otherwise, we can just record the existing variable as an output.
        
        TODO: We can't (yet?) specify a variable as both an output and an input but in different units.
        """
        units = self._get_units_object(units)
        if units.equals(var.get_units()):
            output_var = var
        else:
            if var in self.inputs:
                raise ProtocolError("You can't specify a variable (%s) as output and input in different units!" % var.fullname())
            output_var = self._replace_variable(var, units)
            self.connect_variables(var, output_var)
        self.outputs.add(output_var)
        return output_var
    
    def specify_as_input(self, var, units, copy_initial_value=True):
        """Specify the given variable within the model as a protocol input.
        
        The input is wanted in the given units, which must be added to the model if they don't exist.
        If they differ from its current units, a conversion will be needed, and hence a new version
        of this variable will be added to the protocol component, with a suitable assignment.
        
        The variable that is the input must be set as a modifiable parameter, and any existing definition
        removed.
        
        Note: this method is now only used by legacy low-level tests, not the new language support.
        TODO: remove it!
        """
        if var.get_type() == VarTypes.Mapped:
            raise ProtocolError("Cannot specify a mapped variable (%s) as an input." % var.fullname())
        # Remove any existing definition
        self.remove_definition(var, keep_initial_value=True)
        # Set up the input
        units = self._get_units_object(units)
        if units.equals(var.get_units()):
            input_var = var
            if not hasattr(var, u'initial_value'):
                var.initial_value = u'0'
        else:
            input_var = self._replace_variable(var, units)
            if copy_initial_value:
                if not hasattr(var, u'initial_value'):
                    raise ProtocolError("No initial value available for input " + str(var))
                input_var.initial_value = unicode(self._convert_initial_value(var, units))
            self.del_attr(var, u'initial_value', None)
            # Set all variables connected to the original variable (including itself) to be mapped to the new one
            self._update_connections(var, input_var)
        input_var._set_type(VarTypes.Constant)
        input_var._cml_ok_as_input = True
        self.inputs.add(input_var)
        return input_var

    def _add_converted_variable(self, orig_var, new_units):
        """Add a new version of the given variable with different units.
        
        An assignment will be added making the new variable equal to the old, so that a units conversion will happen.
        We also ensure that the new units are added to the model if needed.
        Used by _add_interpolation, and relies on being in the _add_maths_to_model processing phase.
        
        :returns: the new variable
        """
        new_units = self.add_units(new_units)
        comp = orig_var.component
        new_name = self._uniquify_var_name(orig_var.name, comp)
        new_var = self.add_variable(comp, new_name, new_units)
        cname = comp.name
        assign = mathml_apply.create_new(new_var, 'eq', [cname + ',' + new_var.name, cname + ',' + orig_var.name])
        self.inputs.add(assign)
        self._add_maths_to_model(assign)
        return new_var

    def _replace_variable(self, var, units, allow_existing=False):
        """Replace the given variable with a version in the given units in the protocol component.

        Ensures that the units are added to the model if needed, and transfers the cmeta:id if
        present.  It doesn't transfer the initial_value, since this would break the output variable
        case.

        The new variable will be given a local name equal to the full name of the original, to avoid
        potential conflicts.  If allow_existing is False then it's an error if the variable has
        already been replaced.  If allow_existing is True, then we just reuse the existing replacement.
        """
        units = self.add_units(units)
        new_name = var.fullname(cellml=True)
        comp = self._get_protocol_component()
        try:
            existing_replacement = comp.get_variable_by_name(new_name)
        except KeyError:
            existing_replacement = None
        if existing_replacement:
            if not allow_existing:
                raise ProtocolError("Variable '%s' has already been replaced!" % new_name)
            new_var = existing_replacement
        else:
            new_var = self.add_variable(comp, new_name, units, id=var.cmeta_id)
            self.del_attr(var, u'id', NSS['cmeta'])
#         print 'Replacing', var, 'by', new_var, 'with id', new_var.cmeta_id
        return new_var
    
    def _check_input(self, input):
        """New inputs must not already exist in the model!"""
        if isinstance(input, cellml_units):
            exists = self.model.has_units(input)
        else:
            exists = self.model is getattr(input, 'xml_parent', None)
        if exists and not getattr(input, '_cml_ok_as_input', False):
            raise ProtocolError("Inputs must not already exist in the model. (Input %s exists.)" % repr(input))
        
    def _error_handler(self, errors):
        """Deal with errors found when re-analysing a modified model."""
        raise ProtocolError('Applying protocol created an invalid model:\n  ' + '\n  '.join(map(str, errors)))
    
    def _check_if_output(self, old_var, new_var):
        """A variable is being replaced.  If the original was an output, make the new one instead."""
        if old_var in self.outputs:
            self.outputs.remove(old_var)
            self.outputs.add(new_var)
        if old_var in self._vector_outputs:
            self._vector_outputs.remove(old_var)
            self._vector_outputs.add(new_var)
    
    def _split_all_odes(self):
        """The free variable has been units-converted, so adjust all ODEs to account for this.
        
        We copy all state variables into the protocol component, assign their RHS to a new variable,
        and create a new derivative in the protocol component to which this is assigned.
        
        Also check all equations for occurrences of derivatives on the RHS, and change them to refer
        to the new variable instead.
        
        TODO: How to deal with ODEs added by the protocol?  Will they be picked up?
        TODO: Require explicit call to set free var units even if an output?
        """
        free_var = self._free_var_has_changed
        old_free_var = self.model.find_free_vars()[0]
        self._update_connections(old_free_var, free_var)
        deriv_rhs = {}
        comp = self._get_protocol_component()
        for old_var in self.model.find_state_vars():
            if old_var.component is not comp:
                new_var = self._replace_variable(old_var, old_var.get_units(), allow_existing=True)
                new_var.initial_value = old_var.initial_value
                del old_var.initial_value
                self._check_if_output(old_var, new_var)
                # Add a new variable to assign the RHS to, with units of the original derivative
                deriv_name = self._uniquify_var_name(u'd_%s_d_%s' % (old_var.name, free_var.name), old_var.component)
                orig_ode = old_var.get_all_expr_dependencies()[0]
                orig_rhs_var = self.add_variable(old_var.component, deriv_name, orig_ode.eq.lhs.get_units().extract())
                deriv_rhs[new_var] = orig_rhs_var
                # Add a version of this in the protocol component, with desired units
                desired_units = new_var.get_units().quotient(free_var.get_units())
                mapped_rhs_var = self._replace_variable(orig_rhs_var, desired_units)
                self.connect_variables(orig_rhs_var, mapped_rhs_var)
                # Replace the original ODE with an assignment
                orig_rhs = orig_ode.eq.rhs
                orig_ode.safe_remove_child(orig_rhs)
                self.remove_expr(orig_ode)
                self.add_expr_to_comp(old_var.component,
                                      mathml_apply.create_new(self.model, u'eq',
                                                              [orig_rhs_var.name, orig_rhs]))
                # Create a new ODE in the interface component
                new_ode = mathml_diff.create_new(self.model, free_var.name, new_var.name, mapped_rhs_var.name)
                self.add_expr_to_comp(new_var.component, new_ode)
                new_ode.classify_variables(root=True, dependencies_only=True)
                # Update connections to the state variable
                self._update_connections(old_var, new_var)
            else:
                print >>sys.stderr, "Ignoring state var", old_var
                raise NotImplementedError # TODO
        # Transform references to derivatives
        def xform(expr):
            state_var = expr.diff.dependent_variable.get_source_variable(recurse=True)
            rhs_var = deriv_rhs[state_var]
            # Ensure there's something mapped to it in this component
            rhs_var = self.connect_variables(rhs_var, (expr.component.name, rhs_var.name))
            # Update this expression
            parent = expr.xml_parent
            parent.xml_insert_after(expr, mathml_ci.create_new(parent, rhs_var.name))
            parent.safe_remove_child(expr)
        for expr in self.model.search_for_assignments():
            self._process_operator(list(expr.operands())[1], u'diff', xform)

    def _fix_model_connections(self):
        """Ensure the modified model has all the necessary connections between variables.
        
        Check mathematics for ci elements that refer to variables not defined in that
        component.  These must refer to the variable by its full name (i.e. 'cname,vname').
        These variables will be renamed to use local names by this method, which will
        also create local variables mapped to the relevant source variable if needed.
        
        This needs to take account of the fact that a variable in one nested component
        may need to be connected to a variable in another nested component, and so create
        variables in the parent components to connect the whole thing up.
        """
        for expr in self.model.search_for_assignments():
            for ci_elt in self._find_ci_elts(expr):
                vname = unicode(ci_elt)
                if u',' in vname:
                    cname, vname = self._split_name(vname)
                    comp = expr.component
                    if comp.name != cname:
                        # Check for the special case of the referenced variable having a source in this component already
                        # (ensuring either it has identical units or the expression will be units converted)
                        src_comp = self.model.get_component_by_name(cname)
                        referenced_var = src_comp.get_variable_by_name(vname)
                        src_var = referenced_var.get_source_variable(recurse=True)
                        if src_var.component is comp:
                            # Use the existing var
#                             print 'Using existing var', src_var, 'for reference', unicode(ci_elt), 'in', comp.name
                            if not src_var.get_units().equals(referenced_var.get_units()):
                                self.assignments_to_convert.add(expr)
                            vname = src_var.name
                        else:
#                             print 'Connecting to reference', unicode(ci_elt), 'from component', comp.name
                            self.connect_variables(referenced_var, (comp.name,vname))
                    # Now just rename to be local
                    ci_elt._rename(vname)
    
    def _get_protocol_component(self):
        """Get the protocol component in the model, creating it if necessary.
        
        New variables created just for use by the simulation protocol get put into a
        new 'protocol' component.  If a component with that name already exists,
        underscores will be added to the component name to make it unique.
        """
        if self._protocol_component is None:
            self._protocol_component = self.create_new_component(u'protocol')
        return self._protocol_component
    
    def _split_name(self, full_name):
        """Split a full name into cname,vname, creating the protocol component if needed.
        
        If the full_name doesn't contain a component part, the 'protocol' component
        will be used.
        """
        parts = full_name.split(',')
        if len(parts) == 2:
            cname, vname = parts
        elif len(parts) == 1:
            cname = self._get_protocol_component().name
            vname = full_name
        else:
            raise ValueError("Invalid variable name: " + full_name)
        return cname, vname
        
    def _find_ci_elts(self, expr):
        """Get an iterator over all ci elements on the descendent-or-self axis of the given element."""
        if isinstance(expr, mathml_ci):
            yield expr
        elif hasattr(expr, 'xml_children'):
            # Recurse
            for child in expr.xml_children:
                for ci_elt in self._find_ci_elts(child):
                    yield ci_elt
    
    def _identify_referenced_variables(self, expr, special_name=None, check_optional=False):
        """Figure out which variables are referenced in the given expression, and update ci elements.
        
        The expression should contain names as used in the protocol, i.e. prefixed names giving an
        ontology name for model variables, and bare names for variables added by the protocol.
        Change each ci element to use the full "compname,varname" format.
        A ValueError is raised if any referenced variable doesn't exist.
        However, any reference to special_name is not checked and left as-is.
        Also, names already in fully qualified form are assumed to be ok and left as-is.
        
        If check_optional is True, we don't throw on missing ontology-annotated variables if they're
        in the optional set, but just return False instead.
        """
        all_vars_found = True
        for ci_elt in self._find_ci_elts(expr):
            vname = unicode(ci_elt)
            if vname == special_name:
                continue
            if ',' in vname:
                continue
            if ':' in vname:
                var = self._lookup_ontology_term(vname, check_optional=check_optional)
                if var is None:
                    # It was a missing optional variable
                    all_vars_found = False
                    continue
            else:
                try:
                    var = self._get_protocol_component().get_variable_by_name(vname)
                except KeyError:
                    raise ValueError("The variable name '%s' has not been declared in the protocol" % vname)
            full_name = var.component.name + u',' + var.name
            ci_elt._rename(full_name)
        return all_vars_found
    
    def _lookup_ontology_term(self, prefixed_name, enforce_uniqueness=True, check_optional=False, transitive=False):
        """Find the variable annotated with the given term, if it exists.
        
        The term should be given in prefixed form, with the prefix appearing in the protocol's namespace
        mapping (prefix->uri, as found e.g. at elt.rootNode.xmlns_prefixes).
        
        Will throw ValueError if the variable doesn't exist in the model, or the given term is invalid.
        :param enforce_uniqueness: if True, also ensures there's only one variable with the annotation.
        :param check_optional: if True, we don't throw on missing ontology-annotated variables if they're
        in the optional set, but just return None instead.
        :param transitive: if True, look not just for direct annotations but also for terms belonging to
        the class given by prefixed_name, searching transitively along rdf:type predicates.
        """
        try:
            prefix, _ = prefixed_name.split(':')
        except ValueError:
            raise ValueError("Ontology term '%s' is not a qname - it doesn't have a namespace prefix"
                             % prefixed_name)
        try:
            nsuri = self._protocol_namespaces[prefix]
        except KeyError:
            raise ValueError("The namespace prefix '%s' has not been declared" % prefix)
        vars = self.model.get_variables_by_ontology_term((prefixed_name, nsuri), transitive=transitive)
        if len(vars) == 0:
            if check_optional and prefixed_name in self._optional_vars:
                return None
            else:
                raise ValueError("The ontology term '%s' does not match any variables" % prefixed_name)
        if enforce_uniqueness:
            if len(vars) > 1:
                raise ValueError("The ontology term '%s' matches multiple variables" % prefixed_name)
            else:
                vars = vars[0]
#         print 'Looked up', prefixed_name, 'as', vars
        return vars
    
    def _apply_conversion_rule(self, expr, conv_template, placeholder_name):
        """Apply a units conversion rule defined by self.add_units_conversion_rule.
        
        Modify the given expr in-place, replacing the RHS/expr itself by a copy of conv_template, except
        ci references to placeholder_name are replaced by (a copy of) the original RHS/expr.
        If the expression is a top-level assignment, only the RHS is modified.
        Returns the modified RHS/expression.
        """
#         print '_apply_conv_rule to', element_xpath(expr), 'top-level =', isinstance(expr, mathml_apply) and expr.is_top_level()
        if isinstance(expr, mathml_apply) and expr.is_top_level():
            expr = expr.eq.rhs
        parent = expr.xml_parent
        parent.safe_remove_child(expr)
        new_expr = mathml.clone(conv_template)
        copy_expr = False
        for ci_elt in self._find_ci_elts(new_expr):
            vname = unicode(ci_elt).strip()
            if vname == placeholder_name:
                # Copy the original RHS here, except if it's the first use don't bother copying
                if copy_expr:
                    expr = mathml.clone(expr)
                else:
                    copy_expr = True
                ci_elt.xml_parent.replace_child(ci_elt, expr)
            else:
                # Ensure we have connections needed to get the variable in this component
                cname, local_name = vname.split(',')
                our_cname = parent.component.name
                if cname != our_cname:
                    local_var = self.connect_variables((cname, local_name), (our_cname, local_name))
                    local_name = local_var.name
                ci_elt._rename(local_name)
        parent.xml_append(new_expr)
        return new_expr

    def _add_units_conversions(self):
        """Apply units conversions, in particular 'special' ones, to the protocol component.

        We first check for any variables with 'magic' units and try to give them real units based
        on their definitions.

        Also convert any equations that have been added to the model, even if they don't appear
        in the protocol component, since otherwise we might not do all necessary conversions
        between model & protocol mathematics.
        """
        converter = self.get_units_converter()
        notifier = NotifyHandler(level=logging.WARNING)
        logging.getLogger('units-converter').addHandler(notifier)
        self._fix_magic_units()
        proto_comp = self._get_protocol_component()
        converter.add_conversions_for_component(proto_comp)
        self.assignments_to_convert.update(filter(lambda i: isinstance(i, mathml_apply), self.inputs))
        converter.convert_assignments(self.assignments_to_convert)
        converter.convert_connections(self.connections_made)
        converter.finalize(self._error_handler, check_units=False)
        notifier.flush()
        logging.getLogger('units-converter').removeHandler(notifier)
        if notifier.messages:
            raise ProtocolError("Unable to apply units conversions to the model/protocol interface")
        
    def _add_variable_to_model(self, var):
        """Add or replace a variable in our model.
        
        We don't really do any checking for this case - just add the variable.
        This means that some 'possible' changes don't actually make sense, for
        instance giving a computed variable an initial value will trigger a later
        validation error, unless its definition is also changed to an ODE.
        (To change a variable to a constant, you need to replace its definition
        with a constant expression.)
        
        The one check that is made is that you don't change the interface
        definitions if replacing a variable, otherwise self._fix_model_connections
        could get confused.
        """
        if hasattr(var, 'xml_parent'):
            # It's already been added, e.g. by specify_as_input
            return
        cname, vname = self._split_name(var.name)
        comp = self.model.get_component_by_name(cname)
        try:
            orig_var = comp.get_variable_by_name(vname)
        except KeyError:
            orig_var = None
        if orig_var:
            # We're replacing a variable
            for iface in [u'public', u'private']:
                n = iface + u'_interface'
                assert getattr(orig_var, n, u'none') == getattr(var, n, u'none'), "You are not allowed to change a variable's interfaces"
            # Only keep RDF annotations if the cmeta:id is unchanged
            comp._del_variable(orig_var, keep_annotations=(orig_var.cmeta_id == var.cmeta_id))
        var.name = vname
        comp._add_variable(var)
    
    def _force_evaluate(self, variable):
        """Attempt to force the evaluation of variable, even if it isn't static.
        
        We do a recursive sweep of the variable's dependencies in order to set each variable as temporarily static,
        so that we can use var.get_value() to do the evaluation.  Otherwise the first variable lookup with var's
        definition would give an error if the definition isn't known to be static.
        """
        def process_defn(var, set=True):
            """Set or unset var and its dependencies as static."""
            if set:
                var._set_binding_time(BINDING_TIMES.static, temporary=True)
                # If this is an input with value specified, also set the value temporarily
#                 print 'Checking', var, var.oxmeta_name
                for spec in self._input_specifications:
                    try:
                        if var is self._lookup_ontology_term(spec['prefixed_name']):
#                             print 'Set value for', var, 'using', spec
                            if spec['initial_value']:
                                var.set_value(spec['initial_value'], follow_maps=False)
                            break
                    except ValueError:
                        pass
            else:
                var._unset_binding_time(only_temporary=True)
                var.unset_values()
            defn = var.get_dependencies()
            if defn:
                if isinstance(defn[0], mathml_apply):
                    for ci_elt in self._find_ci_elts(defn[0].eq.rhs):
                        process_defn(ci_elt.variable)
                elif isinstance(defn[0], cellml_variable):
                    process_defn(defn[0])
        try:
            process_defn(variable, set=True)
            value = unicode("%.17g" % variable.get_value())
        finally:
            process_defn(variable, set=False)
#         print 'Evaluated', variable, 'to', value
        return value

    def _add_interpolation(self, expr, lhs, rhs):
        """Add/replace a model equation with linear interpolation from a data file.

        This is a special case of the 'define' construct.  The rhs is an apply of the
        interpolate csymbol with 4 arguments: a data file path (relative to the protocol),
        an independent variable reference, units of the independent variable in the data,
        and units of the lhs values in the data.  The data file contains two columns
        of values for the independent variable and lhs, respectively.  We construct a new
        piecewise definition for the lhs using linear interpolation on these values,
        modifying the rhs of expr to contain this piecewise MathML, so that _add_maths_to_model
        can continue as for normal 'define' constructs.

        :param expr: the new equation
        :param lhs: the left-hand side of the new equation
        :param rhs: the right-hand side of the new equation 
        """
        # Extract operands with paranoia checks
        operands = list(rhs.operands())
        assert len(operands) == 4
        assert operands[0].localName == u'csymbol' and getattr(operands[0], u'definitionURL', u'') == u'https://chaste.cs.ox.ac.uk/nss/protocol/string'
        for i in [1, 2, 3]:
            assert isinstance(operands[i], mathml_ci)
        data_path, indep_var_name, indep_units, dep_units = map(unicode, operands)
        # Load the data
        data_path = os.path.join(self.base, data_path)
        if not os.path.exists(data_path):
            raise ProtocolError("Unable to load data file '%s' for interpolation." % data_path)
        import numpy
        data = numpy.loadtxt(data_path, dtype=float, delimiter=',', ndmin=2, unpack=True)
        assert data.ndim == 2
        if data.shape[0] != 2:
            raise ProtocolError("The data file for an interpolate() must have 2 columns; file '%s' has %d." % (data_path, data.shape[0]))
        # Get units objects so we can do conversions if necessary
        indep_units = self._get_units_object(self.units[indep_units])
        dep_units = self._get_units_object(self.units[dep_units])
        # Check to see whether we can do a special-case optimisation for equally spaced independent variable points
        steps = numpy.diff(data[0])
        if abs(numpy.max(steps) - numpy.min(steps)) < 1e-6 * numpy.max(steps):
            if steps[0] == 0.0:
                raise ProtocolError("The data file for an interpolate() cannot have all zero table steps.")
            # Find the independent variable, which must exist, to determine if we need a conversion
            self._identify_referenced_variables(operands[1])
            indep_var_name = unicode(operands[1])
            indep_var = self.model.get_variable_by_name(*self._split_name(indep_var_name))
            if not indep_var.get_units().equals(indep_units):
                indep_var = self._add_converted_variable(indep_var, indep_units)
                indep_var_name = indep_var.component.name + ',' + indep_var.name
            # Create a 'magic' expression that code generation will comprehend
            new_rhs = mathml_apply.create_new(rhs, 'csymbol', [indep_var_name])
            new_rhs._cml_interp_data = data
            new_rhs._cml_units = UnitsSet([dep_units], expression=new_rhs)
            new_rhs.xml_set_attribute((u'pe:binding_time', NSS[u'pe']), u'dynamic')
            new_rhs.operands().next().xml_set_attribute((u'pe:binding_time', NSS[u'pe']), u'dynamic')
            new_rhs.operator().xml_set_attribute(u'definitionURL', u'https://chaste.cs.ox.ac.uk/nss/protocol/interp')
            # Update list of interpolations on the model for convenience
            self.model._cml_interp_exprs = getattr(self.model, '_cml_interp_exprs', [])
            self.model._cml_interp_exprs.append(new_rhs)
        else:
            # Iterate over data rows to construct a piecewise representation for the interpolation
            eqn = mathml_apply.create_new
            pieces = []
            for row in range(data.shape[1] - 1):
                if data[0,row] == data[0,row+1]:
                    # Vertical jumps are allowed in voltage clamp protocols, but trying to interpolate gives you a divide by zero.
                    # Just skip the jump part, and we'll interpolate correctly at either side.
                    continue
                x1 = (data[0,row], indep_units.name)
                x2 = (data[0,row+1], indep_units.name)
                y1 = (data[1,row], dep_units.name)
                y2 = (data[1,row+1], dep_units.name)
                cond = eqn(rhs, 'and', [eqn(rhs, 'geq', [indep_var_name, x1]),
                                        eqn(rhs, 'leq', [indep_var_name, x2])])
                temp2 = eqn(rhs, 'divide', [eqn(rhs, 'minus', [y2, y1]),
                                            eqn(rhs, 'minus', [x2, x1])])
                temp1 = eqn(rhs, 'times', [temp2, eqn(rhs, 'minus', [indep_var_name, x1])])
                case = eqn(rhs, 'plus', [y1, temp1])
                pieces.append((case, cond))
            new_rhs = mathml_piecewise.create_new(rhs, pieces)
        # Insert and return the new rhs
        expr.replace_child(rhs, new_rhs)
        return new_rhs

    def _get_cn_units(self, elt):
        """Get the units name for a cn element, dereferencing a units_of() construct if present.
        
        Returns None if units_of() references a non-existent optional variable.
        """
        if isinstance(elt, mathml_apply):
            # This is minus a cn; negative constants parse awkwardly due to precedence!
            elt = elt.xml_children[-1]
        assert isinstance(elt, mathml_cn)
        if not hasattr(elt, u'units'):
            raise ProtocolError("All numbers in model equations must have units; '%s' does not." % unicode(elt))
        units = elt.units
        if units.startswith('units_of('):
            # Figure out the units of the variable referred to
            vname = units[9:-1]
            assert ':' in vname
            var = self._lookup_ontology_term(vname, check_optional=True)
            if var is None:
                units = None
            else:
                print 'units_of(): using', var.units, 'from', var, 'for', elt.units
                units = var.units
        return units

    def _add_maths_to_model(self, expr):
        """Add or replace an equation in the model.
        
        This case is more complex than variables, since we may need to change
        the type of variable assigned to, depending on the expression.
        
        Note: variable references within ci elements in the given expression
        should use full names (i.e. cname,vname) or ontology terms.  Any local names will be
        assumed to refer to variables in the protocol component, and modified
        by self._rename_local_variables.  Later, self._fix_model_connections
        will change all references to use local names.
        """
        assert isinstance(expr, mathml_apply)
        assert expr.operator().localName == u'eq', 'Expression is not an assignment'
        lhs, rhs = list(expr.operands())
        if (isinstance(rhs, mathml_apply) and rhs.operator().localName == u'csymbol' and
            getattr(rhs.operator(), u'definitionURL', u'') == u'https://chaste.cs.ox.ac.uk/nss/protocol/interpolate'):
            # Special case for linear interpolation on a data file
            rhs = self._add_interpolation(expr, lhs, rhs)
        orig_defn_kw = u'original_definition' # References the original definition if it appears as a variable on the new RHS
        def handle_unusable_expr():
            """Some optional variables weren't present in the model, so we can't use this equation."""
            print >>sys.stderr, "Warning: optional model variables missing, so not using model interface equation:"
            if hasattr(expr, 'loc'):
                print >>sys.stderr, "  ", expr.loc
            # Check, however, whether it defines a (non-optional) output or the default value for an optional variable, since this is then an error.
            if lhs.localName == u'ci':
                vname = unicode(lhs)
                if u':' in vname:
                    if vname in self._optional_vars:
                        raise ProtocolError("Cannot give a value to optional variable '%s' as not all variables used in its default clause are present." % vname)
                    for output_spec in self._output_specifications:
                        if output_spec['prefixed_name'] == vname and not output_spec['optional']:
                            raise ProtocolError("At least one optional variable required to override the definition of model output '%s' was not found and has no default value." % vname)
            self.inputs.remove(expr)
        if not self._identify_referenced_variables(expr, check_optional=True, special_name=orig_defn_kw):
            return handle_unusable_expr()
        # Determine if the original_definition keyword occurs on the RHS, and check that all numbers have units
        orig_defn_refs = []
        def find_refs(elt):
            if isinstance(elt, mathml_ci):
                if unicode(elt) == orig_defn_kw:
                    orig_defn_refs.append(elt)
            elif isinstance(elt, mathml_cn):
                units = self._get_cn_units(elt)
                if units is None:
                    return handle_unusable_expr()
                elt.units = units
            else:
                for child in elt.xml_element_children():
                    find_refs(child)
        find_refs(lhs)
        if orig_defn_refs:
            raise ProtocolError("Cannot assign to the %s keyword." % orig_defn_kw)
        find_refs(rhs)
        def get_orig_defn(var):
            """Find and store the original definition before removing it."""
            orig_defn = var.get_all_expr_dependencies()
            if orig_defn:
                assert len(orig_defn) == 1
                orig_defn = list(orig_defn[0].operands())[1] # Take the RHS of the defining ODE
            else:
                # It had better be a constant originally with initial_value
                orig_defn = mathml_cn.create_new(var, var.initial_value, var.units)
            return orig_defn
        # Figure out what's on the LHS of the assignment
        if lhs.localName == u'ci':
            # Straight assignment to variable
            cname, vname = self._split_name(unicode(lhs))
            assigned_var = self.model.get_variable_by_name(cname, vname)
            # Check for the special case of "var = var" which signifies clamping to initial value
            clamping = isinstance(rhs, mathml_ci) and unicode(rhs) == unicode(lhs)
            if clamping and not hasattr(assigned_var, u'initial_value'):
                # Hope it's computed, and try to evaluate the definition
                try:
#                     print 'Computing initial value for clamped', unicode(lhs), assigned_var
                    assigned_var.initial_value = self._force_evaluate(assigned_var)
                except:
                    pass # TODO: Check against optional variables in order to give nicer error message
#                     import traceback
#                     traceback.print_exc()
#                     raise ProtocolError("No suitable value found for clamped variable " + unicode(lhs))
            if orig_defn_refs:
                orig_defn = get_orig_defn(assigned_var)
            self.remove_definition(assigned_var, keep_initial_value=clamping)
            if clamping:
                self.inputs.remove(expr) # The equation isn't actually used in this case
            else:
                self.add_expr_to_comp(cname, expr)
                assigned_var._add_dependency(expr)
        else:
            # This had better be an ODE
            assert lhs.localName == u'apply', 'Expression is not a straight assignment or ODE'
            assert lhs.operator().localName == u'diff', 'Expression is not a straight assignment or ODE'
            dep_var = lhs.operands().next()
            assert dep_var.localName == u'ci', 'ODE is malformed'
            cname, dep_var_name = self._split_name(unicode(dep_var))
            dep_var = self.model.get_variable_by_name(cname, dep_var_name)
            if orig_defn_refs:
                orig_defn = get_orig_defn(dep_var)
            self.remove_definition(dep_var, keep_initial_value=True)
            self.add_expr_to_comp(cname, expr)
            indep_name = self._split_name(unicode(lhs.bvar.ci))
            indep_var = self.model.get_variable_by_name(*indep_name)
            dep_var._add_ode_dependency(indep_var, expr)
        # Insert the original definition in the new RHS if needed
        for ref in orig_defn_refs:
            ref.xml_parent.replace_child(ref, orig_defn.clone_self())

    def _filter_assignments(self):
        """Apply protocol outputs to reduce the model size.
        
        Protocol outputs are a list of variable objects that are of interest.
        The assignments used in computing the model should be filtered so that
        only those needed for determining the outputs are used.  This has the
        potential to greatly simplify the model simulation.
        
        The only assignments should be to output variables, or nodes required
        in computing these.  Note that if one of these nodes is a state variable,
        we also require its derivative and the dependencies thereof.
        If there are no outputs specified, we leave the list unchanged.
        
        In addition, any output computed variable should be annotated as a 
        derived quantity, and any output constant annotated as a parameter, to
        ensure they are available for querying.  Other variables should have these
        annotations (and pe:keep) removed.
        """
        all_outputs = self.outputs | self._vector_outputs
        if all_outputs:
            # Remove parts of the model that aren't needed
            needed_nodes = all_outputs.copy()
            needed_nodes.update([input for input in self.inputs
                                 if isinstance(input, cellml_variable)])
            needed_nodes = self.model.calculate_extended_dependencies(needed_nodes,
                                                                      state_vars_depend_on_odes=True)
            for node in self.model.get_assignments()[:]:
                if node not in needed_nodes:
                    if isinstance(node, cellml_variable):
                        node.component._del_variable(node)
                        try:
                            self.inputs.remove(node)
#                             print 'Removed', node, 'from inputs'
                        except:
                            pass
                    elif isinstance(node, mathml_apply):
                        node.xml_parent.xml_remove_child(node)
            # Update connection elements
            for conn in list(getattr(self.model, u'connection', [])):
                comp1 = self.model.get_component_by_name(conn.map_components.component_1)
                comp2 = self.model.get_component_by_name(conn.map_components.component_2)
                any_kept = False
                for mapv in list(conn.map_variables):
                    try:
                        comp1.get_variable_by_name(mapv.variable_1)
                        comp2.get_variable_by_name(mapv.variable_2)
                        any_kept = True
                    except KeyError:
                        # Remove connection
                        conn.xml_remove_child(mapv)
                if not any_kept:
                    self.model.xml_remove_child(conn)
            # Filter assignments list
            new_assignments = filter(lambda node: node in needed_nodes,
                                     self.model.get_assignments())
            self.model._cml_assignments = new_assignments
        # Remove existing annotations
        free_vars, unknown_vars = [], []
        for var in self.model.get_all_variables():
            var.set_pe_keep(False)
            var.set_is_derived_quantity(False)
            var.set_is_modifiable_parameter(False)
            # An algebraic model might have no free vars left, but if there is a single unknown var then treat it as free
            if var.get_type() == VarTypes.Free:
                free_vars.append(var)
            if var.get_type() == VarTypes.Unknown:
                unknown_vars.append(var)
        if len(free_vars) == 0 and len(unknown_vars) == 1:
            unknown_vars[0]._set_type(VarTypes.Free)
        # Add annotations for inputs & outputs
        for var in [input for input in self.inputs if isinstance(input, cellml_variable)]:
            var.set_pe_keep(True)
            if var.get_type() == VarTypes.Constant:
                var.set_is_modifiable_parameter(True)
            var.set_rdf_annotation_from_boolean(('pycml:input-variable', NSS['pycml']), True)
        for var in all_outputs:
            assert isinstance(var, cellml_variable)
            if var.get_type() == VarTypes.Constant:
                var.set_is_modifiable_parameter(True)
            elif var.get_type() in [VarTypes.Computed, VarTypes.Mapped]:
                var.set_is_derived_quantity(True)
            else:
                assert var.get_type() in [VarTypes.State, VarTypes.Free], 'Bad var ' + str(var) + ' of type ' + str(var.get_type())
        for var in self.outputs:
            var.set_is_output_variable(True)
        # Make constant variables with name annotations into model parameters if requested
        if self.model.get_option('expose_named_parameters'):
            for var in cellml_metadata.find_variables(self.model, ('bqbiol:is', NSS['bqbiol'])):
                if var.get_type() == VarTypes.Constant:
                    var.set_is_modifiable_parameter(True)

def apply_protocol_file(doc, proto_file_path):
    """Apply the protocol defined in the given file to a model.
    
    New protocols should be written in the pure XML syntax, for which we use
    Protocol.apply_protocol_file.  However, legacy protocols may be Python code
    with a method apply_protocol(doc) to do the donkey work.
    """
    if proto_file_path[-3:] == '.py':
        import imp
        proto_dir = os.path.dirname(proto_file_path)
        proto_file_name = os.path.basename(proto_file_path)
        proto_module_name = os.path.splitext(proto_file_name)[0]
        (file, pathname, desc) = imp.find_module(proto_module_name, [proto_dir])
        try:
            proto = imp.load_module(proto_module_name, file, pathname, desc)
        finally:
            file.close()
        proto.apply_protocol(doc)
    elif proto_file_path[-4:] == '.xml':
        Protocol.apply_protocol_file(doc, proto_file_path)
    else:
        raise ProtocolError("Unexpected protocol file extension for file: " + proto_file_path)

if __name__ == '__main__':
    # Analyse the supplied protocol file to determine model annotations required
    if len(sys.argv) < 2:
        print "Usage:", sys.argv[0], "<protocol.xml>"
        sys.exit(1)
    Protocol.find_required_annotations(sys.argv[1])
