
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
This module abstracts the interface to RDF metadata about CellML models.
"""

import logging
import types
from cStringIO import StringIO

# We support 2 RDF libraries
RDF = rdflib = None
try:
    import rdflib
except ImportError:
    import RDF

import pycml

# Allowed metadata names, more to come
# TODO #1209: Use a proper ontology!
METADATA_NAMES = frozenset(
    ['membrane_voltage', 'time', 'membrane_capacitance', 'membrane_E_R', 'membrane_stimulus_current', 'membrane_stimulus_current_duration',
     'membrane_stimulus_current_amplitude','membrane_stimulus_current_period','membrane_stimulus_current_offset','membrane_stimulus_current_end',
     'membrane_fast_sodium_current', 'membrane_fast_sodium_current_conductance', 'membrane_fast_sodium_current_m_gate', 'membrane_fast_sodium_current_h_gate',
     'membrane_potassium_current', 'membrane_potassium_current_conductance', 'potassium_channel_n_gate', 
     'leakage_current','membrane_fast_sodium_current_j_gate', 'temperature',
     'membrane_rapid_delayed_rectifier_potassium_current_conductance','membrane_slow_delayed_rectifier_potassium_current_conductance',
     'membrane_inward_rectifier_potassium_current_conductance',
     'rapid_time_dependent_potassium_current_conductance', 'rapid_time_dependent_potassium_current_Xr1_gate',
     'rapid_time_dependent_potassium_current_Xr2_gate', 'slow_time_dependent_potassium_current_conductance',
     'slow_time_dependent_potassium_current_Xs_gate', 'sodium_background_current_conductance',
     'membrane_L_type_calcium_current_conductance', 'membrane_L_type_calcium_current_d_gate', 'membrane_L_type_calcium_current_f_gate', 'L_type_Ca_current_f2_gate',
     'membrane_L_type_calcium_current_fCass_gate', 'calcium_background_current_conductance', 'membrane_transient_outward_current_conductance',
     'membrane_transient_outward_current_s_gate', 'membrane_transient_outward_current_r_gate', 'membrane_sodium_potassium_pump_current_permeability',
     'membrane_sodium_calcium_exchanger_current_maximum', 'calcium_pump_current_conductance', 'membrane_potassium_pump_current_conductance',
     'calcium_dynamics_release_current_maximum', 'calcium_dynamics_leak_current_maximum', 'calcium_leak_current_conductance',
     'calcium_dynamics_uptake_current_maximum',
     'cytosolic_calcium_concentration',
     'state_variable',
     # What follows are some more metadata names agreed by Gary Mirams and Penny Noble
     'cytosolic_potassium_concentration', 'cytosolic_sodium_concentration',
     'cytoplasmic_potassium_concentration', 'cytoplasmic_sodium_concentration',
     'membrane_fast_sodium_current_shift_inactivation', 'membrane_fast_sodium_current_reduced_inactivation', 'membrane_fast_sodium_current_h_gate_tau', 'membrane_fast_sodium_current_j_gate_tau',
     'membrane_L_type_calcium_current', 'membrane_L_type_calcium_current_fCa_gate', 'membrane_L_type_calcium_current_fCa2_gate', 'membrane_L_type_calcium_current_f2_gate', 'membrane_L_type_calcium_current_f2ds_gate', 'membrane_L_type_calcium_current_d2_gate', 'membrane_L_type_calcium_current_f_gate_tau', 'membrane_L_type_calcium_current_f2_gate_tau', 'membrane_L_type_calcium_current_fCa_gate_tau', 'membrane_L_type_calcium_current_fCa2_gate_tau',
     'membrane_rapid_delayed_rectifier_potassium_current', 'membrane_rapid_delayed_rectifier_potassium_current_conductance1', 'membrane_rapid_delayed_rectifier_potassium_current_conductance2',
     'membrane_sodium_calcium_exchanger_current', 'membrane_sodium_calcium_exchanger_current_conductance',
     'membrane_delayed_rectifier_potassium_current', 'membrane_delayed_rectifier_potassium_current_conductance',
     'membrane_transient_outward_current',
     'membrane_inward_rectifier_potassium_current',
     'membrane_slow_delayed_rectifier_potassium_current',
     'calcium_concentration_intracellular',
     'diadicspace_calcium_concentration',
     'SR_release_current',
     'concentration_clamp_onoff',
     'extracellular_potassium_concentration', 'extracellular_calcium_concentration', 'extracellular_sodium_concentration', 'bath_potassium_concentration',
     'SR_leak_current_max', 'SR_release_current_max', 'SR_release_kmcacyt', 'SR_release_kmcads', 'SR_uptake_current_max', 'JSR_calcium_concentration'
])

# Parameters for the stimulus current
STIMULUS_NAMES = frozenset('membrane_stimulus_current_'+ v for v in ['duration', 'amplitude', 'period', 'offset', 'end'])

################################################################################
# The public interface to this module
################################################################################

def remove_model(cellml_model):
    """The given model is being deleted / no longer needed."""
    return _wrapper.remove_model(cellml_model)

def update_serialized_rdf(cellml_model):
    """Ensure the RDF serialized into the given CellML model is up-to-date.
    
    If we have done any metadata processing on the given model, will serialize
    our RDF store into the rdf:RDF element child of the model.
    """
    return _wrapper.update_serialized_rdf(cellml_model)

def create_rdf_node(node_content=None, fragment_id=None):
    """Create an RDF node.
    
    node_content, if given, must either be a tuple (qname, namespace_uri),
    or a string, in which case it is interpreted as a literal RDF node.
    
    Alternatively, fragment_id may be given to refer to a cmeta:id within the
    current model.
    
    If neither are given, a blank node is created.
    """
    return _wrapper.create_rdf_node(node_content, fragment_id)

def create_unique_id(cellml_model, base_id):
    """Create a fragment identifier that hasn't already been used.
    
    If base_id hasn't been used, it will be returned.  Otherwise, underscores will
    be added until a unique id is obtained.
    """
    return _wrapper.create_unique_id(cellml_model, base_id)

def replace_statement(cellml_model, source, property, target):
    """Add a statement to the model, avoiding duplicates.
    
    Any existing statements with the same source and property will first be removed.
    """
    _debug("replace_statement(", source, ",", property, ",", target, ")")
    return _wrapper.replace_statement(cellml_model, source, property, target)

def add_statement(cellml_model, source, property, target):
    """Add a statement to the model."""
    _debug("add_statement(", source, ",", property, ",", target, ")")
    return _wrapper.add_statement(cellml_model, source, property, target)

def remove_statements(cellml_model, source, property, target):
    """Remove all statements matching (source,property,target).
    
    Any of these may be None to match anything.
    """
    _debug("remove_statements(", source, ",", property, ",", target, ")")
    return _wrapper.remove_statements(cellml_model, source, property, target)

def get_target(cellml_model, source, property):
    """Get the target of property from source.
    
    Returns None if no such target exists.  Throws if there is more than one match.
    
    If the target is a literal node, returns its string value.  Otherwise returns an RDF node.
    """
    return _wrapper.get_target(cellml_model, source, property)

def get_targets(cellml_model, source, property):
    """Get a list of all targets of property from source.
    
    If no such targets exist, returns an empty list.
    If property is None, targets of any property will be returned.
    Alternatively if source is None, targets of the given property from any source will be found.
    
    For each target, if it is a literal node then its string value is given.
    Otherwise the list will contain an RDF node.
    """
    return _wrapper.get_targets(cellml_model, source, property)

def find_variables(cellml_model, property, value=None):
    """Find variables in the cellml_model with the given property, and optionally value.
    
    property (and value if given) should be a suitable input for create_rdf_node.
    
    Will return a list of cellml_variable instances.
    """
    _debug("find_variables(", property, ",", value, ")")
    return _wrapper.find_variables(cellml_model, property, value)

def namespace_member(node, nsuri, not_uri_ok=False, wrong_ns_ok=False):
    """Given a URI reference RDF node and namespace URI, return the local part.
    
    Will raise an exception if node is not a URI reference unless not_uri_ok is True.
    Will raise an exception if the node doesn't live in the given namespace, unless
    wrong_ns_ok is True.  In both cases, if the error is suppressed the empty string
    will be returned instead.
    """
    local_part = _wrapper.namespace_member(node, nsuri, not_uri_ok, wrong_ns_ok)
    _debug("namespace_member(", node, ",", nsuri, ") = ", local_part)
    return local_part

################################################################################
# Implementation
################################################################################


def _debug(*args):
    pycml.DEBUG('cellml-metadata', *args)

_must_provide = "\n    Must be implemented by subclasses.\n"

class RdfWrapper(object):
    """Base class for wrappers around particular RDF libraries."""
    
    def __init__(self):
        """Create the RDF handler."""
        # Map from cellml_model instances to RDF stores
        self._models = {}
        
    def _create_new_store(self, cellml_model):
        """Create a new RDF store for the given CellML model.
        The new store will be available as self._models[cellml_model].
        
        Must be implemented by subclasses.
        """
        raise NotImplementedError
    
    def _add_rdf_element(self, cellml_model, rdf_text):
        """Add statements to the model's graph from the given serialized RDF.
        
        Must be implemented by subclasses.
        """
        raise NotImplementedError
    
    def _serialize(self, cellml_model):
        """Serialize the RDF model for this CellML model to XML.
        
        Must be implemented by subclasses.
        """
        raise NotImplementedError

    def get_rdf_from_model(self, cellml_model):
        """Get the RDF graph of the given CellML model.
    
        If this model is already in our map, return the existing RDF store.
        Otherwise, extract metadata from all RDF elements in the cellml_model,
        create a new RDF graph from these, and delete the original elements.
        """
        if not cellml_model in self._models:
            rdf_blocks = cellml_model.xml_xpath(u'//rdf:RDF')
            self._create_new_store(cellml_model)
            for rdf_block in rdf_blocks:
                rdf_text = rdf_block.xml()
                self._add_rdf_element(cellml_model, rdf_text)
                rdf_block.xml_parent.xml_remove_child(rdf_block)
        return self._models[cellml_model]
    
    def remove_model(self, cellml_model):
        if cellml_model in self._models:
            del self._models[cellml_model]
            _debug('Clearing RDF state for model', cellml_model.name)
    remove_model.__doc__ = globals()['remove_model'].__doc__

    def update_serialized_rdf(self, cellml_model):
        if cellml_model in self._models:
            # Paranoia: ensure it doesn't already contain serialized RDF
            rdf_blocks = cellml_model.xml_xpath(u'//rdf:RDF')
            if rdf_blocks:
                pycml.LOG('cellml-metadata', logging.WARNING, 'Removing existing RDF in model.')
                for rdf_block in rdf_blocks:
                    rdf_block.xml_parent.xml_remove_child(rdf_block)
            # Serialize the RDF model into cellml_model.RDF
            rdf_text = self._serialize(cellml_model)
            rdf_doc = pycml.amara.parse(rdf_text)
            cellml_model.xml_append(rdf_doc.RDF)
            # Remove the RDF model
            self.remove_model(cellml_model)
    update_serialized_rdf.__doc__ = globals()['update_serialized_rdf'].__doc__

    def create_rdf_node(self, node_content=None, fragment_id=None):
        raise NotImplementedError
    create_rdf_node.__doc__ = globals()['create_rdf_node'].__doc__ + _must_provide
    
    def create_unique_id(self, cellml_model, base_id):
        while True:
            node = self.create_rdf_node(fragment_id=base_id)
            if not self.get_targets(cellml_model, node, None):
                break
            base_id += u'_'
        return base_id
    create_unique_id.__doc__ = globals()['create_unique_id'].__doc__

    def replace_statement(self, cellml_model, source, property, target):
        raise NotImplementedError
    replace_statement.__doc__ = globals()['replace_statement'].__doc__ + _must_provide

    def remove_statements(self, cellml_model, source, property, target):
        raise NotImplementedError
    remove_statements.__doc__ = globals()['remove_statements'].__doc__ + _must_provide

    def get_target(self, cellml_model, source, property):
        raise NotImplementedError
    get_target.__doc__ = globals()['get_target'].__doc__ + _must_provide

    def get_targets(self, cellml_model, source, property):
        raise NotImplementedError
    get_targets.__doc__ = globals()['get_targets'].__doc__ + _must_provide

    def find_variables(self, cellml_model, property, value=None):
        raise NotImplementedError
    find_variables.__doc__ = globals()['find_variables'].__doc__ + _must_provide

    def namespace_member(self, node, nsuri, not_uri_ok=False, wrong_ns_ok=False):
        raise NotImplementedError
    namespace_member.__doc__ = globals()['namespace_member'].__doc__ + _must_provide


####################################################################################
# Wrapper using the Redland RDF library
####################################################################################

class RedlandWrapper(RdfWrapper):
    """Implements CellML metadata functionality using the Redland RDF library."""

    # Base URI to use for models.  Unfortunately the RDF library won't let
    # us use an empty URI, so we use a dummy URI then strip it out when
    # serializing the RDF.
    _base_uri = 'urn:chaste-pycml:dummy-rdf-base-uri'
    
    def __init__(self):
        """Create the wrapper."""
        assert RDF is not None, "Redland RDF library is not available."
        super(RedlandWrapper, self).__init__()
        self._parser = None

    def _create_new_store(self, cellml_model):
        """Create a new RDF store for the given CellML model.
        The new store will be available as self._models[cellml_model].
        """
        self._models[cellml_model] = RDF.Model()
        if not self._parser:
            self._parser = RDF.Parser()

    def remove_model(self, cellml_model):
        super(RedlandWrapper, self).remove_model(cellml_model)
        if self._parser:
            del self._parser
            self._parser = None
    remove_model.__doc__ = globals()['remove_model'].__doc__
    
    def _add_rdf_element(self, cellml_model, rdf_text):
        """Add statements to the model's graph from the given serialized RDF."""
        self._parser.parse_string_into_model(self._models[cellml_model], rdf_text, self._base_uri)
    
    def _serialize(self, cellml_model):
        """Serialize the RDF model for this CellML model to XML."""
        return self._models[cellml_model].to_string().replace(self._base_uri, '')
    
    def create_rdf_node(self, node_content=None, fragment_id=None):
        if fragment_id:
            node = RDF.Node(uri_string=str(self._base_uri+'#'+fragment_id))
        elif node_content:
            if type(node_content) == types.TupleType:
                qname, nsuri = node_content
                if nsuri[-1] not in ['#', '/']:
                    nsuri = nsuri + '#'
                prefix, local_name = pycml.SplitQName(qname)
                node = RDF.Node(uri_string=str(nsuri+local_name))
            elif type(node_content) in types.StringTypes:
                node = RDF.Node(str(node_content))
            else:
                raise ValueError("Don't know how to make a node from " + str(node_content)
                                 + " of type " + type(node_content))
        else:
            node = RDF.Node()
        return node
    create_rdf_node.__doc__ = globals()['create_rdf_node'].__doc__

    def add_statement(self, cellml_model, source, property, target):
        rdf_model = self.get_rdf_from_model(cellml_model)
        # Add the new statement
        statement = RDF.Statement(subject=source, predicate=property, object=target)
        rdf_model.append(statement)
    add_statement.__doc__ = globals()['add_statement'].__doc__

    def replace_statement(self, cellml_model, source, property, target):
        rdf_model = self.get_rdf_from_model(cellml_model)
        # Check for existing statements
        query = RDF.Statement(subject=source, predicate=property, object=None)
        for statement in rdf_model.find_statements(query):
            del rdf_model[statement]
        # Add the new statement
        statement = RDF.Statement(subject=source, predicate=property, object=target)
        rdf_model.append(statement)
    replace_statement.__doc__ = globals()['replace_statement'].__doc__

    def remove_statements(self, cellml_model, source, property, target):
        rdf_model = self.get_rdf_from_model(cellml_model)
        query = RDF.Statement(subject=source, predicate=property, object=target)
        for statement in rdf_model.find_statements(query):
            del rdf_model[statement]
    remove_statements.__doc__ = globals()['remove_statements'].__doc__

    def get_target(self, cellml_model, source, property):
        rdf_model = self.get_rdf_from_model(cellml_model)
        targets = list(rdf_model.targets(source, property))
        if len(targets) > 1:
            raise ValueError("Too many targets for source " + str(source) + " and property " + str(property))
        elif len(targets) == 1:
            target = targets[0]
        else:
            target = None
        if target and target.is_literal():
            target = target.literal_value['string']
        _debug("get_target(", source, ",", property, ") -> ", "'" + str(target) + "'")
        return target
    get_target.__doc__ = globals()['get_target'].__doc__

    def get_targets(self, cellml_model, source, property):
        rdf_model = self.get_rdf_from_model(cellml_model)
        query = RDF.Statement(source, property, None)
        targets = [stmt.object for stmt in rdf_model.find_statements(query)]
        for i, target in enumerate(targets):
            if target.is_literal():
                targets[i] = target.literal_value['string']
        return targets
    get_targets.__doc__ = globals()['get_targets'].__doc__

    def find_variables(self, cellml_model, property, value=None):
        rdf_model = self.get_rdf_from_model(cellml_model)
        property = self.create_rdf_node(property)
        if value:
            value = self.create_rdf_node(value)
        query = RDF.Statement(None, property, value)
        results = rdf_model.find_statements(query)
        vars = []
        for result in results:
            assert result.subject.is_resource(), "Non-resource annotated."
            uri = str(result.subject.uri)
            if uri.startswith(self._base_uri):
                # Strip the base URI part
                uri = uri[len(self._base_uri):]
            assert uri[0] == '#', "Annotation found on non-local URI"
            var_id = uri[1:] # Strip '#'
            var_objs = cellml_model.xml_xpath(u'*/cml:variable[@cmeta:id="%s"]' % var_id)
            assert len(var_objs) == 1, "Didn't find a unique variable with ID " + var_id
            vars.append(var_objs[0])
        return vars
    find_variables.__doc__ = globals()['find_variables'].__doc__
    
    def namespace_member(self, node, nsuri, not_uri_ok=False, wrong_ns_ok=False):
        if not node.is_resource():
            if not_uri_ok:
                return ""
            else:
                raise ValueError("Cannot extract namespace member for a non-URI RDF node.")
        uri = str(node.uri)
        if uri.startswith(nsuri):
            return uri[len(nsuri):]
        elif wrong_ns_ok:
            return ""
        else:
            raise ValueError("Node is not in correct namespace.")
    namespace_member.__doc__ = globals()['namespace_member'].__doc__

####################################################################################
# Wrapper using the RDFLib library
####################################################################################

class RdflibWrapper(RdfWrapper):
    """Implements CellML metadata functionality using the RDFLib library."""
    def __init__(self):
        """Create the wrapper."""
        assert rdflib is not None, "RDFLib library is not available."
        super(RdflibWrapper, self).__init__()
        # Cope with differences in API?
        self.Graph = rdflib.ConjunctiveGraph
        self.URIRef = rdflib.URIRef
        self.Literal = rdflib.Literal
        self.BNode = rdflib.BNode
        self.Namespace = rdflib.Namespace

    def _create_new_store(self, cellml_model):
        """Create a new RDF store for the given CellML model.
        The new store will be available as self._models[cellml_model].
        """
        self._models[cellml_model] = self.Graph()
    
    def _add_rdf_element(self, cellml_model, rdf_text):
        """Add statements to the model's graph from the given serialized RDF."""
        g = self.Graph()
        g.parse(StringIO(rdf_text))
        rdf_model = self._models[cellml_model]
        for stmt in g:
            rdf_model.add(stmt)
    
    def _serialize(self, cellml_model):
        """Serialize the RDF model for this CellML model to XML."""
        return self._models[cellml_model].serialize()
    
    def create_rdf_node(self, node_content=None, fragment_id=None):
        if fragment_id:
            node = self.URIRef(str('#'+fragment_id))
        elif node_content:
            if type(node_content) == types.TupleType:
                qname, nsuri = node_content
                if nsuri[-1] not in ['#', '/']:
                    nsuri = nsuri + '#'
                ns = self.Namespace(nsuri)
                prefix, local_name = pycml.SplitQName(qname)
                node = ns[local_name]
            elif type(node_content) in types.StringTypes:
                node = self.Literal(node_content)
            else:
                raise ValueError("Don't know how to make a node from " + str(node_content)
                                 + " of type " + type(node_content))
        else:
            node = self.BNode()
        return node
    create_rdf_node.__doc__ = globals()['create_rdf_node'].__doc__

    def add_statement(self, cellml_model, source, property, target):
        rdf_model = self.get_rdf_from_model(cellml_model)
        rdf_model.add((source, property, target))
    add_statement.__doc__ = globals()['add_statement'].__doc__

    def replace_statement(self, cellml_model, source, property, target):
        rdf_model = self.get_rdf_from_model(cellml_model)
        rdf_model.set((source, property, target))
    replace_statement.__doc__ = globals()['replace_statement'].__doc__

    def remove_statements(self, cellml_model, source, property, target):
        rdf_model = self.get_rdf_from_model(cellml_model)
        rdf_model.remove((source, property, target))
    remove_statements.__doc__ = globals()['remove_statements'].__doc__

    def get_target(self, cellml_model, source, property):
        rdf_model = self.get_rdf_from_model(cellml_model)
        try:
            target = rdf_model.value(subject=source, predicate=property, any=False)
        except rdflib.exceptions.UniquenessError:
            raise ValueError("Too many targets for source " + str(source) + " and property " + str(property))
        if isinstance(target, rdflib.Literal):
            target = str(target)
        _debug("get_target(", source, ",", property, ") -> ", "'" + str(target) + "'")
        return target
    get_target.__doc__ = globals()['get_target'].__doc__

    def get_targets(self, cellml_model, source, property):
        rdf_model = self.get_rdf_from_model(cellml_model)
        targets = list(rdf_model.objects(subject=source, predicate=property))
        for i, target in enumerate(targets):
            if isinstance(target, rdflib.Literal):
                targets[i] = str(target)
        return targets
    get_targets.__doc__ = globals()['get_targets'].__doc__

    def find_variables(self, cellml_model, property, value=None):
        rdf_model = self.get_rdf_from_model(cellml_model)
        property = self.create_rdf_node(property)
        if value:
            value = self.create_rdf_node(value)
        vars = []
        for result in rdf_model.subjects(property, value):
            assert isinstance(result, self.URIRef), "Non-resource annotated."
            uri = str(result)
            assert uri[0] == '#', "Annotation found on non-local URI"
            var_id = uri[1:] # Strip '#'
            var_objs = cellml_model.xml_xpath(u'*/cml:variable[@cmeta:id="%s"]' % var_id)
            assert len(var_objs) == 1, "Didn't find a unique variable with ID " + var_id
            vars.append(var_objs[0])
        return vars
    find_variables.__doc__ = globals()['find_variables'].__doc__

    def namespace_member(self, node, nsuri, not_uri_ok=False, wrong_ns_ok=False):
        if not isinstance(node, self.URIRef):
            if not_uri_ok:
                return ""
            else:
                raise ValueError("Cannot extract namespace member for a non-URI RDF node.")
        if node.startswith(nsuri):
            return node[len(nsuri):]
        elif wrong_ns_ok:
            return ""
        else:
            raise ValueError("Node is not in correct namespace.")
    namespace_member.__doc__ = globals()['namespace_member'].__doc__

####################################################################################
# Finally, instantiate a suitable wrapper instance
####################################################################################

if rdflib:
    _wrapper = RdflibWrapper()
else:
    _wrapper = RedlandWrapper()

