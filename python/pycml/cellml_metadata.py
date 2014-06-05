
"""Copyright (c) 2005-2014, University of Oxford.
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
This module abstracts the interface to RDF metadata about CellML models.

Public methods on a single instance of the RdfProcessor class below are also
exposed as module-level functions, and these should typically be called by users.
"""

import logging
import types
from cStringIO import StringIO

# We now only support rdflib for RDF processing
import rdflib

import pycml

# Allowed metadata names, more to come
# TODO #2547: Use a proper ontology!
METADATA_NAMES = frozenset(
    ['state_variable', 'time', 'temperature', 
     # =====================================================
     # Cardiac-Specific Labels
     # =====================================================
     # VOLTAGES / POTENTIALS
     'membrane_voltage', 
     # These are intended to refer to outer cell membrane, new tags should be 
     # introduced for SR membrane, mitochondrial membranes etc.
     'sodium_reversal_potential', 'potassium_reversal_potential', 
     'calcium_reversal_potential', 'chloride_reversal_potential',
     # Membrane properties
     'membrane_capacitance', 'membrane_E_R',
     # =====================================================
     # Stimulus Current
     # =====================================================
     'membrane_stimulus_current', 
        'membrane_stimulus_current_duration',
        'membrane_stimulus_current_amplitude',
        'membrane_stimulus_current_period',
        'membrane_stimulus_current_offset',
        'membrane_stimulus_current_end',
     # =====================================================
     # IONIC CONCENTRATIONS
     # =====================================================
     # basic 'intracellular' and 'extracellular'
     'extracellular_potassium_concentration', 'extracellular_calcium_concentration', 
     'extracellular_sodium_concentration', 'extracellular_chloride_concentration',
     'cytosolic_calcium_concentration','cytosolic_potassium_concentration', 
     'cytosolic_sodium_concentration','cytosolic_chloride_concentration',
     # in Calcium subsystem SR = sarcoplasmic reticulum
     'SR_calcium_concentration', # Some models have just the SR,
       'JSR_calcium_concentration', # Other models divide it into the Junctional SR (near RyRs),
       'NSR_calcium_concentration', # and the Network SR (rest of SR),
     'diadicspace_calcium_concentration', # Some models also have a separate diadic sub-space (cytosol between JSR and t-tubules).
     # Others
     'bath_potassium_concentration',
     # =====================================================
     # CURRENTS
     # =====================================================
     # HISTORIC metadata - only for early models lacking components.
     'membrane_potassium_current', 
        'membrane_potassium_current_conductance', 
        'potassium_channel_n_gate', 
    'membrane_delayed_rectifier_potassium_current', 
        'membrane_delayed_rectifier_potassium_current_conductance',
     'rapid_time_dependent_potassium_current_conductance', 
        'rapid_time_dependent_potassium_current_Xr1_gate',
        'rapid_time_dependent_potassium_current_Xr2_gate', 
     'slow_time_dependent_potassium_current_conductance',
        'slow_time_dependent_potassium_current_Xs_gate', 
     'membrane_slow_inward_current',
        'membrane_slow_inward_current_conductance',
     'leakage_current',
     # MODERN metadata - labels all new models should be able to use.
     # ========================================================================
     # SODIUM CURRENTS
     # ========================================================================
     # I Na (fast)
     'membrane_fast_sodium_current', 
        'membrane_fast_sodium_current_conductance', 
        'membrane_fast_sodium_current_m_gate', 
        'membrane_fast_sodium_current_h_gate',
            'membrane_fast_sodium_current_h_gate_tau',
        'membrane_fast_sodium_current_j_gate', 
            'membrane_fast_sodium_current_j_gate_tau',
        'membrane_fast_sodium_current_shift_inactivation', 'membrane_fast_sodium_current_reduced_inactivation', 
     # I_Na_L (late or persistent)
     'membrane_persistent_sodium_current', 
        'membrane_persistent_sodium_current_conductance', 
     # I Na,b (background)
     'membrane_background_sodium_current',
        'membrane_background_sodium_current_conductance',
     # ========================================================================
     # Potassium currents
     # ========================================================================
     # I Kr
     'membrane_rapid_delayed_rectifier_potassium_current', 
        'membrane_rapid_delayed_rectifier_potassium_current_conductance',
        'membrane_rapid_delayed_rectifier_potassium_current_conductance1', 
        'membrane_rapid_delayed_rectifier_potassium_current_conductance2',
     # I Ks
     'membrane_slow_delayed_rectifier_potassium_current',   
        'membrane_slow_delayed_rectifier_potassium_current_conductance',
        'membrane_slow_delayed_rectifier_potassium_current_xs1_gate_tau', # (really scaling factor for tau)
        'membrane_slow_delayed_rectifier_potassium_current_xs2_gate_tau', # (really scaling factor for tau)
     # I_Kur
     'membrane_ultrarapid_delayed_rectifier_potassium_current',
        'membrane_ultrarapid_delayed_rectifier_potassium_current_conductance',
     # I_Kss or Iss   
     'membrane_non_inactivating_steady_state_potassium_current',
        'membrane_non_inactivating_steady_state_potassium_current_conductance',   
     # I K1
     'membrane_inward_rectifier_potassium_current',     
        'membrane_inward_rectifier_potassium_current_conductance',
     # I to, sometimes fast and slow components, sometimes not.
     'membrane_transient_outward_current',
        'membrane_transient_outward_current_conductance',
        'membrane_fast_transient_outward_current',
            'membrane_transient_outward_current_r_gate',
            'membrane_fast_transient_outward_current_conductance',
        'membrane_slow_transient_outward_current',
            'membrane_transient_outward_current_s_gate',    
            'membrane_slow_transient_outward_current_conductance',       
        'membrane_transient_outward_current_time_independent_rectification_gate_constant',
     # I katp
     'membrane_atp_dependent_potassium_current',
        'membrane_atp_dependent_potassium_current_conductance',  
     # I K,b (background current / leak)
     'membrane_background_potassium_current',
        'membrane_background_potassium_current_conductance',
     # ========================================================================
     # Mixed Currents
     # ========================================================================   
     # I f (funny current)
     # Generally, but not always, this is formulated as I_f = I_f_Na + I_f_K, with separate conductances.
     # if it isn't try and figure out which ionic species is being modelled for tagging, and give two tags if necessary.
     'membrane_hyperpolarisation_activated_funny_current',
        'membrane_hyperpolarisation_activated_funny_current_single_gate',
        'membrane_hyperpolarisation_activated_funny_current_potassium_component',
           'membrane_hyperpolarisation_activated_funny_current_potassium_component_conductance',
        'membrane_hyperpolarisation_activated_funny_current_sodium_component', 
           'membrane_hyperpolarisation_activated_funny_current_sodium_component_conductance', 
     # ICaL conductance of non-calcium ions:
     # Things here are getting a bit confusing, have to be careful as tags may have different
     # effects in different models. i.e. does main 'membrane_L_type_calcium_current_conductance'
     # scale all these as well, or are they treated as completely separate ion currents (as per O'Hara).
     'membrane_L_type_calcium_channel_sodium_current',
     'membrane_L_type_calcium_channel_sodium_current_conductance',
     'membrane_L_type_calcium_channel_potassium_current',
     'membrane_L_type_calcium_channel_potassium_current_conductance',
     # ========================================================================
     # CALCIUM CURRENTS
     # ========================================================================
     # I CaL
     'membrane_L_type_calcium_current', 
        'membrane_L_type_calcium_current_conductance', 
        'membrane_L_type_calcium_current_d_gate', 
        'membrane_L_type_calcium_current_f_gate', 
        'membrane_L_type_calcium_current_fCass_gate',
        'membrane_L_type_calcium_current_fCa_gate', 
        'membrane_L_type_calcium_current_fCa2_gate', 
        'membrane_L_type_calcium_current_f2_gate', 
        'membrane_L_type_calcium_current_f2ds_gate', 
        'membrane_L_type_calcium_current_d2_gate', 
        'membrane_L_type_calcium_current_f_gate_tau', 
        'membrane_L_type_calcium_current_f2_gate_tau', 
        'membrane_L_type_calcium_current_fCa_gate_tau', 
        'membrane_L_type_calcium_current_fCa2_gate_tau',
        'membrane_L_type_calcium_current_d_gate_power_tau',
     # I Ca,b (background)
     'membrane_background_calcium_current',
        'membrane_background_calcium_current_conductance', 
     # ========================================================================
     # Calcium subsystem parameters - needs tidying up.
     # ========================================================================
     'SR_release_current', # a.k.a. Jrel or RyR channel current 
     'SR_uptake_current',# a.k.a. Jup or SERCA current
     'SR_leak_current', 
     'SR_leak_current_max', 'SR_release_current_max', 'SR_uptake_current_max', 
     'SR_release_kmcacyt', 'SR_release_kmcads', 
     'calcium_dynamics_release_current_maximum', 'calcium_dynamics_leak_current_maximum', 
     'calcium_leak_current_conductance', 'calcium_dynamics_uptake_current_maximum',
     # ========================================================================     
     # Pumps and Exchangers
     # ========================================================================
      # I NCX
      'membrane_sodium_calcium_exchanger_current', 
        'membrane_sodium_calcium_exchanger_current_conductance', # a.k.a. permeability
      'SR_sodium_calcium_exchanger_current', 
        'SR_sodium_calcium_exchanger_current_conductance', # a.k.a. permeability
      # INaK
      'membrane_sodium_potassium_pump_current',
          'membrane_sodium_potassium_pump_current_permeability', # often INaK_max
      # Ip,Ca
      'membrane_calcium_pump_current',
         'membrane_calcium_pump_current_conductance', # a.k.a. permeability
      # Ip,K 
      'membrane_potassium_pump_current',
         'membrane_potassium_pump_current_conductance',      # a.k.a. permeability
      # Penny and Alan, protocol-specific stuff (to be replaced by Functional Curation in the end)
      'concentration_clamp_onoff',
])

# Parameters for the stimulus current
STIMULUS_NAMES = frozenset('membrane_stimulus_current_'+ v for v in ['duration', 'amplitude', 'period', 'offset', 'end'])


def _debug(*args):
    pycml.DEBUG('cellml-metadata', *args)


class RdfProcessor(object):
    """Implements CellML metadata functionality using the RDFLib library."""
    def __init__(self):
        """Create the wrapper."""
        # Map from cellml_model instances to RDF stores
        self._models = {}
        # Cope with differences in API between library versions
        rdflib_major_version = int(rdflib.__version__[0])
        if rdflib_major_version >= 3:
            self.Graph = rdflib.Graph
        else:
            self.Graph = rdflib.ConjunctiveGraph

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
        """The given model is being deleted / no longer needed."""
        if cellml_model in self._models:
            del self._models[cellml_model]
            _debug('Clearing RDF state for model', cellml_model.name)

    def update_serialized_rdf(self, cellml_model):
        """Ensure the RDF serialized into the given CellML model is up-to-date.
        
        If we have done any metadata processing on the given model, will serialize
        our RDF store into the rdf:RDF element child of the model.
        """
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

    def create_rdf_node(self, node_content=None, fragment_id=None):
        """Create an RDF node.
        
        node_content, if given, must either be a tuple (qname, namespace_uri),
        or a string, in which case it is interpreted as a literal RDF node.
        
        Alternatively, fragment_id may be given to refer to a cmeta:id within the
        current model.
        
        If neither are given, a blank node is created.
        """
        if fragment_id:
            node = rdflib.URIRef(str('#'+fragment_id))
        elif node_content:
            if type(node_content) == types.TupleType:
                qname, nsuri = node_content
                if nsuri[-1] not in ['#', '/']:
                    nsuri = nsuri + '#'
                ns = self.Namespace(nsuri)
                prefix, local_name = pycml.SplitQName(qname)
                node = ns[local_name]
            elif type(node_content) in types.StringTypes:
                node = rdflib.Literal(node_content)
            else:
                raise ValueError("Don't know how to make a node from " + str(node_content)
                                 + " of type " + type(node_content))
        else:
            node = rdflib.BNode()
        return node

    def create_unique_id(self, cellml_model, base_id):
        """Create a fragment identifier that hasn't already been used.
        
        If base_id hasn't been used, it will be returned.  Otherwise, underscores will
        be added until a unique id is obtained.
        """
        while True:
            node = self.create_rdf_node(fragment_id=base_id)
            if not self.get_targets(cellml_model, node, None):
                break
            base_id += u'_'
        return base_id

    def add_statement(self, cellml_model, source, property, target):
        """Add a statement to the model."""
        _debug("add_statement(", source, ",", property, ",", target, ")")
        rdf_model = self.get_rdf_from_model(cellml_model)
        rdf_model.add((source, property, target))

    def replace_statement(self, cellml_model, source, property, target):
        """Add a statement to the model, avoiding duplicates.
        
        Any existing statements with the same source and property will first be removed.
        """
        _debug("replace_statement(", source, ",", property, ",", target, ")")
        rdf_model = self.get_rdf_from_model(cellml_model)
        rdf_model.set((source, property, target))

    def remove_statements(self, cellml_model, source, property, target):
        """Remove all statements matching (source,property,target).
        
        Any of these may be None to match anything.
        """
        _debug("remove_statements(", source, ",", property, ",", target, ")")
        rdf_model = self.get_rdf_from_model(cellml_model)
        rdf_model.remove((source, property, target))

    def get_target(self, cellml_model, source, property):
        """Get the target of property from source.
        
        Returns None if no such target exists.  Throws if there is more than one match.
        
        If the target is a literal node, returns its string value.  Otherwise returns an RDF node.
        """
        rdf_model = self.get_rdf_from_model(cellml_model)
        try:
            target = rdf_model.value(subject=source, predicate=property, any=False)
        except rdflib.exceptions.UniquenessError:
            raise ValueError("Too many targets for source " + str(source) + " and property " + str(property))
        if isinstance(target, rdflib.Literal):
            target = str(target)
        _debug("get_target(", source, ",", property, ") -> ", "'" + str(target) + "'")
        return target

    def get_targets(self, cellml_model, source, property):
        """Get a list of all targets of property from source.
        
        If no such targets exist, returns an empty list.
        If property is None, targets of any property will be returned.
        Alternatively if source is None, targets of the given property from any source will be found.
        
        For each target, if it is a literal node then its string value is given.
        Otherwise the list will contain an RDF node.
        """
        rdf_model = self.get_rdf_from_model(cellml_model)
        targets = list(rdf_model.objects(subject=source, predicate=property))
        for i, target in enumerate(targets):
            if isinstance(target, rdflib.Literal):
                targets[i] = str(target)
        return targets

    def find_variables(self, cellml_model, property, value=None):
        """Find variables in the cellml_model with the given property, and optionally value.
        
        property (and value if given) should be a suitable input for create_rdf_node.
        
        Will return a list of cellml_variable instances.
        """
        _debug("find_variables(", property, ",", value, ")")
        rdf_model = self.get_rdf_from_model(cellml_model)
        property = self.create_rdf_node(property)
        if value:
            value = self.create_rdf_node(value)
        vars = []
        for result in rdf_model.subjects(property, value):
            assert isinstance(result, rdflib.URIRef), "Non-resource annotated."
            uri = str(result)
            assert uri[0] == '#', "Annotation found on non-local URI"
            var_id = uri[1:] # Strip '#'
            var_objs = cellml_model.xml_xpath(u'*/cml:variable[@cmeta:id="%s"]' % var_id)
            assert len(var_objs) == 1, "Didn't find a unique variable with ID " + var_id
            vars.append(var_objs[0])
        return vars
    
    def get_all_rdf(self, cellml_model):
        """Return an iterator over all RDF triples in the model."""
        rdf_model = self.get_rdf_from_model(cellml_model)
        for triple in rdf_model:
            yield triple

    def namespace_member(self, node, nsuri, not_uri_ok=False, wrong_ns_ok=False):
        """Given a URI reference RDF node and namespace URI, return the local part.
        
        Will raise an exception if node is not a URI reference unless not_uri_ok is True.
        Will raise an exception if the node doesn't live in the given namespace, unless
        wrong_ns_ok is True.  In both cases, if the error is suppressed the empty string
        will be returned instead.
        """
        local_part = ""
        if not isinstance(node, rdflib.URIRef):
            if not not_uri_ok:
                raise ValueError("Cannot extract namespace member for a non-URI RDF node.")
        if node.startswith(nsuri):
            local_part = node[len(nsuri):]
        elif not wrong_ns_ok:
            raise ValueError("Node is not in correct namespace.")
        _debug("namespace_member(", node, ",", nsuri, ") = ", local_part)
        return local_part

####################################################################################
# Finally, instantiate a single processor instance and expose its methods
####################################################################################

_instance = RdfProcessor()

for attr in dir(_instance):
    if attr[0] != '_':
        meth = getattr(_instance, attr)
        if isinstance(meth, types.MethodType):
            globals()[attr] = meth
