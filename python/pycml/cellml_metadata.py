
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
This module abstracts the interface to RDF metadata about CellML models.

The RdfProcessor class below pretends to be the module itself, so all its properties
are available at module-level, and these should typically be called by users.

It also provides sets METADATA_NAMES and STIMULUS_NAMES, which contain the local names
of terms in the ontology that can annotate variables, and the subset of those names
which define properties of the stimulus current (but not the current itself), respectively.
"""

import logging
import os
import sys
import types
from cStringIO import StringIO

# We now only support rdflib for RDF processing
import rdflib


def __init__(module):
    # Import pycml here, to avoid circular import surprises
    import pycml
    module.pycml = pycml


class RdfProcessor(object):
    """Implements CellML metadata functionality using the RDFLib library."""
    def __init__(self, name):
        """Create the wrapper."""
        # Magic for pretending to be a module
        self._module = sys.modules[name]
        sys.modules[name] = self
        self._initializing = True
        # Map from cellml_model instances to RDF stores
        self._models = {}
        # Oxford metadata will be loaded lazily
        self._metadata_names = self._stimulus_names = self._ontology = None
        # Cope with differences in API between library versions
        rdflib_major_version = int(rdflib.__version__[0])
        if rdflib_major_version >= 3:
            self.Graph = rdflib.Graph
            self.Node = rdflib.term.Node
        else:
            self.Graph = rdflib.ConjunctiveGraph
            self.Node = rdflib.Node.Node

    def __getattribute__(self, name):
        """Provide access to real module-level variables as though they're class properties."""
        # call module.__init__ after import introspection is done
        baseget = super(RdfProcessor, self).__getattribute__
        module = baseget('_module')
        if baseget('_initializing') and not name[:2] == '__' == name[-2:]:
            setattr(self, '_initializing', False)
            __init__(module)
        try:
            return baseget(name)
        except AttributeError:
            return getattr(module, name)

    def _debug(*args):
        pycml.DEBUG('cellml-metadata', *args)

    def _load_ontology(self):
        """Load the Oxford metadata ontology the first time it's needed."""
        pycml_path = os.path.dirname(os.path.realpath(__file__))
        oxmeta_ttl = os.path.join(pycml_path, 'oxford-metadata.ttl')
        oxmeta_rdf = os.path.join(pycml_path, 'oxford-metadata.rdf')

        g = self._ontology = self.Graph()
        # We allow a difference in modification time of 10s, so we don't get confused when checking out!
        if os.stat(oxmeta_ttl).st_mtime > os.stat(oxmeta_rdf).st_mtime + 10.0:
            # Try to regenerate RDF/XML version of ontology
            try:
                g.parse(oxmeta_ttl, format='turtle')
            except Exception, e:
                print >> sys.stderr, 'Unable to convert metadata from Turtle format to RDF/XML.'
                print >> sys.stderr, 'Probably you need to upgrade rdflib to version 4.\nDetails of error:'
                raise
            g.serialize(oxmeta_rdf, format='xml')
        else:
            # Just parse the RDF/XML version
            g.parse(oxmeta_rdf, format='xml')
        
        annotation_terms = list(g.subjects(rdflib.RDF.type, rdflib.URIRef(pycml.NSS['oxmeta']+u'Annotation')))
        self._metadata_names = frozenset(map(lambda node: self.namespace_member(node, pycml.NSS['oxmeta']), annotation_terms))
        
        # Parameters for the stimulus current
        self._stimulus_names = frozenset(filter(lambda name: name.startswith('membrane_stimulus_current_'), self._metadata_names))

    @property
    def METADATA_NAMES(self):
        """Fake a module-level constant as a property for lazy loading."""
        if self._metadata_names is None:
            self._load_ontology()
        return self._metadata_names

    @property
    def STIMULUS_NAMES(self):
        """Fake a module-level constant as a property for lazy loading."""
        if self._stimulus_names is None:
            self._load_ontology()
        return self._stimulus_names

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
            self._debug('Clearing RDF state for model', cellml_model.name)

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
        
        node_content, if given, must either be a self.Node instance, a tuple (qname, namespace_uri),
        or a string, in which case it is interpreted as a literal RDF node.
        
        Alternatively, fragment_id may be given to refer to a cmeta:id within the current model.
        
        If neither are given, a blank node is created.
        """
        if fragment_id:
            node = rdflib.URIRef(str('#'+fragment_id))
        elif node_content:
            if type(node_content) == types.TupleType:
                qname, nsuri = node_content
                if nsuri[-1] not in ['#', '/']:
                    nsuri = nsuri + '#'
                ns = rdflib.Namespace(nsuri)
                prefix, local_name = pycml.SplitQName(qname)
                node = ns[local_name]
            elif type(node_content) in types.StringTypes:
                node = rdflib.Literal(node_content)
            elif isinstance(node_content, self.Node):
                node = node_content
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
        self._debug("add_statement(", source, ",", property, ",", target, ")")
        rdf_model = self.get_rdf_from_model(cellml_model)
        rdf_model.add((source, property, target))

    def replace_statement(self, cellml_model, source, property, target):
        """Add a statement to the model, avoiding duplicates.
        
        Any existing statements with the same source and property will first be removed.
        """
        self._debug("replace_statement(", source, ",", property, ",", target, ")")
        rdf_model = self.get_rdf_from_model(cellml_model)
        rdf_model.set((source, property, target))

    def remove_statements(self, cellml_model, source, property, target):
        """Remove all statements matching (source,property,target).
        
        Any of these may be None to match anything.
        """
        self._debug("remove_statements(", source, ",", property, ",", target, ")")
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
        self._debug("get_target(", source, ",", property, ") -> ", "'" + str(target) + "'")
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
        self._debug("find_variables(", property, ",", value, ")")
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
            assert len(var_objs) > 0, "Didn't find any variable with ID '" + var_id + "' when dereferencing annotation"
            assert len(var_objs) == 1, "Found " + str(len(var_objs)) + " variables with ID '" + var_id + "' when dereferencing annotation - IDs should be unique"
            vars.append(var_objs[0])
        return vars
    
    def transitive_subjects(self, term):
        """Transitively generate subjects connected to term by the rdf:type (a) property in the ontology."""
        if self._ontology is None:
            self._load_ontology()
        term = self.create_rdf_node(term)
        return self._ontology.transitive_subjects(rdflib.RDF.type, term)
    
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
        self._debug("namespace_member(", node, ",", nsuri, ") = ", local_part)
        return local_part

####################################################################################
# Instantiate a processor instance that pretends to be this module
####################################################################################

p = RdfProcessor(__name__)

if __name__ == '__main__':
    # Just load the ontology to trigger TTL -> RDF/XML conversion if needed
    p._load_ontology()
