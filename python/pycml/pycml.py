
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

# Initial work on a Python tool for processing CellML files.
# Eventual features:
#  - Featureful & draconian validation
#  - Apply various automatic optimisations, including PE & LUT
#  - Easy code generation
#  - Work around common problems in models, such as 0/0
#  - As in Memfem, convert from alpha&beta form to tau&inf

# This module contains code common to both validation and transformation.

# Ideas:
#  - CellML has the potential for introducing scripting functionality.
#    This may be a good way of incorporating LUT into the data model.
#    A special component could represent the table generation, and a
#    scripted function the lookup.
#    Alternatively we could add attributes in an extension namespace
#    to specify LUT parameters (var name, min, max, step).
#  - Add attributes in an extension namespace to represent binding
#    time annotations.
#  - We could do units conversions in a separate processing pass,
#    adding in extra mathematics to the CellML.

# Turn off DeprecationWarnings
import warnings
warnings.simplefilter('ignore', DeprecationWarning)

# Pythonic XML bindings
import amara
from amara import binderytools as bt

from Ft.Xml import SplitQName
from xml.dom import Node # For nodeType values

import copy
import itertools
import math
import operator
import sys
import types
from cStringIO import StringIO

from utilities import *
from _enum import Enum # Pythonic enums

import cellml_metadata # Handle RDF metadata for CellML

processors = None
def import_processors():
    """Lazy import."""
    global processors
    if processors is None:
        import processors


__version__ = "$Revision$"[11:-2]


######################################################################
#                               Logging                              #
######################################################################
import logging

# Default config for root logger
logging.basicConfig(level=logging.CRITICAL,
                    format="%(name)s: %(levelname)s: %(message)s",
                    stream=sys.stderr)
logging.getLogger().handlers[0].setLevel(logging.CRITICAL)

# Extra logging levels
# This level is a warning according to the spec, but an unrecoverable
# condition as far as translation is concerned.
logging.WARNING_TRANSLATE_ERROR = logging.WARNING + 5
logging.addLevelName(logging.WARNING_TRANSLATE_ERROR, 'WARNING')


# We specify some namespace prefixes; others are picked
# up automatically.  These are the standard namespaces we
# expect to see in CellML documents; a warning will be given
# if others are found.
NSS = {u'm'  : u'http://www.w3.org/1998/Math/MathML',
       u'cml': u'http://www.cellml.org/cellml/1.0#',
       # Our extensions; URIs will probably change?
       u'pe': u'https://chaste.comlab.ox.ac.uk/cellml/ns/partial-evaluation#',
       u'lut': u'https://chaste.comlab.ox.ac.uk/cellml/ns/lookup-tables',
       u'solver': u'https://chaste.comlab.ox.ac.uk/cellml/ns/solver-info',
       u'oxmeta': u'https://chaste.comlab.ox.ac.uk/cellml/ns/oxford-metadata#',
       u'pycml': u'https://chaste.comlab.ox.ac.uk/cellml/ns/pycml#',
       u'proto': u'https://chaste.cs.ox.ac.uk/nss/protocol/0.1#',
       # Metadata-related
       u'cmeta'  : u"http://www.cellml.org/metadata/1.0#",
       u'rdf'    : u"http://www.w3.org/1999/02/22-rdf-syntax-ns#",
       u'dc'     : u"http://purl.org/dc/elements/1.1/",
       u'dcterms': u"http://purl.org/dc/terms/",
       u'bqs'    : u"http://www.cellml.org/bqs/1.0#",
       u'vCard'  : u"http://www.w3.org/2001/vcard-rdf/3.0#",
       u'cg'     : u"http://www.cellml.org/metadata/graphs/1.0#",
       u'cs'     : u"http://www.cellml.org/metadata/simulation/1.0#",
       u'csub'   : u"http://www.cellml.org/metadata/custom_subset/1.0#",
       u'bqbiol' : u"http://biomodels.net/biology-qualifiers/",
       # Temporary documentation namespace
       u'doc' : u"http://cellml.org/tmp-documentation"
       }


# Variable classifications
VarTypes = Enum('Unknown', 'Free', 'State', 'MaybeConstant', 'Constant',
                'Computed', 'Mapped')

# Elements in the CellML subset of MathML
CELLML_SUBSET_ELTS = frozenset(
    ['math', 'cn', 'sep', 'ci', 'apply', 'piecewise', 'piece', 'otherwise',
     'eq', 'neq', 'gt', 'lt', 'geq', 'leq',
     'plus', 'minus', 'times', 'divide', 'power', 'root', 'abs',
     'exp', 'ln', 'log', 'floor', 'ceiling', 'factorial',
     'and', 'or', 'not', 'xor',
     'diff', 'degree', 'bvar', 'logbase',
     'sin', 'cos', 'tan', 'sec', 'csc', 'cot',
     'sinh', 'cosh', 'tanh', 'sech', 'csch', 'coth',
     'arcsin', 'arccos', 'arctan', 'arcsec', 'arccsc', 'arccot',
     'arcsinh', 'arccosh', 'arctanh', 'arcsech', 'arccsch', 'arccoth',
     'true', 'false', 'notanumber', 'pi', 'infinity', 'exponentiale',
     'semantics', 'annotation', 'annotation-xml'])

# Binding times for BTA
BINDING_TIMES = Enum('static', 'dynamic')


######################################################################
#                      Helpful utility functions                     #
######################################################################

def make_xml_binder():
    """
    Create a specialised binder, given some mappings from element names
    to python classes, and setting namespace prefixes.
    """
    binder = amara.bindery.binder(prefixes=NSS)
    binder.set_binding_class(NSS[u'cml'], "model", cellml_model)
    binder.set_binding_class(NSS[u'cml'], "component", cellml_component)
    binder.set_binding_class(NSS[u'cml'], "variable", cellml_variable)
    binder.set_binding_class(NSS[u'cml'], "units", cellml_units)
    binder.set_binding_class(NSS[u'cml'], "unit", cellml_unit)
    for mathml_elt in ['math', 'degree', 'logbase', 'otherwise',
                       'diff', 'plus', 'minus', 'times', 'divide',
                       'exp', 'ln', 'log', 'abs', 'power', 'root',
                       'leq', 'geq', 'lt', 'gt', 'eq', 'neq',
                       'rem',
                       'ci', 'cn', 'apply', 'piecewise', 'piece',
                       'sin', 'cos', 'tan', 'arcsin', 'arccos', 'arctan',
                       'csymbol']:
        exec "binder.set_binding_class(NSS[u'm'], '%s', mathml_%s)" % (mathml_elt, mathml_elt)
    binder.set_binding_class(NSS[u'm'], "and_", mathml_and)
    binder.set_binding_class(NSS[u'm'], "or_", mathml_or)
    return binder

def amara_parse_cellml(source, uri=None, prefixes=None):
    """Parse a CellML source with default rules and bindings."""
    binder = make_xml_binder()
    rules = [bt.ws_strip_element_rule(u'*')]
    return amara_parse(source, rules=rules, binderobj=binder)

def check_append_safety(elt):
    """
    Check whether elt is safe to make a child, i.e. that it isn't
    already a child elsewhere.
    """
    assert getattr(elt, 'next_elem', None) is None
    parent = getattr(elt, 'xml_parent', None)
    if parent:
        assert elt not in parent.xml_children

class element_base(amara.bindery.element_base):
    """
    Base element class to allow me to set certain attributes on my instances
    that are Python objects rather than unicode strings.
    """
    def __init__(self):
        self.xml_attributes = {} # Amara should really do this!
        super(element_base, self).__init__()
    
    def __delattr__(self, key):
        """
        Bypass Amara's __delattr__ for attribute names that start with _cml_
        """
        if key.startswith('_cml_'):
            del self.__dict__[key]
        else:
            amara.bindery.element_base.__delattr__(self, key)
    
    def __setattr__(self, key, value):
        """
        Bypass Amara's __setattr__ for attribute names that start with _cml_
        """
        if key.startswith('_cml_'):
            self.__dict__[key] = value
        else:
            amara.bindery.element_base.__setattr__(self, key, value)
    
    @property
    def rootNode(self):
        p = self.parentNode
        if p:
            return p.rootNode
        elif isinstance(self, mathml):
            raise ValueError('MathML element with no parent!')
        return self

    @property
    def cmeta_id(self):
        """Get the value of the cmeta:id attribute, or the empty string if not set."""
        return self.getAttributeNS(NSS['cmeta'], u'id')

    def xml_remove_child_at(self, index=-1):
        """
        Remove child object at a given index
        index - optional, 0-based index of child to remove (defaults to the last child)
        """
        obj = self.xml_children[index]
        if isinstance(obj, unicode):
            del self.xml_children[index]
        else:
            # Remove references to the object
            # Probably a slow way to go about this
            for attr, val in self.__dict__.items():
                if not (attr.startswith('xml') or
                        attr.startswith('_cml_') or
                        attr in self.xml_ignore_members):
                    next = getattr(val, 'next_elem', None)
                    if val == obj:
                        del self.__dict__[attr]
                        if next: self.__dict__[attr] = next
                    while next:
                        prev, val = val, next
                        next = getattr(val, 'next_elem', None)
                        if val == obj:
                            prev.next_elem = next
                            break
            del self.xml_children[index]

    def xml_doc(self):
        msg = []
        xml_attrs = []
        if hasattr(self, 'xml_attributes'):
            msg.append('Object references based on XML attributes:')
            for apyname in self.xml_attributes:
                local, ns = self.xml_attributes[apyname]
                if ns:
                    source_phrase = " based on '{%s}%s' in XML"%(ns, local)
                else:
                    source_phrase = " based on '%s' in XML"%(local)
                msg.append(apyname+source_phrase)
                xml_attrs.append(apyname)
        msg.append('Object references based on XML child elements:')
        for attr, val in self.__dict__.items():
            if not (attr.startswith('xml') or
                    attr.startswith('_cml_') or
                    attr in self.xml_ignore_members):
                if attr not in xml_attrs:
                    count = len(list(getattr(self, attr)))
                    if count == 1:
                        count_phrase = " (%s element)"%count
                    else:
                        count_phrase = " (%s elements)"%count
                    local, ns = val.localName, val.namespaceURI
                    if ns:
                        source_phrase = " based on '{%s}%s' in XML"%(ns, local)
                    else:
                        source_phrase = " based on '%s' in XML"%(local)
                    msg.append(attr+count_phrase+source_phrase)
        return u'\n'.join(msg)
    
    @property
    def xml_properties(self):
        """
        Return a dictionary whose keys are Python properties on this
        object that represent XML attributes and elements, and whose vaues
        are the corresponding objects (a subset of __dict__)
        """
        properties = {}
        for attr in self.__dict__:
            if (not (attr.startswith('xml')
                     or attr.startswith('_cml_')
                     or attr in self.xml_ignore_members)):
                properties[attr] = self.__dict__[attr]
        return properties

# Add some improved/new methods to all bindings
def add_methods_to_amara():
    def getAttributeNS(self, ns, local, default=u""):
        """
        Get the value of an attribute specified by namespace and localname.
    
        Optionally can also pass a default value if the attribute
        doesn't exist (defaults to the empty string).
        """
        attrs = getattr(self, 'xml_attributes', {})
        keys = [ (ns_, SplitQName(qname)[1]) 
                   for _, (qname, ns_) in attrs.items() ]
        values = [ unicode(getattr(self, attr))
                   for attr, (qname, ns_) in attrs.items() ]
        attr_dict = dict(zip(keys, values))
        return attr_dict.get((ns, local), default)
    
    def xml_element_children(self, elt=None):
        """Return an iterable over child elements of this element."""
        if elt is None:
            elt = self
        for child in elt.xml_children:
            if getattr(child, 'nodeType', None) == Node.ELEMENT_NODE:
                yield child
    
    def safe_remove_child(self, child, parent=None):
        """Remove a child element from parent in such a way that it can safely be added elsewhere."""
        if parent is None: parent = self
        parent.xml_remove_child(child)
        child.next_elem = None
        
    def replace_child(self, old, new, parent=None):
        """Replace child old of parent with new."""
        if parent is None: parent = self
        parent.xml_insert_after(old, new)
        self.safe_remove_child(old, parent)
    
    import new
    for method in ['getAttributeNS', 'xml_element_children', 'safe_remove_child', 'replace_child']:
        meth = new.instancemethod(locals()[method], None, amara.bindery.element_base)
        setattr(amara.bindery.element_base, method, meth)

add_methods_to_amara()

class comment_base(amara.bindery.comment_base):
    """An iterable version of comment nodes."""
    def __init__(self, body=None):
        amara.bindery.comment_base.__init__(self, body)
    def __iter__(self):
        return unitary_iterator(self)

######################################################################
#                          CellML elements                           #
######################################################################

class cellml_model(element_base):
    """
    Specialised class for the model element of a CellML document.
    Adds methods for collecting and reporting validation errors, etc.
    """
    
    def __init__(self):
        element_base.__init__(self)
        self._cml_validation_errors = []
        self._cml_validation_warnings = []
        self._cml_variables = {}
        self._cml_components = {}
        self._cml_units = {}
        self._cml_standard_units = {}
        self._cml_units_map = {}
        # Topologically sorted assignments list
        self._cml_assignments = []
    
    def __del__(self):
        self.clean_up()
        
    def clean_up(self):
        """Try to get the RDF library to clean up nicely."""
        cellml_metadata.remove_model(self)
        
    def get_config(self, config_attr=None):
        """Get the configuration store if it exists, or an attribute thereof."""
        config = getattr(self.xml_parent, '_cml_config', None)
        if config_attr:
            config = getattr(config, config_attr, None)
        return config
    
    def get_option(self, option_name):
        """Get the value of a command-line option."""
        config = getattr(self.xml_parent, '_cml_config', None)
        return config and getattr(config.options, option_name)

    def get_component_by_name(self, compname):
        """Return the component object that has name `compname'."""
        return self._cml_components[compname]
    
    def get_variable_by_name(self, compname, varname):
        """
        Return the variable object with name `varname' in component
        `compname'.
        """
        try:
            var = self._cml_variables[(compname, varname)]
        except KeyError, e:
            if compname == u'':
                if self.component.ignore_component_name:
                    compname = self.component.name
                else:
                    try:
                        compname, varname = cellml_variable.split_name(varname)
                    except ValueError:
                        raise e
                var = self._cml_variables[(compname, varname)]
            else:
                raise e
        return var
    
    def get_variable_by_oxmeta_name(self, name, throw=True):
        """
        Get the unique variable in this model with the given Oxford metadata
        name annotation.
        
        If throw is True, will raise ValueError if there is no such variable,
        or more than 1 match.  If throw is False, returns None in these cases.
        """
        vars = cellml_metadata.find_variables(self,
                                             ('bqbiol:is', NSS['bqbiol']),
                                             ('oxmeta:'+str(name), NSS['oxmeta']))
        if len(vars) == 1:
            var = vars[0]
        elif throw:
            raise ValueError('"%s" does not name a unique variable (matches: %s)'
                             % (name, str(vars)))
        else:
            var = None
        return var
    
    def get_variables_by_ontology_term(self, term, transitive=False):
        """Return a list of variables annotated with the given ontology term.
        
        The annotations have the same form as for oxmeta name annotations (see get_variable_by_oxmeta_name).
        However, here we are not restricted to namespace, and no check is done on the number of results returned.
        
        :param term: must be a (prefixed_name, nsuri) tuple.
        :param transitive: if True, look not just for direct annotations but also for terms belonging to
        the class given by prefixed_name, searching transitively along rdf:type predicates.
        """
        assert isinstance(term, tuple)
        assert len(term) == 2
        vars = []
        if transitive:
            terms = cellml_metadata.transitive_subjects(term)
        else:
            terms = [term]
        for term in terms:
            vars.extend(cellml_metadata.find_variables(self, ('bqbiol:is', NSS['bqbiol']), term))
            vars.extend(cellml_metadata.find_variables(self, ('bqbiol:isVersionOf', NSS['bqbiol']), term))
        return vars
    
    def get_variable_by_cmeta_id(self, cmeta_id):
        """
        Get the unique variable in this model with the given cmeta:id attribute value.
        """
        vars = self.xml_xpath(u'cml:component/cml:variable[@cmeta:id="%s"]' % cmeta_id)
        if len(vars) != 1:
            raise ValueError('"%s" does not ID a unique variable (matches: %s)'
                             % (cmeta_id, str(vars)))
        return vars[0]
    
    def get_all_variables(self):
        """Return an iterator over the variables in the model."""
        for comp in getattr(self, u'component', []):
            for var in getattr(comp, u'variable', []):
                yield var

    def _add_variable(self, var, varname, compname):
        """Add a new variable to the model."""
        assert (compname, varname) not in self._cml_variables
        self._cml_variables[(compname, varname)] = var

    def _del_variable(self, varname, compname):
        """Remove a variable from the model."""
        del self._cml_variables[(compname, varname)]

    def _add_component(self, comp, special=False):
        """Add a new component to the model."""
        if special:
            comp.xml_parent = self
        else:
            self.xml_append(comp)
        self._cml_components[comp.name] = comp

    def _del_component(self, comp):
        """Remove the given component from the model."""
        self.xml_remove_child(comp)
        del self._cml_components[comp.name]

    def validation_error(self, errmsg, level=logging.ERROR):
        """Log a validation error message.
        
        Message should be a unicode string.
        """
        self._cml_validation_errors.append(errmsg)
        logging.getLogger('validator').log(level, errmsg.encode('UTF-8'))
    def get_validation_errors(self):
        """Return the list of all errors found (so far) while validating this model."""
        return self._cml_validation_errors
    def validation_warning(self, errmsg, level=logging.WARNING):
        """Log a validation warning message.
        
        Message should be a unicode string.
        """
        self._cml_validation_warnings.append(errmsg)
        logging.getLogger('validator').log(level, errmsg.encode('UTF-8'))
    def get_validation_warnings(self):
        """Return the list of all warnings found (so far) while validating this model.
        """
        return self._cml_validation_warnings
    def _report_exception(self, e, show_xml_context):
        """Report an exception e as a validation error or warning.
        
        If show_xml_context is True, display the XML of the context
        of the exception as well.
        """
        e.show_xml_context = show_xml_context
        if e.warn:
            self.validation_warning(unicode(e), level=e.level)
        else:
            self.validation_error(unicode(e), level=e.level)

    def validate(self, xml_context=False,
                 invalid_if_warnings=False,
                 warn_on_units_errors=False,
                 check_for_units_conversions=False,
                 assume_valid=False, **ignored_kwargs):
        """Validate this model.

        Assumes that RELAX NG validation has been done.  Checks rules
        3.4.2.2, 3.4.3.2, 3.4.3.3, 3.4.5.2, 3.4.5.3, 3.4.5.4, 3.4.6.2, 3.4.6.3, 3.4.6.4,
        4.4.2.1, 4.4.3.2, 4.4.4, 5.4.1.2, 5.4.2.2, 6.4.2.5, 6.4.3.2, and 6.4.3.3
        in the CellML 1.0 spec, and performs units checking.

        Note that if some checks fail, most of the remaining checks
        will not be performed.  Hence when testing a model validate
        repeatedly until it passes.

        If xml_context is True, then the failing MathML tree will be
        displayed with every units error.

        If check_for_units_conversions is True, then generate a warning if
        units conversions will be needed.
        
        If assume_valid is True then fewer checks will be done - only
        what is required to set up the data structures needed for model
        transformation.

        Returns True iff the model validates.
        When invalid_if_warnings is True the model will fail to validate
        if there are any warnings, as well as if there are any errors.
        """
        self._validate_component_hierarchies()

        # Rule 5.4.2.2: units definitions may not be circular.
        # Also checks 5.4.1.2: no duplicate units names.
        if not assume_valid:
            for unit in self.get_all_units():
                self._check_unit_cycles(unit)
            DEBUG('validator', 'Checked for units cycles')
            # Check rule 3.4.3.3 too.
            self._check_variable_units_exist()

        if not self._cml_validation_errors:
            self._check_variable_mappings() # This sets up source variable links
        if not self._cml_validation_errors and not assume_valid:
            self._check_connection_units(check_for_units_conversions)

        # Rules 4.4.2.1 and 4.4.3.2: check name references in mathematics
        if not self._cml_validation_errors:
            assignment_exprs = self.search_for_assignments()
            if not assume_valid:
                for expr in assignment_exprs:
                    self._check_maths_name_references(expr, xml_context)

        # Rule 4.4.4: mathematical expressions may only modify
        # variables belonging to the current component.
        if not self._cml_validation_errors and not assume_valid:
            self._check_assigned_vars(assignment_exprs, xml_context)

        # Warn if mathematics outside the CellML subset is used.
        if not self._cml_validation_errors and not assume_valid:
            math_elts = self.xml_xpath(self.math_xpath_1 + u' | ' + self.math_xpath_2)
            self._check_cellml_subset(math_elts)

        # Classify variables and check for circular equations.
        # Does a topological sort of all equations in the process.
        # TODO: Handle reactions properly.
        if not self._cml_validation_errors:
            self._classify_variables(assignment_exprs, xml_context)
            self._order_variables(assignment_exprs, xml_context)

        # Appendix C.3.6: Equation dimension checking.
        if not self._cml_validation_errors and (not assume_valid or check_for_units_conversions):
            self._check_dimensional_consistency(assignment_exprs,
                                                xml_context,
                                                warn_on_units_errors,
                                                check_for_units_conversions)
        
        # Warn if unknown namespaces are used, just in case.
        unknown_nss = set(self.rootNode.xml_namespaces.keys()).difference(set(NSS.values()))
        if unknown_nss:
            self.validation_warning(u'Unrecognised namespaces used:\n  ' +
                                    u'\n  '.join(list(unknown_nss)))

        # Return validation result
        return not self._cml_validation_errors and \
               (not invalid_if_warnings or not self._cml_validation_warnings)

    def _validate_component_hierarchies(self):
        """Check Rule 6.4.3.2 (4): hierarchies must not be circular.
        
        Builds all the hierarchies, and checks for cycles.
        In the process, we also check the other rules in 6.4.3, and 6.4.2.5.
        """
        # First, we find the hierarchies that are defined.
        hiers = set()
        rels = []
        for group in getattr(self, u'group', []):
            local_hiers = set()
            for rel in getattr(group, u'relationship_ref', []):
                rels.append(rel)
                reln = rel.relationship
                ns = rel.xml_attributes[u'relationship'][1]
                name = getattr(rel, u'name', None)
                hier = (reln, ns, name)
                if hier in local_hiers:
                    self.validation_error("A group element must not contain two or more relationship_ref"
                                          " elements that define a relationship attribute in a common"
                                          " namespace with the same value and that have the same name"
                                          " attribute value (which may be non-existent) (6.4.2.5)."
                                          " Relationship '%s' name '%s' in namespace '%s' is repeated."
                                          % (reln, name or '', ns))
                local_hiers.add(hier)
                hiers.add(hier)
        # Now build & check each hierarchy
        for hier in hiers:
            self.build_component_hierarchy(hier[0], hier[1], hier[2], rels=rels)
        DEBUG('validator', 'Checked component hierachies')
    
    def _check_variable_mappings(self):
        """Check Rules 3.4.{5,6}: check variable mappings and interfaces are sane."""
        # First check connection elements and build mappings dict
        self.build_name_dictionaries()
        connected_components = set()
        for connection in getattr(self, u'connection', []):
            comps = frozenset([connection.map_components.component_1, connection.map_components.component_2])
            if comps in connected_components:
                self.validation_error("Each map_components element must map a unique pair of components "
                                      "(3.4.5.4). The pair ('%s', '%s') is repeated." % tuple(comps))
            connected_components.add(comps)
            self._validate_connection(connection)
        # Now check for variables that should receive a value but don't
        for comp in getattr(self, u'component', []):
            for var in getattr(comp, u'variable', []):
                for iface in [u'private_interface', u'public_interface']:
                    if getattr(var, iface, u'none') == u'in':
                        try:
                            var.get_source_variable()
                        except TypeError:
                            # No source variable found
                            self.validation_error("Variable '%s' has a %s attribute with value 'in', "
                                                  "but no component exports a value to that variable."
                                                  % (var.fullname(), iface))
        DEBUG('validator', 'Checked variable mappings')
    
    def _validate_connection(self, conn):
        """Validate the given connection element.
        
        Check that the given connection object defines valid mappings
        between variables, according to rules 3.4.5 and 3.4.6.
        """
        # Check we are allowed to connect these components
        try:
            comp1 = self.get_component_by_name(conn.map_components.component_1)
        except KeyError:
            self.validation_error("Connections must be between components defined in the current model "
                                  "(3.4.5.2). There is no component '%s'." % conn.map_components.component_1)
            return
        try:
            comp2 = self.get_component_by_name(conn.map_components.component_2)
        except KeyError:
            self.validation_error("Connections must be between components defined in the current model "
                                  "(3.4.5.3). There is no component '%s'." % conn.map_components.component_2)
            return
        if comp1 is comp2:
            self.validation_error("A connection must link two different components (3.4.5.4). "
                                  "The component '%s' is being connected to itself." % comp1.name)
            return
        # Get the parent of each component in the encapsulation hierarchy
        par1, par2 = comp1.parent(), comp2.parent()
        # The two components must either be siblings (maybe top-level) or parent & child.
        if not (par1 == comp2 or par2 == comp1 or par1 == par2):
            self.validation_error(u' '.join([
                'Connections are only permissible between sibling',
                'components, or where one is the parent of the other.\n',
                comp1.name,'and',comp2.name,'are unrelated.']))
            return
        # Now check each variable mapping
        for mapping in conn.map_variables:
            try:
                var1 = self.get_variable_by_name(comp1.name, mapping.variable_1)
            except KeyError:
                self.validation_error("A variable mapping must be between existing variables (3.4.6.2). "
                                      "Variable '%s' doesn't exist in component '%s'."
                                       % (mapping.variable_1, comp1.name))
                continue
            try:
                var2 = self.get_variable_by_name(comp2.name, mapping.variable_2)
            except KeyError:
                self.validation_error("A variable mapping must be between existing variables (3.4.6.2). "
                                      "Variable '%s' doesn't exist in component '%s'."
                                       % (mapping.variable_2, comp2.name))
                continue
            errm, e = ['Interface mismatch mapping',var1.fullname(),'and',var2.fullname(),':\n'], None
            if par1 == par2:
                # Siblings, so should have differing public interfaces
                if not hasattr(var1, 'public_interface'):
                    e = 'missing public_interface attribute on ' + \
                        var1.fullname() + '.'
                elif not hasattr(var2, 'public_interface'):
                    e = 'missing public_interface attribute on ' + \
                        var2.fullname() + '.'
                elif var1.public_interface == var2.public_interface:
                    e = 'public_interface attributes are identical.'
                else:
                    if var1.public_interface == 'in':
                        var1._set_source_variable(var2)
                    else:
                        var2._set_source_variable(var1)
            else:
                if par2 == comp1:
                    # Component 1 is the parent of component 2
                    var1, var2 = var2, var1
                # Now var2 is in the parent component, and var1 in the child
                if not hasattr(var1, 'public_interface'):
                    e = var1.fullname()+' missing public_interface.'
                elif not hasattr(var2, 'private_interface'):
                    e = var2.fullname()+' missing private_interface.'
                elif var1.public_interface == var2.private_interface:
                    e = 'relevant interfaces have identical values.'
                else:
                    if var1.public_interface == 'in':
                        var1._set_source_variable(var2)
                    else:
                        var2._set_source_variable(var1)
            # If there was an error, log it
            if e:
                errm.append(e)
                self.validation_error(u' '.join(errm))
    
    def _check_variable_units_exist(self):
        """Check rule 3.4.3.3: that the units declared for variables exist."""
        for var in self.get_all_variables():
            try:
                var.get_units()
            except KeyError:
                self.validation_error("The value of the units attribute on a variable must be either "
                                      "one of the standard units or the name of a unit defined in the "
                                      "current component or model (3.4.3.3). The units '%s' on the "
                                      "variable '%s' in component '%s' do not exist."
                                      % (var.units, var.name, var.component.name))

    def _check_connection_units(self, check_for_units_conversions=False):
        """Check that the units of mapped variables are dimensionally consistent.
        
        If check_for_units_conversions is True we also warn if they are not equivalent,
        since much processing software may not be able to handle that case.
        """
        for conn in getattr(self, u'connection', []):
            comp1 = self.get_component_by_name(conn.map_components.component_1)
            comp2 = self.get_component_by_name(conn.map_components.component_2)
            for mapping in conn.map_variables:
                var1 = self.get_variable_by_name(comp1.name, mapping.variable_1)
                var2 = self.get_variable_by_name(comp2.name, mapping.variable_2)
                # Check the units
                u1 = var1.get_units()
                u2 = var2.get_units()
                if not u1 == u2:
                    if not u1.dimensionally_equivalent(u2):
                        self.validation_error(u' '.join([
                            var1.fullname(), 'and', var2.fullname(),
                            'are mapped, but have dimensionally inconsistent units.']))
                    elif check_for_units_conversions:
                        self.validation_warning(
                            u' '.join(['Warning: mapping between', var1.fullname(), 'and',
                                       var2.fullname(), 'will require a units conversion.']),
                            level=logging.WARNING_TRANSLATE_ERROR)

    def _check_maths_name_references(self, expr, xml_context=False):
        """Check rules 4.4.2.1 and 4.4.3.2: name references in mathematics."""
        if isinstance(expr, mathml_ci):
            # Check variable exists
            try:
                _ = expr.variable
            except KeyError:
                self._report_exception(
                        MathsError(expr,
                                   "The content of a MathML ci element must match the name of a variable "
                                   "in the enclosing component, once whitespace normalisation has been "
                                   "performed (4.4.2.1). Variable '%s' does not exist in component '%s'."
                                   % (unicode(expr).strip(), expr.component.name)),
                        xml_context)
        elif isinstance(expr, mathml_cn):
            # Check units exist
            try:
                expr.get_units()
            except KeyError:
                self._report_exception(
                        MathsError(expr,
                                   "Units on a cn element must be standard or defined in the current "
                                   "component or model (4.4.3.2). Units '%s' are not defined in "
                                   "component '%s'." % (expr.units, expr.component.name)),
                        xml_context)
        else:
            # Recurse
            for child in expr.xml_element_children():
                self._check_maths_name_references(child, xml_context)

    def _check_assigned_vars(self, assignments, xml_context=False):
        """Check Rule 4.4.4: mathematical expressions may only modify
        variables belonging to the current component.
        """
        for expr in assignments:
            try:
                expr.check_assigned_var()
            except MathsError, e:
                self._report_exception(e, xml_context)
        DEBUG('validator', 'Checked variable assignments')

    def _check_cellml_subset(self, math_elts, root=True):
        """Warn if MathML outside the CellML subset is used."""
        for elt in math_elts:
            if not elt.localName in CELLML_SUBSET_ELTS and \
                   elt.namespaceURI == NSS[u'm']:
                self.validation_warning(u' '.join([
                    u'MathML element', elt.localName,
                    u'is not in the CellML subset.',
                    u'Some tools may not be able to process it.']))
            self._check_cellml_subset(self.xml_element_children(elt), False)
        if root:
            DEBUG('validator', 'Checked for CellML subset')

    def _classify_variables(self, assignment_exprs, xml_context=False):
        """Determine the type of each variable.
        
        Note that mapped vars must have already been classified by
        self._check_variable_mappings, and the RELAX NG schema ensures
        that a variable cannot be both Mapped and MaybeConstant.
        
        Builds the equation dependency graph in the process.
        """
        # Classify those vars that might be constants,
        # i.e. have an initial value assigned.
        for var in self.get_all_variables():
            if hasattr(var, u'initial_value'):
                var._set_type(VarTypes.MaybeConstant)
        # Now classify by checking usage in mathematics, building
        # an equation dependency graph in the process.
        for expr in assignment_exprs:
            try:
                expr.classify_variables(root=True)
            except MathsError, e:
                self._report_exception(e, xml_context)
        # Unused vars still classified as MaybeConstant are constants
        for var in self.get_all_variables():
            if var.get_type() == VarTypes.MaybeConstant:
                var._set_type(VarTypes.Constant)
        DEBUG('validator', 'Classified variables')
            
    def _order_variables(self, assignment_exprs, xml_context=False):
        """Topologically sort the equation dependency graph.
        
        This orders all the assignment expressions in the model, to
        allow procedural code generation.  It also checks that equations
        are not cyclic (we don't support DAEs).
        """
        self.clear_assignments() # Ensure we start from a clean slate
        try:
            self._cml_sorting_variables_stack = []
            for var in self.get_all_variables():
                if var.get_colour() == DFS.White:
                    self.topological_sort(var)
            for expr in assignment_exprs:
                if expr.get_colour() == DFS.White:
                    self.topological_sort(expr)
        except MathsError, e:
            self._report_exception(e, xml_context)
        DEBUG('validator', 'Topologically sorted variables')

    math_xpath_1 = u'cml:component/m:math'
    math_xpath_2 = u'cml:component/cml:reaction/cml:variable_ref/cml:role/m:math'
    apply_xpath_1 = u'/m:apply[m:eq]'
    apply_xpath_2 = u'/m:semantics/m:apply[m:eq]'
    
    def search_for_assignments(self, comp=None):
        """Search for assignment expressions in the model's mathematics.
        
        If comp is supplied, will only return assignments in that component.
        """
        assignments_xpath = u' | '.join([self.math_xpath_1 + self.apply_xpath_1,
                                         self.math_xpath_1 + self.apply_xpath_2,
                                         self.math_xpath_2 + self.apply_xpath_1,
                                         self.math_xpath_2 + self.apply_xpath_2])
        if comp:
            assignments_xpath = assignments_xpath.replace(u'component',
                                                          u'component[@name="%s"]' % comp.name)
        return self.xml_xpath(assignments_xpath)
    
    def build_name_dictionaries(self, rebuild=False):
        """
        Create dictionaries mapping names of variables and components to
        the objects representing them.
        Dictionary keys for variables will be
          (component_name, variable_name).
        
        If rebuild is True, clear the dictionaries first.
        """
        if rebuild:
            self._cml_variables = {}
            self._cml_components = {}
        if not self._cml_components:
            for comp in getattr(self, u'component', []):
                if comp.name in self._cml_components:
                    self.validation_error("Component names must be unique within a model (3.4.2.2)."
                                          " The name '%s' is repeated." % comp.name)
                self._cml_components[comp.name] = comp
                for var in getattr(comp, u'variable', []):
                    key = (comp.name, var.name)
                    if key in self._cml_variables:
                        self.validation_error("Variable names must be unique within a component (3.4.3.2)."
                                              " The name '%s' is repeated in component '%s'."
                                              % (var.name, comp.name))
                    self._cml_variables[key] = var

    def build_component_hierarchy(self, relationship, namespace=None, name=None, rels=None):
        """
        Create all the parent-child links for the given component hierarchy.

        relationship gives the type of the hierarchy. If it is not one of the
        CellML types (i.e. encapsulation or containment) then the namespace URI
        must be specified.  Multiple non-encapsulation hierarchies of the same
        type can be specified by giving the name argument.
        """
        key = (relationship, namespace, name)
        # Set all components to have no parent or children, under this hierarchy
        for comp in getattr(self, u'component', []):
            comp._clear_hierarchy(key)
        self.build_name_dictionaries()
        # Find nodes defining this hierarchy
        if rels is None:
            rels = self.xml_xpath(u'cml:group/cml:relationship_ref')
        groups = []
        for rel in rels:
            # NB: The first test below relies on there only being one
            # attr with localname of 'relationship' on each rel.
            # So let's check this...
            if hasattr(rel, u'relationship_'):
                self.validation_error(u' '.join([
                    'relationship_ref element has multiple relationship',
                    'attributes in different namespaces:\n'] + 
                        map(lambda qn,ns: '('+qn+','+ns+')',
                            rel.xml_attributes.values())))
                return
            if rel.relationship == relationship and \
               rel.xml_attributes[u'relationship'][1] == namespace and \
               getattr(rel, u'name', None) == name:
                # It's in the hierarchy
                groups.append(rel.xml_parent)
        # Now build all the parent links for this hierarchy
        def set_parent(p, crefs):
            for cref in crefs:
                # Find component cref refers to
                try:
                    c = self.get_component_by_name(cref.component)
                except KeyError:
                    self.validation_error("A component_ref element must reference a component in the current"
                                          " model (6.4.3.3). Component '%s' does not exist." % cref.component)
                    return
                # Set c's parent to p
                c._set_parent_component(key, p)
                if hasattr(cref, 'component_ref'):
                    # If we already have children under this hierarchy it's a validation error
                    if c._has_child_components(key):
                        self.validation_error("In a given hierarchy, only one component_ref element "
                                              "referencing a given component may contain children (6.4.3.2)."
                                              " Component '%s' has children in multiple locations." % c.name)
                    # Set parent of c's children to c
                    set_parent(c, cref.component_ref)
                elif p is None and namespace is None:
                    self.validation_error("Containment and encapsulation relationships must be hierarchical"
                                          " (6.4.3.2). Potentially top-level component '%s' has not been"
                                          " given children." % cref.component)
        for group in groups:
            set_parent(None, group.component_ref)
        if self._cml_validation_errors:
            return
        # Check for a cycle in the hierarchy (rule 6.4.3.2 (4)).
        # Note that since we have already ensured that no component is a parent in more than one location,
        # nor is any component a child more than once, so the only possibility for a cycle is if one of
        # the components referenced as a child of a group element is also referenced as a (leaf) descendent
        # of one of its children.  We check for this by following parent links backwards.
        def has_cycle(root_comp, cur_comp):
            if cur_comp is None:
                return False
            elif cur_comp is root_comp:
                return True
            else:
                return has_cycle(root_comp, cur_comp.parent(reln_key=key))
        for group in groups:
            for cref in group.component_ref:
                # Find component cref refers to
                c = self.get_component_by_name(cref.component)
                if has_cycle(c, c.parent(reln_key = key)):
                    n, ns = name or "", namespace or ""
                    self.validation_error("The '%s' relationship hierarchy with name '%s' and namespace"
                                          " '%s' has a cycle" % (relationship, n, ns))
        return

    def topological_sort(self, node):
        """
        Do a topological sort of all assignment expressions and variables
        in the model.

        node should be an expression or variable object that inherits from
        Colourable and has methods get_dependencies, get_component
        """
        node.set_colour(DFS.Gray)
        # Keep track of gray variables, for reporting cycles
        if isinstance(node, cellml_variable):
            self._cml_sorting_variables_stack.append(node.fullname())
        elif node.is_ode():
            # This is an expression defining an ODE; the name on the
            # stack will look something like d(V)/d(t)
            n1, n2 = map(lambda v: v.fullname(), node.assigned_variable())
            self._cml_sorting_variables_stack.append(u'd'+n1+u'/d'+n2)
        # Visit children in the dependency graph
        for dep in node.get_dependencies():
            if type(dep) == types.TupleType:
                # This is an ODE dependency, so get the defining expression
                dep = dep[0].get_ode_dependency(dep[1], node)
            if dep.get_colour() == DFS.White:
                self.topological_sort(dep)
            elif dep.get_colour() == DFS.Gray:
                # We have a cyclic dependency
                if isinstance(dep, cellml_variable):
                    i = self._cml_sorting_variables_stack.index(dep.fullname())
                elif node.is_ode():
                    n1, n2 = map(lambda v: v.fullname(),
                                 dep.assigned_variable())
                    i = self._cml_sorting_variables_stack.index(
                        u'd'+n1+u'/d'+n2)
                else:
                    # Since any variable always depends on a mathematical
                    # expression defining it, the only case where the
                    # expression is gray before the corresponding variable
                    # (apart from an ODE, dealt with above) is if a tree
                    # started at an expression.  Hence the whole stack
                    # is in the cycle.
                    i = 0
                varnames = self._cml_sorting_variables_stack[i:]
                self.validation_error(u' '.join([
                    u'There is a cyclic dependency involving the following',
                    u'variables:', u','.join(varnames)]))
        # Finish this node, and add it to the appropriate sorted list
        node.set_colour(DFS.Black)
        self._add_sorted_assignment(node)
        # Pop the gray variables stack
        if (isinstance(node, cellml_variable) or node.is_ode()):
            self._cml_sorting_variables_stack.pop()
        return
    
    def _add_sorted_assignment(self, a):
        """
        During the topological sort, add a finished assignment to the
        list.  This list can then be executed in order to simulate the
        model.

        The element added can either be a MathML expression
        representing an assignment, or a CellML variable element,
        indicating an assignment due to a variable mapping.
        """
        self._cml_assignments.append(a)
    def _remove_assignment(self, a):
        """Remove the given assignment from our list.

        This method is used by the partial evaluator."""
        self._cml_assignments.remove(a)
    def get_assignments(self):
        """
        Return a sorted list of all the assignments in the model.

        Assignments can either be instances of cellml_variable, in
        which case they represent a variable mapping, or instances of
        mathml_apply, representing top-level assignment expressions.
        """
        return self._cml_assignments
    def clear_assignments(self):
        """Clear the assignments list."""
        self._cml_assignments = []

    def do_binding_time_analysis(self):
        """Perform a binding time analysis on the model's mathematics.

        This requires variables to have been classified and a
        topological sort of the mathematics to have been performed.

        Variables and top-level expressions are processed in the order
        given by the topological sort, hence only a single pass is
        necessary.

        Variables are classified based on their type:
          State, Free -> dynamic
          Constant -> static
          Mapped -> binding time of source variable
          Computed -> binding time of defining expression

        Expressions are dealt with by recursively annotating
        subexpressions.  See code in the MathML classes for details.
        """
        for item in self.get_assignments():
            if isinstance(item, cellml_variable):
                # Set binding time based on type
                item._get_binding_time()
            else:
                # Compute binding time recursively
                item._get_binding_time()

    def _check_dimensional_consistency(self, assignment_exprs,
                                       xml_context=False,
                                       warn_on_units_errors=False,
                                       check_for_units_conversions=False):
        """Appendix C.3.6: Equation dimension checking."""
        self._cml_conversions_needed = False
        # Check dimensions
        for expr in assignment_exprs:
            try:
                expr.get_units()
            except UnitsError, e:
                if warn_on_units_errors:
                    e.warn = True
                    e.level = logging.WARNING
                self._report_exception(e, xml_context)
        # Check if units conversions will be needed
        if check_for_units_conversions and not self._cml_validation_errors:
            boolean = self.get_units_by_name('cellml:boolean')
            for expr in assignment_exprs:
                try:
                    DEBUG('validator', "Checking units in", element_xpath(expr), expr.component.name)
                    expr._set_in_units(boolean, no_act=True)
                except UnitsError:
                    pass
            # Warn if conversions used
            if self._cml_conversions_needed:
                self.validation_warning('The mathematics in this model require units conversions.',
                                        level=logging.WARNING)
        DEBUG('validator', 'Checked units')

    def _check_unit_cycles(self, unit):
        """Check for cyclic units definitions.

        We do this by doing a depth-first search from unit.
        """
        if unit.get_colour() != DFS.White:
            # Allow self.validate to call us without colour check
            return
        unit.set_colour(DFS.Gray)
        # Get the object unit is defined in
        parent = unit.xml_parent or self
        for u in getattr(unit, u'unit', []):
            # Explore units that this unit is defined in terms of
            try:
                v = parent.get_units_by_name(u.units)
            except KeyError:
                self.validation_error("The value of the units attribute on a unit element must be taken"
                                      " from the dictionary of standard units or be the name of a"
                                      " user-defined unit in the current component or model (5.4.2.2)."
                                      " Units '%s' are not defined." % u.units)
                continue
            if v.get_colour() == DFS.White:
                self._check_unit_cycles(v)
            elif v.get_colour() == DFS.Gray:
                # We have a cycle
                self.validation_error("Units %s and %s are in a cyclic units definition"
                                      % (unit.name, v.name))
        unit.set_colour(DFS.Black)

    def _build_units_dictionary(self):
        """Create a dictionary mapping units names to objects, for all units definitions in this element."""
        # Standard units
        std_units = self._cml_standard_units
        def make(name, bases):
            return cellml_units.create_new(self, name, bases, standard=True)
        # SI base units & dimensionless
        base_units = [u'ampere', u'candela',  u'dimensionless', u'kelvin',
                      u'kilogram', u'metre', u'mole', u'second']
        base_units.append(u'#FUDGE#') # Used for PE of naughty models
        for units in base_units:
            std_units[units] = make(units, [])
        # Special cellml:boolean units
        boolean = make(u'cellml:boolean', [])
        std_units[u'cellml:boolean'] = boolean
        # Convenience derived units
        gram = make('gram', [{'units': 'kilogram', 'multiplier': '0.001'}])
        litre = make('litre', [{'multiplier': '1000', 'prefix': 'centi',
                                'units': 'metre', 'exponent': '3'}])
        # SI derived units
        radian = make('radian', [{'units': 'metre'},
                                 {'units': 'metre', 'exponent': '-1'}])
        steradian = make('steradian', [{'units': 'metre', 'exponent': '2'},
                                       {'units': 'metre', 'exponent': '-2'}])
        hertz = make('hertz', [{'units': 'second', 'exponent': '-1'}])
        newton = make('newton', [{'units': 'metre'},
                                 {'units': 'kilogram'},
                                 {'units': 'second', 'exponent': '-2'}])
        pascal = make('pascal', [{'units': 'newton'},
                                 {'units': 'metre', 'exponent': '-2'}])
        joule = make('joule', [{'units': 'newton'},
                               {'units': 'metre'}])
        watt = make('watt', [{'units': 'joule'},
                             {'units': 'second', 'exponent': '-1'}])
        coulomb = make('coulomb', [{'units': 'second'},
                                   {'units': 'ampere'}])
        volt = make('volt', [{'units': 'watt'},
                             {'units': 'ampere', 'exponent': '-1'}])
        farad = make('farad', [{'units': 'coulomb'},
                               {'units': 'volt', 'exponent': '-1'}])
        ohm = make('ohm', [{'units': 'volt'},
                           {'units': 'ampere', 'exponent': '-1'}])
        siemens = make('siemens', [{'units': 'ampere'},
                                   {'units': 'volt', 'exponent': '-1'}])
        weber = make('weber', [{'units': 'volt'},
                               {'units': 'second'}])
        tesla = make('tesla', [{'units': 'weber'},
                               {'units': 'metre', 'exponent': '-2'}])
        henry = make('henry', [{'units': 'weber'},
                               {'units': 'ampere', 'exponent': '-1'}])
        celsius = make('celsius', [{'units': 'kelvin', 'offset': '-273.15'}])
        lumen = make('lumen', [{'units': 'candela'},
                               {'units': 'steradian'}])
        lux = make('lux', [{'units': 'lumen'},
                           {'units': 'metre', 'exponent': '-2'}])
        becquerel = make('becquerel', [{'units': 'second', 'exponent': '-1'}])
        gray = make('gray', [{'units': 'joule'},
                             {'units': 'kilogram', 'exponent': '-1'}])
        sievert = make('sievert', [{'units': 'joule'},
                                   {'units': 'kilogram', 'exponent': '-1'}])
        katal = make('katal', [{'units': 'second', 'exponent': '-1'},
                               {'units': 'mole'}])
        for units in [becquerel, celsius, coulomb, farad, gram, gray, henry,
                      hertz, joule, katal, litre, lumen, lux, newton, ohm,
                      pascal, radian, siemens, sievert, steradian, tesla,
                      volt, watt, weber]:
            std_units[units.name] = units
        # American spellings
        std_units[u'meter'] = std_units[u'metre']
        std_units[u'liter'] = std_units[u'litre']
        # User-defined units
        model_units = self._cml_units
        model_units.update(std_units)
        if hasattr(self, u'units'):
            for units in self.units:
                if units.name in model_units:
                    self.validation_error("Units names must be unique within the parent component or model,"
                                          " and must not redefine the standard units (5.4.1.2)."
                                          " The units definition named '%s' in the model is a duplicate." % units.name)
                model_units[units.name] = units
        # Update units hashmap
        for u in model_units.itervalues():
            self._add_units_obj(u)
    
    def get_standard_units(self):
        """Get a dictionary mapping the names of the standard CellML units to their definitions."""
        if not self._cml_standard_units:
            self._build_units_dictionary()
        return self._cml_standard_units

    def get_all_units(self):
        """Get a list of all units objects, including the standard units."""
        if not self._cml_units:
            self._build_units_dictionary()
        units = self._cml_units.values()
        for comp in getattr(self, u'component', []):
            units.extend(comp.get_all_units())
        return units

    def get_units_by_name(self, uname):
        """Return an object representing the element that defines the units named `uname'."""
        if not self._cml_units:
            self._build_units_dictionary()
        try:
            return self._cml_units[uname]
        except KeyError:
            raise KeyError("Units '%s' are not defined in the current component or model." % uname)

    def add_units(self, name, units):
        """Add an entry in our units dictionary for units named `name' with element object `units'."""
        if not self._cml_units:
            self._build_units_dictionary()
        assert name not in self._cml_units
        self._cml_units[name] = units
        self._add_units_obj(units)
        return

    def has_units(self, units):
        """Test whether a given units definition appears in the model."""
        return units in self._cml_units.itervalues()

    def _add_units_obj(self, units):
        """Add a units object into the global hashmap."""
        if not units.uniquify_tuple in self._cml_units_map:
            self._cml_units_map[units.uniquify_tuple] = units
        return
    
    def _get_units_obj(self, units):
        """Unique-ify this units object.

        If an object with the same definition already exists, return that.
        Otherwise return the given units object.

        'Same definition' is based on the cellml_units.uniquify_tuple
        property, which in turn is based partly on the generated name
        which would be given to these units, since that really *must*
        be unique in generated models.
        """
        return self._cml_units_map.get(units.uniquify_tuple, units)
    
    def _is_new_units_obj(self, units):
        """Have these units been generated already?

        i.e. is a units object with this definition in our map?
        """
        return units.uniquify_tuple not in self._cml_units_map

    def add_units_conversions(self):
        """Add explicit units conversion mathematics where necessary."""
        import_processors()
        processors.UnitsConverter(self).add_all_conversions()

    def find_state_vars(self):
        """Return a list of the state variable elements in this model."""
        state_vars = []
        for var in self.get_all_variables():
            if var.get_type() == VarTypes.State:
                state_vars.append(var)
        return state_vars

    def find_free_vars(self):
        """Return a list of the free variable elements in this model."""
        free_vars = []
        for var in self.get_all_variables():
            if var.get_type() == VarTypes.Free:
                free_vars.append(var)
        return free_vars
    
    def calculate_extended_dependencies(self, nodes, prune=[],
                                        prune_deps=[],
                                        state_vars_depend_on_odes=False,
                                        state_vars_examined=set()):
        """Calculate the extended dependencies of the given nodes.

        Recurse into the dependency graph, in order to construct a
        set, for each node in nodes, of all the nodes on which it
        depends, either directly or indirectly.

        Each node IS included in its own dependency set.

        If prune is specified, it should be a set of nodes for which
        we won't include their dependencies or the nodes themselves.
        This is useful e.g. for pruning variables required for calculating
        a stimulus if the stimulus is being provided by another method.
        prune_deps is similar: dependencies of these nodes will be excluded,
        but the nodes themselves will be included if asked for.
        
        If state_vars_depend_on_odes is True, then considers state variables
        to depend on the ODE defining them.

        Requires the dependency graph to be acyclic.

        Return the union of all the dependency sets.
        """
        deps = set()
        for node in nodes:
            if node in prune or (isinstance(node, mathml_apply) and
                                 isinstance(node.operator(), mathml_eq) and
                                 isinstance(node.eq.lhs, mathml_ci) and
                                 node.eq.lhs.variable in prune):
                continue
            if type(node) == tuple:
                # This is an ODE dependency, so get the defining expression instead.
                ode = True
                orig_node = node
                node = node[0].get_ode_dependency(node[1])
                if orig_node in prune_deps:
                    # Include the defining expression, but skip its dependencies
                    deps.add(node)
                    continue
                free_var = node.eq.lhs.diff.independent_variable
            else:
                ode = False
            deps.add(node)
            if node in prune_deps:
                # Skip dependencies of this node
                continue
            nodedeps = set(node.get_dependencies())
            if ode and not node._cml_ode_has_free_var_on_rhs:
                # ODEs depend on their independent variable.  However,
                # when writing out code we don't want to pull the free
                # variable in just for this, as the compiler then
                # gives us unused variable warnings.
                nodedeps.remove(free_var)
            if (state_vars_depend_on_odes and isinstance(node, cellml_variable)
                and node.get_type() == VarTypes.State
                and node not in state_vars_examined):
                nodedeps.update(node.get_all_expr_dependencies())
                state_vars_examined.add(node)
            deps.update(self.calculate_extended_dependencies(nodedeps,
                                                             prune=prune,
                                                             prune_deps=prune_deps,
                                                             state_vars_depend_on_odes=state_vars_depend_on_odes,
                                                             state_vars_examined=state_vars_examined))
        return deps

    def is_self_excitatory(self):
        """Determine whether this model is self-excitatory,
        i.e. does not require an external stimulus.
        """
        meta_id = self.cmeta_id
        if not meta_id:
            return False
        property = cellml_metadata.create_rdf_node(('pycml:is-self-excitatory', NSS['pycml']))
        source = cellml_metadata.create_rdf_node(fragment_id=meta_id)
        return cellml_metadata.get_target(self, source, property) == 'yes'

    def xml(self, stream=None, writer=None, **wargs):
        """Serialize back to XML.
        If stream is given, output to stream.
        If writer is given, use it directly.
        If neither a stream nor a writer is given, return the output text
        as a Python string (not Unicode) encoded as UTF-8.

        This overrides Amara's method, in order to force declaration of
        various namespaces with prefixes on this element, and to ensure
        the RDF annotations are up-to-date.

        See base class docs for possible keyword arguments.
        """
        extra_namespaces = {u'cellml': NSS[u'cml'],
                            u'pe': NSS[u'pe'],
                            u'lut': NSS[u'lut']}
        
        # Update RDF block if necessary
        cellml_metadata.update_serialized_rdf(self)
        
        temp_stream = None
        close_document = 0
        if not writer:
            #Change the default to *not* generating an XML decl
            if not wargs.get('omitXmlDeclaration'):
                wargs['omitXmlDeclaration'] = u'yes'
            if stream:
                writer = amara.bindery.create_writer(stream, wargs)
            else:
                temp_stream = StringIO()
                writer = amara.bindery.create_writer(temp_stream, wargs)

            writer.startDocument()
            close_document = 1
        writer.startElement(self.nodeName, self.namespaceURI,
                            extraNss=extra_namespaces)
        if hasattr(self, 'xml_attributes'):
            for apyname in self.xml_attributes:
                aqname, ans = self.xml_attributes[apyname]
                val = self.__dict__[apyname]
                writer.attribute(aqname, val, ans)
        for child in self.xml_children:
            if isinstance(child, unicode):
                writer.text(child)
            else:
                child.xml(writer=writer)
        writer.endElement(self.nodeName, self.namespaceURI)
        if close_document:
            writer.endDocument()
        return temp_stream and temp_stream.getvalue()
    

class cellml_component(element_base):
    """
    Specialised component class, with additional helper methods.
    """

    def __init__(self):
        element_base.__init__(self)
        self._cml_parents = {}
        self._cml_children = {}
        self._cml_units = None
        self._cml_created_by_pe = False
    
    @property
    def ignore_component_name(self):
        """Whether to not include the component name in the full names of contained variables."""
        return self._cml_created_by_pe or self.name == u''
    
    
    def parent(self, relationship=u'encapsulation', namespace=None, name=None, reln_key=None):
        """Find the parent of this component in the given hierarchy.
        
        We default to the encapsulation hierarchy.

        relationship gives the type of the hierarchy. If it is not one
        of the CellML types (i.e. encapsulation or containment) then
        the namespace URI must be specified.  Multiple non-encapsulation
        hierarchies of the same type can be specified by giving the name
        argument.

        Results are cached for efficiency.
        """
        key = reln_key or (relationship, namespace, name)
        if not key in self._cml_parents:
            assert(reln_key is None)
            self.xml_parent.build_component_hierarchy(relationship, namespace, name)
        return self._cml_parents[key]

    def _clear_hierarchy(self, reln_key):
        """Unset our parent & children in the given hierarchy."""
        self._cml_parents[reln_key] = None
        self._cml_children[reln_key] = []

    def _set_parent_component(self, reln_key, parent):
        """Set the parent of this component in the relationship hierarchy indexed by reln_key to parent.
        Trigger a validation error if we already have a parent in this hierarchy.
        Also add ourselves to parent's children.
        """
        if not reln_key in self._cml_parents or self._cml_parents[reln_key] is None:
            # Only set parent if we don't already have one
            self._cml_parents[reln_key] = parent
        else:
            self.xml_parent.validation_error("In a given hierarchy, a component may not be a child more"
                                             " than once (6.4.3.2). Component '%s' has multiple parents."
                                             % self.name)
        if not parent is None:
            parent._add_child_component(reln_key, self)
    
    def _add_child_component(self, reln_key, child):
        """Add child to our list of children in the relationship hierarchy indexed by reln_key."""
        if not reln_key in self._cml_children:
            self._cml_children[reln_key] = []
        self._cml_children[reln_key].append(child)
    
    def _has_child_components(self, reln_key):
        """Determine whether we have any children in the given relationship hierarchy."""
        return self._cml_children.get(reln_key, []) != []


    def _build_units_dictionary(self):
        """Create a dictionary mapping units names to objects, for all units definitions in this element."""
        self._cml_units = {}
        for units in getattr(self, u'units', []):
            if units.name in self._cml_units:
                self.validation_error("Units names must be unique within the parent component (5.4.1.2)."
                                      " The name '%s' in component '%s' is duplicated."
                                      % (units.name, self.name))
            try:
                if self.xml_parent.get_units_by_name(units.name).standard == u'yes':
                    self.validation_error("Units definitions must not redefine the standard units (5.4.1.2)."
                                          " The name '%s' in component '%s' is not allowed."
                                          % (units.name, self.name))
            except:
                pass
            self._cml_units[units.name] = units
            self.xml_parent._add_units_obj(units)
    def get_units_by_name(self, uname):
        """Return an object representing the element that defines the units named `uname'."""
        if self._cml_units is None:
            self._build_units_dictionary()
        if uname in self._cml_units:
            # Units are defined in this component
            return self._cml_units[uname]
        else:
            # Look up units in model element instead
            return self.xml_parent.get_units_by_name(uname)
    def add_units(self, name, units):
        """Add an entry in our units dictionary for units named `name' with element object `units'."""
        if self._cml_units is None:
            self._build_units_dictionary()
        self._cml_units[name] = units
        self.xml_parent._add_units_obj(units)
        return

    def get_all_units(self):
        """Get a list of all units objects defined in this component."""
        if self._cml_units is None:
            self._build_units_dictionary()
        return self._cml_units.values()

    def get_variable_by_name(self, varname):
        """Return the variable object with name `varname' in this component."""
        return self.xml_parent.get_variable_by_name(self.name, varname)

    def _add_variable(self, var):
        """Add a variable to this component."""
        # Add element
        self.xml_append(var)
        # Add to dictionary
        self.xml_parent._add_variable(var, var.name, self.name)
        return
    def _del_variable(self, var, keep_annotations=False):
        """Remove a variable from this component."""
        if not keep_annotations:
            # Remove metadata about the variable
            var.remove_rdf_annotations()
        # Remove the element
        self.xml_remove_child(var)
        # Remove from dictionary
        self.xml_parent._del_variable(var.name, self.name)
        return

    @staticmethod
    def create_new(elt, name):
        """Create a new component with the given name."""
        new_comp = elt.xml_create_element(u'component', NSS[u'cml'],
                                          attributes={u'name': unicode(name)})
        return new_comp


class cellml_variable(Colourable, element_base):
    """
    Class representing CellML <variable> elements.
    """
    def __init__(self):
        super(cellml_variable, self).__init__()
        self.clear_dependency_info()
        return
    
    def clear_dependency_info(self):
        """Clear the type, dependency, etc. information for this variable.
        
        This allows us to re-run the type & dependency analysis for the model.
        """
        # The type of this variable is not yet known
        self._cml_var_type = VarTypes.Unknown
        self._cml_source_var = None
        self._cml_value = {}
        self._cml_binding_time = None
        # Dependency graph edges
        self._cml_depends_on = []
        self._cml_depends_on_ode = {}
        self._cml_usage_count = 0
        self.clear_colour()

    def __hash__(self):
        """Hashing function for variables.

        Hash is based on hashing the full name of the variable, as
        this must be unique within a model.  Unfortunately, when we do
        partial evaluation, the full name changes!
        
        TODO: do we need a hash function?
        """
        return hash(self.fullname(cellml=True))

    def fullname(self, cellml=False, debug=False):
        """
        Return the full name of this variable, i.e.
        '(component_name,variable_name)'.

        If cellml is given as True, return the name in a form compatible with
        the CellML spec instead, i.e. component_name__variable_name, unless
        the component has its ignore_component_name property set, in which case
        just use variable_name.
        
        If debug is True, also give information about the variable object.
        """
        if hasattr(self, 'xml_parent'):
            parent_name = self.xml_parent.name
            ignore_component = self.component.ignore_component_name
        else:
            parent_name = '*orphan*'
            ignore_component = True
        if cellml:
            if ignore_component:
                vn = self.name
            else:
                vn = parent_name + u'__' + self.name
        else:
            vn = u'(' + parent_name + u',' + self.name + u')'
        if debug:
            vn = '%s@0x%x' % (vn, id(self))
        return vn

    @staticmethod
    def split_name(varname):
        """Split a variable name as given by cellml_variable.fullname into constituent parts.
        
        Returns a tuple (component name, local variable name).  If the component name cannot
        be identified, it will be returned as the empty string.
        """
        try:
            if varname[0] == u'(':
                cname, vname = varname[1:-1].split(u',')
            else:
                cname, vname = varname.split(u'__', 1)
        except ValueError:
            cname, vname = u'', varname
        return cname, vname

    @staticmethod
    def get_variable_object(model, varname):
        """Return the variable object that has name varname.
        
        This method tries to handle various forms of fully qualified variable name, i.e. names which
        include the name of the component the variable occurs in, including those created by
        CellMLTranslator.code_name.
        """
        varname = unicode(varname)
        if varname[:4] == 'var_':
            varname = varname[4:]
        cname, vname = cellml_variable.split_name(varname)
        if len(model.component) == 1 and cname == model.component.name:
            var = model.component.get_variable_by_name(vname)
        else:
            try:
                var = model.get_variable_by_name(cname, vname)
            except KeyError, e:
                try:
                    if cname:
                        vname = cname + u'__' + vname
                    var = model.component.get_variable_by_name(vname)
                except KeyError:
                    raise e
        return var
    
    def __str__(self):
        return 'cellml_variable' + self.fullname()
    
    def __repr__(self):
        return '<cellml_variable %s @ 0x%x>' % (self.fullname(), id(self))

    def get_component(self):
        return self.xml_parent
    component = property(get_component)

    @property
    def model(self):
        return self.component.xml_parent
    
    def get_units(self):
        """Get the cellml_units object giving this variable's units."""
        return self.component.get_units_by_name(self.units)

    def _add_dependency(self, dep):
        """Add a dependency of this variable.
        
        This could be an expression defining it, or a variable it's mapped from.
        Triggers a validation error if we already have another dependency,
        since a variable can't be defined in more than one way.
        """
        if self._cml_depends_on:
            if not dep in self._cml_depends_on:
                # Multiple dependencies. TODO: Give more info.
                raise MathsError(dep, u' '.join([
                    u'The variable',self.fullname(),
                    u'gets its value from multiple locations.']))
        else:
            self._cml_depends_on.append(dep)
        return
    
    def get_dependencies(self):
        """
        Return the list of things this variable depends on.
        """
        return self._cml_depends_on
    
    def _add_ode_dependency(self, independent_var, expr):
        """Add a dependency of this variable as the dependent variable in an ODE.
        
        independent_var is the corresponding independent variable, and expr is the
        expression defining the ODE.
        Triggers a validation error if the same ODE is defined by multiple expressions.
        """
        if independent_var in self._cml_depends_on_ode:
            if self._cml_depends_on_ode[independent_var] != expr:
                # Multiple definitions.  TODO: Give more info.
                raise MathsError(expr, u''.join([
                    u'There are multiple definitions of the ODE d',self.fullname(),
                    u'/d',independent_var.fullname()]))
        else:
            self._cml_depends_on_ode[independent_var] = expr
        return
    
    def _update_ode_dependency(self, free_var, defn):
        """Update an ODE dependency due to partial evaluation.
        
        When the PE processes the LHS of a derivative equation, it alters the independent
        variable to reference its source directly.  If this new source wasn't originally
        our independent variable's source (e.g. due to a units conversion expression) then
        this will break lookups of the ODE via _get_ode_dependency, so we need to update
        the reference.
        """
        if self.get_type() == VarTypes.Mapped:
            return self.get_source_variable()._update_ode_dependency(free_var, defn)
        self._cml_depends_on_ode.clear()
        self._cml_depends_on_ode[free_var] = defn
    
    def get_ode_dependency(self, independent_var, context=None):
        """
        Return the expression defining the ODE giving the derivative of this
        variable w.r.t. independent_var.
        Triggers a validation error if the ODE has not been defined.
        """
        if self.get_type() == VarTypes.Mapped:
            return self.get_source_variable().get_ode_dependency(independent_var, context=context)
        free_vars = self._cml_depends_on_ode.keys()
        free_vars = dict(zip(map(lambda v: v.get_source_variable(recurse=True), free_vars),
                             free_vars))
        independent_src = independent_var.get_source_variable(recurse=True)
        if not independent_src in free_vars:
            raise MathsError(context or self, u''.join([
                u'The ODE d',self.fullname(),u'/d',independent_var.fullname(),
                u'is used but not defined.']))
        return self._cml_depends_on_ode[free_vars[independent_src]]
    
    def get_all_expr_dependencies(self):
        """Return all expressions this variable depends on, either directly or as an ODE."""
        deps = filter(lambda d: isinstance(d, mathml_apply), self._cml_depends_on)
        deps.extend(self._cml_depends_on_ode.values())
        return deps

    def _set_source_variable(self, src_var):
        """Set this to be a mapped variable which imports its value from src_var.
        A validation error is generated if we are already mapped.
        """
        if not self._cml_source_var is None:
            # Mapping already exists
            model = self.model
            debug = model.get_option('debug')
            model.validation_error(u' '.join([
                'A variable with interface "in" may only obtain its',
                'value from one location.\nBoth',
                self._cml_source_var.fullname(debug=debug), 'and',
                src_var.fullname(debug=debug), 'are exported to', self.fullname(debug=debug)]))
        else:
            self._cml_source_var = src_var
            self._set_type(VarTypes.Mapped)
            self._add_dependency(src_var)
        return

    def get_source_variable(self, recurse=False):
        """
        Assuming this variable imports its value, return the variable
        from which we obtain our value.
        If our value is determined within this component, raise a
        TypeError.

        If recurse is set to True, recursively follow mappings until
        a non-mapped variable is found.
        """
        if self._cml_source_var is None:
            if recurse:
                src = self
            else:
                raise TypeError(u' '.join([
                    'Variable', self.fullname(), u'is not a mapped variable',
                    '(i.e. does not obtain its value from another component).'
                    ]))
        else:
            src = self._cml_source_var
            if recurse:
                src = src.get_source_variable(recurse=True)
        return src
    
    def _set_type(self, var_type, _orig=None):
        """Update the type of this variable.
        
        The caller should check that the update makes sense.
        
        If this variable already has type Mapped, then update the type of
        our source variable, instead.
        """
        # If this is becoming a state variable, increment its usage count
        if var_type is VarTypes.State and not self._cml_var_type is VarTypes.State:
            self._cml_usage_count += 1
        if self._cml_var_type == VarTypes.Mapped and not _orig is self:
            # Guard against infinite loops, since we haven't done a cyclic dependency check yet
            if _orig is None: _orig = self
            self.get_source_variable()._set_type(var_type, _orig=_orig)
        else:
            self._cml_var_type = var_type
        return

    def get_type(self, follow_maps=False):
        """Return the type of this variable.
        
        If follow_maps is True and the value of this variable is imported
        from another component, then return the type of the variable we
        get our value from instead.
        """
        if follow_maps and self._cml_var_type == VarTypes.Mapped:
            return self.get_source_variable().get_type(follow_maps=True)
        return self._cml_var_type

    def _used(self):
        """Note this variable as being used in an expression.
        
        Keep track of the usage count.
        If this is a mapped variable, note its source as being used,
        as well.

        Note that if a variable is used in 2 branches of a conditional
        then this counts as 2 uses.
        """
        self._cml_usage_count += 1
        if self._cml_var_type == VarTypes.MaybeConstant:
            self._cml_var_type = VarTypes.Constant
        elif self._cml_var_type == VarTypes.Mapped:
            self.get_source_variable()._used()
        return

    def get_usage_count(self):
        """
        Return the number of times this variable is used in an expression.
        """
        return self._cml_usage_count

    def _decrement_usage_count(self, follow_maps=True):
        """Decrement our usage count."""
        DEBUG('partial-evaluator', "Dec usage for", self.fullname())
        assert self._cml_usage_count > 0
        self._cml_usage_count -= 1
        if follow_maps and self._cml_var_type == VarTypes.Mapped:
            self.get_source_variable()._decrement_usage_count()
        # Note in the model if a usage count has decreased to 1, in
        # order to repeat the partial evaluation loop.
        if self._cml_usage_count == 1:
            model = self.xml_parent.xml_parent
            model._pe_repeat = u'yes'
        return

    def add_rdf_annotation(self, property, target, allow_dup=False):
        """Add an RDF annotation about this variable.
        
        property must be a tuple (qname, namespace_uri).
        target may either be a tuple as above, or a unicode string, in which
        case it is interpreted as a literal RDF node.
        
        If the variable does not already have a cmeta:id, one will be created
        for it with value self.fullname(cellml=True).

        The actual RDF will be added to the main RDF block in the model, which
        will be created if it does not exist.  Any existing annotations of this
        variable with the same property will be removed, unless allow_dup=True.
        """
        meta_id = self.cmeta_id
        if not meta_id:
            # Create ID for this variable, so we can refer to it in RDF
            meta_id = cellml_metadata.create_unique_id(self.model, unicode(self.fullname(cellml=True)))
            self.xml_set_attribute((u'cmeta:id', NSS['cmeta']), meta_id)
        property = cellml_metadata.create_rdf_node(property)
        target = cellml_metadata.create_rdf_node(target)
        source = cellml_metadata.create_rdf_node(fragment_id=meta_id)
        if allow_dup:
            cellml_metadata.add_statement(self.model, source, property, target)
        else:
            cellml_metadata.replace_statement(self.model, source, property, target)

    def get_rdf_annotation(self, property):
        """Get an RDF annotation about this variable.
        
        property must be a tuple (qname, namespace_uri).
        
        Will return the unique annotation found with source being this variable's id,
        and the given property.  If no annotation is found (or if the variable does
        not have a cmeta:id), returns None.  Throws an error if more than one
        annotation is found.
        """
        meta_id = self.cmeta_id
        if not meta_id:
            return None
        property = cellml_metadata.create_rdf_node(property)
        source = cellml_metadata.create_rdf_node(fragment_id=meta_id)
        return cellml_metadata.get_target(self.model, source, property)
    
    def get_rdf_annotations(self, property):
        """Get all RDF annotations about this variable that use the given property.
        
        property must be a tuple (qname, namespace_uri).
        
        Will return all annotations found with source being this variable's id,
        and the given property.  If no annotation is found (or if the variable does
        not have a cmeta:id), returns the empty list
        """
        meta_id = self.cmeta_id
        if not meta_id:
            return []
        property = cellml_metadata.create_rdf_node(property)
        source = cellml_metadata.create_rdf_node(fragment_id=meta_id)
        return cellml_metadata.get_targets(self.model, source, property)

    def remove_rdf_annotations(self, property=None):
        """Remove all RDF annotations about this variable.
        
        If property is given, only remove annotations with the given property.
        """
        meta_id = self.cmeta_id
        if meta_id:
            DEBUG('cellml-metadata', "Removing RDF annotations for", self, "with id", meta_id)
            source = cellml_metadata.create_rdf_node(fragment_id=meta_id)
            if property:
                property = cellml_metadata.create_rdf_node(property)
            cellml_metadata.remove_statements(self.model, source, property, None)

    def set_rdf_annotation_from_boolean(self, property, is_yes):
        """Set an RDF annotation as 'yes' or 'no' depending on a boolean value."""
        if is_yes:
            val = 'yes'
        else:
            val = 'no'
        self.add_rdf_annotation(property, val)

    def _set_binding_time(self, bt, temporary=False):
        """Set the binding time of this variable.
        
        Options are members of the BINDING_TIMES Enum.
        
        If temporary is True, then we're temporarily overriding the normal setting,
        so save any existing annotation for later replacement.
        """
        #print "set bt", self, bt, temporary
        assert bt in BINDING_TIMES
        if temporary and not hasattr(self, '_cml_saved_bt'):
            self._cml_saved_bt = self.get_rdf_annotation(('pe:binding_time', NSS[u'pe']))
        self.add_rdf_annotation(('pe:binding_time', NSS[u'pe']), str(bt))
        self._cml_binding_time = bt
        return
    
    def _unset_binding_time(self, only_temporary=False):
        """Unset any stored binding time.
        
        If the stored BT was a temporary setting, replace it with the original value.
        """
        #print "unset bt", self, only_temporary
        self._cml_binding_time = None
        if hasattr(self, '_cml_saved_bt'):
            if self._cml_saved_bt:
                self._set_binding_time(getattr(BINDING_TIMES, self._cml_saved_bt))
            else:
                self.remove_rdf_annotations(('pe:binding_time', NSS[u'pe']))
            del self._cml_saved_bt
        elif not only_temporary:
            self.remove_rdf_annotations(('pe:binding_time', NSS[u'pe']))
        
    def _get_binding_time(self):
        """Return the binding time of this variable, as a member of
        the BINDING_TIMES Enum.

        This method will (try to) compute & cache the binding time if
        necessary.

        Variables are classified based on their type:
          State, Free -> dynamic
          Constant -> static
          Mapped -> binding time of source variable
          Computed -> binding time of defining expression
        """
        if self._cml_binding_time is None:
            # Check for an annotation setting the BT
            bt_annotation = self.get_rdf_annotation(('pe:binding_time', NSS[u'pe']))
            if bt_annotation:
                bt = getattr(BINDING_TIMES, bt_annotation)
                DEBUG('partial-evaluator', "BT var", self.fullname(), "is annotated as", bt)
            elif self.pe_keep:
                # This variable must appear in the specialised model
                bt = BINDING_TIMES.dynamic
                DEBUG('partial-evaluator', "BT var", self.fullname(), "is kept")
            else:
                # Compute BT based on our type
                t = self.get_type()
                DEBUG('partial-evaluator', "BT var", self.fullname(), "type", str(t))
                if t in [VarTypes.State, VarTypes.Free, VarTypes.Unknown]:
                    bt = BINDING_TIMES.dynamic
                elif t == VarTypes.Constant:
                    bt = BINDING_TIMES.static
                elif t == VarTypes.Mapped:
                    bt = self.get_source_variable()._get_binding_time()
                elif t == VarTypes.Computed:
                    bt = self._cml_depends_on[0]._get_binding_time()
                else:
                    raise TypeError("Unexpected variable type " + str(t) +
                                    " of variable " + self.fullname() +
                                    " in BTA.")
                DEBUG('partial-evaluator', "BT var", self.fullname(), "is", bt)
            self._set_binding_time(bt)
        else:
            bt = self._cml_binding_time
        return bt
    
    def is_statically_const(self, ignore_annotations=False):
        """Determine loosely if this variable is considered constant.
        
        Checks if we're Constant, or Computed with a static binding time (or
        of unknown type).
        
        If ignore_annotations is True, will ignore cached binding time values and
        pe:keep annotations.  It instead finds all variables we depend on, directly or
        indirectly, and gives a dynamic result iff any is a state or free variable.
        """
        result = False
        t = self.get_type()
        if t in [VarTypes.Constant, VarTypes.Unknown]:
            result = True
        elif t == VarTypes.Computed:
            if ignore_annotations:
                dependencies = self.model.calculate_extended_dependencies([self])
                result = True
                for node in dependencies:
                    if isinstance(node, cellml_variable) and node.get_type() in [VarTypes.State, VarTypes.Free]:
                        result = False
                        break
            else:
                result = self._get_binding_time() == BINDING_TIMES.static
        elif t == VarTypes.Mapped:
            result = self.get_source_variable().is_statically_const(ignore_annotations)
        return result

    def set_value(self, value, ode=None, follow_maps=True):
        """Set the value of this variable.

        Expects a floating point or boolean value.

        If ode is given, it should be an instance of cellml_variable.
        In this case, we're setting the value of d(self)/d(ode).

        If this is a mapped variable, assign the value to its source
        variable instead, unless follow_maps is set to False
        """
        #print "set_value", self, ode, value
        if follow_maps and self.get_type() == VarTypes.Mapped:
            self.get_source_variable().set_value(value, ode=ode)
        else:
            assert type(value) in [types.FloatType, types.BooleanType]
            self._cml_value[ode] = float(value)
        return
    def unset_values(self):
        """Unset all values for this variable set with set_value."""
        #print "unset_values", self
        if self.get_type() == VarTypes.Mapped:
            self.get_source_variable().unset_values()
        self._cml_value.clear()
    def get_value(self, ode=None):
        """Return the value of this variable.

        If a value has been set with set_value(), return that.
        Otherwise, the behaviour depends on the type of this variable.
        If it is Mapped, return the value of the source variable.
        If it is Constant or State, return its initial value.
        Otherwise, raise an EvaluationError.

        If ode is given, it should be an instance of cellml_variable.
        In this case, we're getting the value of d(self)/d(ode).
        """
        # TODO: Might want to alter the handling of initial value for an ODE?
        if ode in self._cml_value:
            val = self._cml_value[ode]
        elif self.get_type() == VarTypes.Mapped:
            val = self.get_source_variable().get_value(ode=ode)
        elif ode is None and self.get_type() in [VarTypes.Unknown,
                                                 VarTypes.State, VarTypes.Constant]:
            if hasattr(self, 'initial_value'):
                val = float(self.initial_value)
            else:
                raise EvaluationError("Variable " + self.fullname() + " has no initial value set.")
        elif self.get_type() == VarTypes.Computed and self._get_binding_time() == BINDING_TIMES.static:
            # Evaluate the defining expression
            val = self._cml_depends_on[0].evaluate()
        else:
            raise EvaluationError("Unable to find a suitable value for variable " + self.fullname())
        return val

    @property
    def is_modifiable_parameter(self):
        """Whether this variable is a parameter that should be modifiable at run-time."""
        return (self.get_type() == VarTypes.Constant and
                self.get_rdf_annotation(('pycml:modifiable-parameter', NSS['pycml'])) == 'yes')
    def set_is_modifiable_parameter(self, is_param):
        """Set method for the is_modifiable_parameter property.
        
        We need a separate method for this to bypass Amara's property setting checks.
        """
        if is_param and self.get_type() != VarTypes.Constant:
            raise ValueError("A non-constant variable (%s) cannot be set as a parameter" % (self.fullname(),))
        self.set_rdf_annotation_from_boolean(('pycml:modifiable-parameter', NSS[u'pycml']), is_param)

    @property
    def is_derived_quantity(self):
        """Whether this variable should be included in reports of derived quantities."""
        return self.get_rdf_annotation(('pycml:derived-quantity', NSS['pycml'])) == 'yes'
    def set_is_derived_quantity(self, is_dq):
        """Set method for the is_derived_quantity property.
        
        We need a separate method for this to bypass Amara's property setting checks.
        """
        self.set_rdf_annotation_from_boolean(('pycml:derived-quantity', NSS[u'pycml']), is_dq)

    @property
    def is_output_variable(self):
        """Whether a protocol has requested this variable as a model output."""
        return self.get_rdf_annotation(('pycml:output-variable', NSS['pycml'])) == 'yes'
    def set_is_output_variable(self, is_ov):
        """Set method for the is_output_variable property.
        
        We need a separate method for this to bypass Amara's property setting checks.
        """
        self.set_rdf_annotation_from_boolean(('pycml:output-variable', NSS[u'pycml']), is_ov)

    @property
    def pe_keep(self):
        """Whether PE should retain this variable in the specialised model."""
        return (self.get_rdf_annotation(('pe:keep', NSS[u'pe'])) == 'yes' or
                self.is_modifiable_parameter or
                self.is_derived_quantity or
                self.is_output_variable)
    def set_pe_keep(self, keep):
        """Set method for the pe_keep property.
        
        We need a separate method for this to bypass Amara's property setting checks.
        """
        self.set_rdf_annotation_from_boolean(('pe:keep', NSS[u'pe']), keep)

    @property
    def oxmeta_name(self):
        """The canonical name of this variable, as given by Oxford metadata.

        Returns the empty string if no annotation is given.
        """
        for annotation in self.get_rdf_annotations(('bqbiol:is', NSS['bqbiol'])):
            name = cellml_metadata.namespace_member(annotation, NSS['oxmeta'], wrong_ns_ok=True)
            if name:
                return name
        # No suitable annotation found
        return ""
    def set_oxmeta_name(self, name):
        """Set method for the oxmeta_name property.
        
        Sets a bqbiol:is RDF annotation with the name.
        
        We need a separate method for this to bypass Amara's property setting checks.
        """
        self.add_rdf_annotation(('bqbiol:is', NSS['bqbiol']), ('oxmeta:'+name, NSS['oxmeta']))

    def _reduce(self, update_usage=False):
        """Reduce this dynamic variable that is being kept.

        If it has a static source, then turn it into a constant.
        Otherwise make it depend directly on an appropriate source:
        either one of its source variables that is also being kept,
        or the ultimate defining expression.

        If update_usage is True then this variable is used in an equation,
        so reduce the usage of any source variables it no longer depends on.
        """
        assert self.pe_keep
        src = self.get_source_variable(recurse=True)
        if src._get_binding_time() is BINDING_TIMES.static:
            # Become a constant
            self._cml_var_type = VarTypes.Constant
            # Set our value
            value = unicode(src.get_value())
            self.initial_value = value
            # Fix up usage counts
            if update_usage:
                self.get_source_variable()._decrement_usage_count()
            # We now have no dependencies
            self._cml_depends_on = []
            self._cml_source_var = None
        elif self.get_type() == VarTypes.Mapped:
            # Manually recurse down maps to find a suitable source
            src = self.get_source_variable()
            while not src.pe_keep and src.get_type() == VarTypes.Mapped and src.get_usage_count() == 1:
                src = src.get_source_variable()
            if src.pe_keep or src.get_usage_count() > 1:
                # Depend directly on src
                self._cml_depends_on = [src]
                self._cml_source_var = src
                # Fix up usage counts
                if update_usage:
                    src._used()
                    self.get_source_variable()._decrement_usage_count()
            else: # src.get_type() != VarTypes.Mapped
                # This variable is the only reference to the ultimate defining
                # expression, so become computed.
                self._cml_var_type = VarTypes.Computed
                defn = src.get_dependencies()[0]
                assert isinstance(defn, mathml_apply)
                ## Move the definition to this component
                #defn._unset_cached_links()
                #defn.xml_parent.xml_remove_child(defn)
                #self.component.math.xml_append(defn)
                # Schedule the LHS of the defining expression for update
                defn._cml_assigns_to = self
                defn._pe_process = u'retarget'
                # Fix up usage counts
                if update_usage:
                    self.get_source_variable()._decrement_usage_count()
                # Depend directly on the assignment expression
                self._cml_depends_on = [defn]
                self._cml_source_var = None
    
    @staticmethod
    def create_new(elt, name, units, id=None, initial_value=None, interfaces={}):
        """Create a new <variable> element with the given name and units.
        
        Optionally id, initial_value, and interfaces may also be given.
        
        elt may be any existing XML element.
        """
        attrs = {(u'units', None): unicode(units),
                 (u'name', None): unicode(name)}
        if id:
            attrs[(u'cmeta:id', NSS[u'cmeta'])] = unicode(id)
        if initial_value is not None and initial_value != u'':
            attrs[(u'initial_value', None)] = unicode(initial_value)
        for iface, val in interfaces.items():
            attrs[(iface + u'_interface', None)] = unicode(val)
        new_elt = elt.xml_create_element(u'variable', NSS[u'cml'], attributes=attrs)
        return new_elt

class UnitsSet(set):
    """A set of cellml_units objects.

    This class behaves like a normal set, but also has additional
    methods for operations specific to sets of <units> objects:
      simplify - allow for multiplication of sets of units
      dimensionally_equivalent - compare 2 sets of units for dimensional equivalence
      description - describe the units in this set

    All units in the set (normally) must be dimensionally equivalent.  The exception is when dealing with
    Functional Curation protocols which can defined units conversion rules for non-scaling cases.  We can
    then have sets of alternative units in different dimensions, and the get_consistent_set method helps
    with selecting from these.
    """
    def __new__(cls, iterable=[], expression=None):
        """
        Work around annoyance in set implementation of python 2.4.2c1 and on.
        setobject.c checks for keyword arguments in its __new__ instead of its
        __init__, so we get an error
        'TypeError: set() does not take keyword arguments' if we don't do this.
        """
        return super(UnitsSet, cls).__new__(cls, iterable)
    
    def __init__(self, iterable=[], expression=None):
        super(UnitsSet, self).__init__(iterable)
        self._expression = expression
        self._sources = {}
        return

    def copy(self):
        """Do a shallow copy of this UnitsSet."""
        new_set = super(UnitsSet, self).copy()
        new_set._expression = self._expression
        new_set._sources = {}
        for units, src_list in self._sources.iteritems():
            new_set._sources[units] = copy.copy(src_list)
        return new_set
    
    def get_consistent_set(self, desired_units):
        """Extract a subset of the units in this set that are dimensionally equivalent.
        
        When dealing with potential non-scaling conversions, an expression may have potential units that are
        not within the same dimension.  However, when deciding on the units for the expression most operations
        need to handle a set of options that *are* within the same dimension, and this operation supports that.

        Given a cellml_units object for the desired units of the expression, this method first tries to create
        a UnitsSet containing just our members in the same dimension.  If we have no members in the desired
        dimension, then we check that all members of this set are in the same dimension, and return the set
        itself - there should never be more than 2 different dimensions to choose from (I hope!).
        """
        new_set = UnitsSet([], expression=self._expression)
        for units in self:
            if units.dimensionally_equivalent(desired_units):
                new_set.add(units)
                new_set._sources[units] = copy.copy(self._sources[units])
        if not new_set:
            rep_u = self.extract()
            for units in self:
                if not units.dimensionally_equivalent(rep_u):
                    raise ValueError("Unexpected dimensional variation in UnitsSet; " + rep_u.description() + " and " + units.description()
                                     + " do not match.")
            new_set = self
        return new_set

    def equals(self, other):
        """Test whether the units in the set are equal to those in another set."""
        try:
            equal = self.extract(check_equality=True).equals(other.extract(check_equality=True))
        except ValueError:
            equal = False
        return equal

    def extract(self, check_equality=False):
        """Extract a representative element from this set.

        This is intended to be used to get the cellml_units object from a singleton set.
        
        If check_equality is True, check that all members of this set have the same multiplicative factor.
        """
        representative = iter(self).next()
        if check_equality:
            for u in self:
                if not u._rel_error_ok(u.expand().get_multiplicative_factor(),
                                       representative.expand().get_multiplicative_factor(),
                                       1e-6):
                    raise ValueError("UnitsSet equality check failed")
                if u.is_simple() and not u._rel_error_ok(u.expand().get_offset(),
                                                         representative.expand().get_offset(),
                                                         1e-6):
                    raise ValueError("UnitsSet equality check failed")
        return representative
    
    def set_expression(self, expr):
        """Store a reference to the expression that has these units."""
        self._expression = expr
        return

    def _add_source(self, units, src_units_set, src_units):
        """Add source units for the given units.

        This stores references to how units were arrived at when doing a
        simplify.  It manages lists of (src_units_set, src_units) pairs for
        each units definition in this set.

        If src_units_set has no associated expression, then it is
        considered to be a temporary object, created for example whilst
        doing an n-ary times operation.  In this case, we add its list of
        sources for src_units to our sources list instead.
        """
        if not units in self._sources:
            self._sources[units] = []
        if not hasattr(src_units_set, '_expression'):
            print self.description(), units.description()
            print src_units_set.description(), src_units.description()
        if src_units_set._expression:
            self._sources[units].append((src_units_set, src_units))
        else:
            try:
                self._sources[units].extend(src_units_set._sources[src_units])
            except KeyError:
                # No sources list found.  This can happen if we do
                # dimensionless.simplify(...), for example, to evaluate powers.
                pass
        return
    
    def _get_sources(self, units):
        """Return the sources list for the given units."""
        return self._sources.get(units, [])

    def get_expression(self):
        """Return an expression that has these units."""
        return self._expression
    
    def simplify(self, other_units=None, other_exponent=1):
        """Simplify the units in this set.

        Each cellml_units object in this set is simplified, and a new
        UnitsSet returned with the results.

        If other_units is not None, then products of units are
        calculated.  The returned set is essentially the cartesian
        product of the input sets under the simplify operator, i.e.
        u1.simplify(other_units=u2, other_exponent=other_exponent)
        will be called for each member u1 of self and each member u2
        of other_units (if other_units is a UnitsSet; otherwise
        u2=other_units).
        """
        result_set = UnitsSet()
        for units in self:
            if other_units is None:
                res_u = units.simplify()
                result_set.add(res_u)
                result_set._add_source(res_u, self, units)
            else:
                if isinstance(other_units, cellml_units):
                    other_units = UnitsSet([other_units])
                for u in other_units:
                    res_u = units.simplify(u, other_exponent)
                    result_set.add(res_u)
                    result_set._add_source(res_u, self, units)
                    result_set._add_source(res_u, other_units, u)
        return result_set

    def dimensionally_equivalent(self, other_units):
        """Check for dimensional equivalence between sets of units.

        Since all units in each set should be equivalent, we just compare
        an arbitrary member from each set.

        other_units may be a single cellml_units instance, in which case we
        compare an arbitrary member of self to it.
        """
        u1 = self.extract()
        u2 = other_units.extract()
        return u1.dimensionally_equivalent(u2)

    def description(self):
        """Describe these units.

        Shows the descriptions of each member, as though this were a set of
        unicode strings.  If multiple members have the same description,
        only one instance is shown.  If only one description is being shown,
        then the curly brackets are not added.
        """
        desc = list(set(u.description() for u in self))
        desc.sort()
        if len(desc) > 1:
            desc = u'{' + u','.join(desc) + u'}'
        else:
            desc = desc[0]
        return desc

class cellml_units(Colourable, element_base):
    """
    Specialised units class.
    Contains useful methods for defining the standard units dictionary,
    and checking for units consistency.

    After being defined, a units definition should be regarded as being
    immutable, as should individual <unit> elements, so the expansion
    and simplification methods create new objects, which don't really
    live in the document tree (they're not a child of any element in the
    model), although they do have their xml_parent attribute set.

    Note that <unit> elements must not be made a child of more than one
    <units> element, otherwise the linked lists get tangled up.
    """
    
    def __init__(self):
        super(cellml_units, self).__init__()
        self._cml_expanded = None
        self._cml_simplified = {}
        self._cml_generated = False
        self._cml_quotients = {}
        self._cml_hash = None
        return

    def __repr__(self):
        return '<cellml_units %s @ 0x%x>' % (self.name, id(self))

    def __contains__(self, item):
        """Prevent the default implementation using Amara's __iter__ to give unexpected results."""
        return False

    @property
    def _hash_tuple(self):
        """Generate a tuple used as the basis for our hash value."""
        if self._cml_hash is None:
            hash_list = [self.is_base_unit()]
            if not self._cml_generated:
                hash_list.append(self.name)
            hash_list.extend(list(getattr(self, u'unit', [])))
            hash_tup = tuple(hash_list)
            self._cml_hash_tup = hash_tup
            self._cml_hash = hash(hash_tup)
        return self._cml_hash_tup

    @property
    def uniquify_tuple(self):
        """For avoiding duplicate identical units definitions.

        Based on description(cellml=True), since that is what we really want to be unique.
        Also includes offset information, since that is omitted in the name given by description.
        """
        l = [self.description(cellml=True)]
        if self.is_simple():
            # Include offset information
            l.append((self.unit.get_units_element().uniquify_tuple,
                      self.unit.get_offset()))
        return tuple(l)

    def __hash__(self):
        """Generate a hash for these units.

        Hashes a tuple, the first element of which is the result of self.is_base_unit(),
        the second of which is our name if we are not auto-generated,
        and the remaining elements of which are our <unit> elements.
        """
        if self._cml_hash is None:
            _ = self._hash_tuple
        return self._cml_hash

    def __cmp__(self, other):
        """Compare 2 units objects, by comparing their hashes.

        Means using units as dictionary keys is more sane."""
        if isinstance(other, cellml_units):
            if hash(self) == hash(other):
                return cmp(self._cml_hash_tup, other._cml_hash_tup)
            else:
                return cmp(hash(self), hash(other))
        else:
            return super(cellml_units, self).__cmp__(other)

# The following currently causes infinite loops/recursions
##    def __eq__(self, other):
##        """Compare 2 units elements for equality.
##        return self.equals(other)
##    def __ne__(self, other):
##        """Inverse of self.__eq__(other)."""
##        return not self.__eq__(other)

    def _rel_error_ok(self, value1, value2, tol):
        """Test if the relative difference of 2 values is within tol."""
        if abs(value1) == 0.0:
            return abs(value2) < tol
        else:
            return (abs(value1-value2)/abs(value1)) < tol

    def equals(self, other):
        """Compare 2 units elements for equality.

        Two units are deemed equal if they are both dimensionally equivalent and have the same
        multiplicative factor (to within a relative tolerance of 10^-6).
        
        If both are simple units, they must also have the same offset.
        """
        equal = isinstance(other, cellml_units) and \
                self.dimensionally_equivalent(other) and \
                self._rel_error_ok(self.expand().get_multiplicative_factor(),
                                   other.expand().get_multiplicative_factor(),
                                   1e-6)
        if equal and self.is_simple():
            equal = self._rel_error_ok(self.get_offset(),
                                       other.get_offset(),
                                       1e-6)
        return equal

    @property
    def model(self):
        return self.rootNode.model

    def extract(self, check_equality=False):
        """Return these units.

        Used for interface compatibility with UnitsSet."""
        return self

    def copy(self):
        """Return a new UnitsSet containing this cellml_units object.

        Used for interface compatibility with UnitsSet, where the method performs a shallow copy.
        """
        return UnitsSet([self])

    def description(self, force=False, cellml=False):
        """Return a human-readable name for these units.
        
        The name will be descriptive and based on the consituent <unit> elements, e.g. 'volt per second^2'

        By default, if these are user-defined units, then return self.name.  Set force to True to override this behaviour.

        Set cellml to True to get a description that is also a valid CellML identifier.
        """
        if self.is_base_unit():
            desc = self.name
        elif not force and not self._cml_generated:
            desc = self.name
        else:
            descs, per_descs = [], []
            # Multiplier
            m = self.get_multiplier()
            if m != 1:
                descs.append(unicode(m))
            # Constituent units
            dimensionless = self.get_units_by_name('dimensionless')
            for unit in self.unit:
                if unit.get_units_element() is dimensionless:
                    continue
                desc = [getattr(unit, u'prefix_', u''), unit.get_units_element().name]
                e = unit.get_exponent()
                if int(e) == e:
                    # Cast to integer so name looks nicer.
                    e = int(e)
                if abs(e) != 1:
                    desc.extend(['^', str(abs(e))])
                desc = u''.join(desc)
                if e < 0:
                    per_descs.append(desc)
                else:
                    descs.append(desc)
            # Sort unit descriptions for readability
            descs.sort()
            descs = u' '.join(descs)
            per_descs.sort()
            desc = u' per '.join([descs] + per_descs)
            if not desc:
                desc = u'dimensionless'
            elif descs:
                desc = u"'" + desc + u"'"
            else:
                # Remove unwanted space from the start
                desc = u"'" + desc[1:] + u"'"
            # Convert to CellML identifier?
            if cellml:
                desc = desc.replace(u"'", u"").replace(u"^", u"")
                desc = desc.replace(u"*", u"").replace(u".", u"_")
                desc = desc.replace(u" ", u"_")
        return desc

    def get_units_by_name(self, uname):
        """Return an object representing the element that defines the units named `uname'."""
        # Look up units in our parent model or component element instead
        return self.xml_parent.get_units_by_name(uname)

    def expand(self):
        """Expand to a product of powers of base units.
        
        Expand this units definition according to the algorithm given in appendix C.3.4 of the CellML spec.
        
        Caches and returns the resulting <units> object, which will be
        newly created if there are any changes made, or this object if not.
        """
        if self._cml_expanded is None:
            # Do the expansion
            if self.is_base_unit():
                # We are a base unit, so stop the recursion; the result is us
                self._cml_expanded = self
            elif self.is_simple():
                # Simple units definition, i.e. only one <unit> element,
                # exponent=1, and referenced units are simple or base.
                # Expand the units referenced by our <unit> element
                expanded_child = self.unit.get_units_element().expand()
                # Get our multiplicative factor & offset as numbers
                m1p1 = self.unit.get_multiplicative_factor()
                o1 = self.unit.get_offset()
                if expanded_child.is_base_unit():
                    # New multiplier, etc. are just ours
                    m_new = m1p1
                    o_new = o1
                    # Referenced units are expanded_child
                    ref_obj_new = expanded_child
                else:
                    # Get the multiplier & offset of our expanded child
                    m2p2 = expanded_child.unit.get_multiplicative_factor()
                    o2 = expanded_child.unit.get_offset()
                    # Calculate new multiplier, etc. per Equation (11)
                    m_new = m1p1*m2p2
                    o_new = o1 + o2/m1p1
                    # Referenced units are those referenced by expanded_child
                    # These will be base units
                    ref_obj_new = expanded_child.unit.get_units_element()
                # Create the new units & unit elements
                attrs = {u'name': self.name}
                self._cml_expanded = self.xml_create_element(u'units',
                                                             NSS[u'cml'],
                                                             attributes=attrs)
                self._cml_expanded._cml_generated = True
                attrs = {u'units': ref_obj_new.name,
                         u'multiplier': unicode(m_new),
                         u'offset': unicode(o_new)}
                unit = self.xml_create_element(u'unit', NSS[u'cml'],
                                               attributes=attrs)
                # Manually set the reference object for unit, since we have
                # it handy
                unit._set_units_element(ref_obj_new)
                self._cml_expanded.xml_append(unit)
            else:
                # Complex units definition, i.e. multiple <unit> elements,
                # or non-unitary exponent, at some point within this defn.
                # Create the new units element so we can add children to it
                attrs = {u'name': self.name}
                exp_units = self.xml_create_element(u'units', NSS[u'cml'],
                                                    attributes=attrs)
                exp_units._cml_generated = True
                # Compute m_t (Equation (18)) expanding referenced
                # units as we go
                m_t = 1
                for unit in self.unit:
                    m_t *= unit.get_multiplicative_factor() # * is assoc. :)
                    # We'll need the exponent to modify units referenced
                    # by this reference
                    e = unit.get_exponent()
                    # Expand referenced units
                    exp_u = unit.get_units_element().expand()
                    if exp_u.is_base_unit():
                        # Create reference to exp_u
                        attrs = {u'units': exp_u.name,
                                 u'exponent': unicode(e)}
                        u = self.xml_create_element(u'unit', NSS[u'cml'],
                                                    attributes=attrs)
                        exp_units.xml_append(u)
                    else:
                        # Process units referenced by exp_u, which will be
                        # base units.
                        for u in exp_u.unit:
                            m_t *= u.get_multiplicative_factor() ** e
                            attrs = {u'units': u.units,
                                     u'exponent': unicode(
                                u.get_exponent() * e)}
                            u_new = u.xml_create_element(u'unit',
                                                         NSS[u'cml'],
                                                         attributes=attrs)
                            exp_units.xml_append(u_new)
                # Set the multiplier for the expanded units to m_t.
                # Since a <units> elements doesn't have a multiplier
                # attribute, we set it on the first <unit> element.
                # Note that all the <unit> elements have been newly created,
                # so have an implicit multiplier of 1 currently.
                # Note also that the fact that each <unit> has an exponent
                # doesn't matter, since the exponent doesn't affect the
                # multiplier.
                # Alan pointed out a bug: if the <unit> elements are
                # in non-canonical order then we can get identical
                # units (according to the intended canonical form)
                # comparing as non-equal, because expansion put the
                # multiplicative factor on different <unit> elements
                # in each case.  Hence we have to put the factor on
                # the first <unit> element according to a canonical
                # sorting order.
                # TODO: Be a bit cleverer about this?  In the base unit case
                # above we could do exp_units.xml_append(unit.clone()) and
                # not have its multiplicative factor contributing to m_t.
                # Would then need to multiply here, since first unit may
                # already have multiplier != 1.
                # Alternative idea from Alan: I have just noticed that
                # when expanding a complex units definition based on a
                # complex units definition, we can get away with not
                # sorting the unit elements. All that is required is
                # to create a new unit element which type is
                # dimensionless and value of its multiplier attribute
                # m*10^(p*e). This avoids having to sort things and
                # propagating the multiplier attribute...
                first_unit = sorted(exp_units.unit, key=lambda u: u.units)[0]
                first_unit.multiplier = unicode(m_t)
                # Cache result
                self._cml_expanded = exp_units
            # Set parent of the expanded units to be our parent
            self._cml_expanded.xml_parent = self.xml_parent
        # Returned the cached result
        return self._cml_expanded

    def is_base_unit(self):
        """Return True iff this is a base units definition."""
        return getattr(self, u'base_units', u'no') == u'yes'
    def is_simple(self):
        """Return True iff this is a simple units definition.

        Units are simple if:
          there is 1 <unit> element
          the exponent is omitted or has value 1.0
          the referenced units are simple or base units
        """
        simple = False
        if not self.is_base_unit() and len(self.unit) == 1:
            if self.unit.get_exponent() == 1.0:
                u = self.unit.get_units_element()
                if u.is_base_unit() or u.is_simple():
                    simple = True
        return simple

    def get_multiplicative_factor(self):
        """Return the multiplicative factor of this units definition.

        The multiplicative factor of a units definition can be defined as
        the product of the multiplicative factors of its unit children.
        """
        m = reduce(operator.mul,
                   map(lambda unit: unit.get_multiplicative_factor(),
                       getattr(self, u'unit', [])),
                   1)
        return m

    def get_multiplier(self):
        """Return the multiplier of this units definition.

        The multiplier of a units definition can be defined as the product
        of the multipliers of its unit children.
        """
        return reduce(operator.mul,
                      map(lambda unit: unit.get_multiplier(),
                          getattr(self, u'unit', [])),
                      1)

    def get_offset(self):
        """Return the offset associated with this units definition.

        If these are simple units, return the offset on our unit child.
        Otherwise, return 0.
        """
        if self.is_simple():
            o = self.unit.get_offset()
        else:
            o = 0
        return o
    
    @staticmethod
    def create_new(parent, name, bases, add_to_parent=False, standard=False):
        """Create a new units definition element.

        It requires either a cellml_model or cellml_component element
        to be passed as the parent for the definition.  If
        add_to_parent is set to true the new units element will be
        appended to parent's children.

        The bases parameter should be a list of dictionaries suitable
        for use as the keyword arguments of cellml_units._based_on.
        If the list is empty it will be defined as a base unit.
        """
        # Create the units element
        attrs = {u'name': unicode(name)}
        if not bases:
            attrs[u'base_units'] = u'yes'
        if standard:
            attrs[u'standard'] = u'yes'
        u = parent.xml_create_element(u'units', NSS[u'cml'], attributes=attrs)
        if add_to_parent:
            parent.xml_append(u)
        else:
            u.xml_parent = parent # Hack so units lookups work
        # Add references to units we're based on
        for basis in bases:
            u._based_on(**basis)
        return u

    def _based_on(self, units=None, prefix=None, exponent=None,
                  multiplier=None, offset=None):
        """Convenience function for defining new units."""
        for v in ['units', 'prefix', 'exponent', 'multiplier', 'offset']:
            # Coerce from str to unicode
            exec "if type(%s) == str: %s = unicode(%s)" % (v,v,v)
            # Check type
            exec "assert(%s is None or type(%s) == unicode)" % (v,v)
        assert(not hasattr(self, 'base_units') or self.base_units == u'no')
        attrs = {u'units': units}
        if offset:
            # Simple units definition
            assert(not hasattr(self, 'unit'))
            attrs[u'offset'] = offset
            if not prefix is None: attrs[u'prefix'] = prefix
            if not exponent is None: attrs[u'exponent'] = exponent
        else:
            # Complex units definition
            if not prefix is None: attrs[u'prefix'] = prefix
            if not exponent is None: attrs[u'exponent'] = exponent
            if not multiplier is None: attrs[u'multiplier'] = multiplier
        self.xml_append(self.xml_create_element(u'unit', NSS[u'cml'],
                                                attributes=attrs))
        return

    # Class attribute, so all <units> elements share the same value.
    units_name_counter = [0]
    def simplify(self, other_units=None, other_exponent=1):
        """Simplify these units.
        
        Create a new <units> element representing a simplified version of
        this element.  This implements the algorithm of appendix C.3.1.  It
        is however slightly different, in order to allow for units
        conversions rather than just preserving dimensional equivalence.

        If other_units is not None, then produce a simplified version of
        the product of these units and other_units.  other_units are
        considered to be raised to the power of other_exponent (so a
        quotient can be performed by passing other_exponent=-1).  Note that
        other_exponent should have numeric type.

        If other_units is a UnitsSet, then we construct a UnitsSet
        containing self, and call the simplify method on that, thus
        returning a UnitsSet instance.

        Multiplicative factors on the original <unit> objects are
        maintained on the generated references, by taking the product of
        all factors from <unit> objects that contribute to the new <unit>
        object.  Note that this means that we may need to retain a
        reference to dimensionless with a multiplier, for example if the
        quotient of centimetres by metres is taken.

        Offsets are only permitted on simple units, so may not appear where
        there are multiple <unit> elements.  Hence when a product of units
        is taken, any offsets are discarded.
        """
        if isinstance(other_units, UnitsSet):
            # Use the set version of simplify
            return UnitsSet([self]).simplify(other_units, other_exponent)
        
        # Check for result in cache (see #1714)
        if (other_units, other_exponent) in self._cml_simplified:
            return self._cml_simplified[(other_units, other_exponent)]

        # Make a list of all <unit> elements to be included
        units = []
        if self.is_base_unit():
            # We are a base unit, so invent a reference to ourselves
            u = self.xml_create_element(u'unit', NSS[u'cml'],
                                        attributes={u'units': self.name})
            u.xml_parent = self # Hack, but it might need a parent...
            u._set_units_element(self)
            units.append(u)
        else:
            units = list(self.unit)
        our_unit_elements = frozenset(units)
        other_unit_elements = None
        if not other_units is None:
            if other_units.is_base_unit():
                # other_units are base units, so invent a reference
                attrs = {u'units': other_units.name,
                         u'exponent': unicode(other_exponent)}
                u = self.xml_create_element(u'unit', NSS[u'cml'],
                                            attributes=attrs)
                u.xml_parent = other_units
                u._set_units_element(other_units)
                units.append(u)
                if other_exponent == 1:
                    other_unit_elements = frozenset([u])
            else:
                if other_exponent == 1:
                    units.extend(list(other_units.unit))
                    other_unit_elements = frozenset(other_units.unit)
                else:
                    for unit in other_units.unit:
                        # Need to create a new <unit> element with different
                        # exponent and multiplier
                        u = unit.clone()
                        u.exponent = unicode(other_exponent *
                                             u.get_exponent())
                        u.multiplier = unicode(u.get_multiplier() **
                                               other_exponent)
                        units.append(u)
        # Sort <unit> elements according to the <units> objects they
        # reference
        dimensionless = self.get_units_by_name(u'dimensionless')
        d = {dimensionless: []}
        for unit in units:
            obj = unit.get_units_element()
            if obj in d:
                d[obj].append(unit)
            else:
                d[obj] = [unit]
        # Collapse equivalent units references into a single reference.
        # That is, all references to the same object get collapsed into a new
        # reference to that object with different exponent, etc.
        for obj in d.keys():
            if obj != dimensionless and len(d[obj]) > 1:
                # Sum the exponents
                expt = sum(map(lambda u: u.get_exponent(), d[obj]))
                # If exponents cancel, replace with ref to dimensionless
                if expt == 0:
                    attrs = {u'units': u'dimensionless'}
                else:
                    attrs = {u'units': d[obj][0].units,
                             u'exponent': unicode(expt)}
                # Compute the multiplier for the new unit reference, as
                # the product of multiplicative factors on the originals
                m = reduce(operator.mul,
                           map(lambda u: u.get_multiplicative_factor(),
                               d[obj]))
                attrs[u'multiplier'] = unicode(m)
                # Create a new reference
                new = self.xml_create_element(u'unit', NSS[u'cml'],
                                              attributes=attrs)
                new.xml_parent = self
                if expt == 0:
                    # Note an extra reference to dimensionless...
                    d[dimensionless].append(new)
                    # ...and remove the references to obj from d
                    del d[obj]
                else:
                    d[obj] = new
            elif obj != dimensionless and d[obj]:
                # d[obj] is a singleton list.  Create a new reference and
                # store it instead (a new reference is needed to avoid
                # altering the linked list from self.unit).
                d[obj] = d[obj][0].clone()
        # Note d must have at least one key, namely dimensionless
        # If dimensionless is referenced in addition to other units,
        # remove the references to dimensionless.
        # But remember the multipliers!
        m = reduce(operator.mul,
                   map(lambda u: u.get_multiplicative_factor(),
                       d[dimensionless]),
                   1)
        if m == 1:
            del d[dimensionless]
        else:
            # Retain a single reference to dimensionless, with the
            # combined multiplier.
            d[dimensionless] = d[dimensionless][0].clone()
            d[dimensionless].multiplier = unicode(m)
##            print "Keeping d'less ref with m =",m,"from",
##            print self.description(),
##            if other_units:
##                print "and",other_units.description()
##            else:
##                print
        if not d:
            # The units definition only referenced dimensionless, and
            # the multiplier was 1, so we can replace it by dimensionless
            # itself
            new_units = dimensionless
        else:
            # Avoid creating new <units> elements unnecessarily
            unit_elements = set(d.values())
            if unit_elements == our_unit_elements:
                new_units = self
            elif unit_elements == other_unit_elements:
                new_units = other_units
            else:
                # Create new <units> element
                self.units_name_counter[0] += 1
                uname = u'___units_' + str(self.units_name_counter[0])
##                print ".simplify", uname
                new_units = self.xml_create_element(
                    u'units', NSS[u'cml'], attributes={u'name': uname})
                new_units._cml_generated = True
                # Set its parent to be ours or other_units'
                new_units.xml_parent = self._best_parent(other_units)
                # Add <unit> children
                for unit in unit_elements:
                    new_units.xml_append(unit)
                if self.model._is_new_units_obj(new_units):
                    # Add name to units dictionary
##                    print "Adding",uname,hash(new_units),new_units.description()
                    new_units.xml_parent.add_units(uname, new_units)
        # Cache and return result
        new_units = self.model._get_units_obj(new_units)
        self._cml_simplified[(other_units, other_exponent)] = new_units
        return new_units

    def _best_parent(self, other_units):
        """Return a suitable parent for units formed from self and other_units.

        If either constituent is in a component, that component should be the
        parent, otherwise the model should be.
        """
        p1, p2 = self.xml_parent, other_units and other_units.xml_parent
        if p2 and p1 != p2 and isinstance(p1, cellml_model):
            p1 = p2
        return p1

    def quotient(self, other_units):
        """Form the quotient of two units definitions.

        This method does not simplify the resulting units.
        Quotient units will be cached.
        """
        if not other_units in self._cml_quotients:
            # Create new <units> element
            self.units_name_counter[0] += 1
            uname = u'___units_' + str(self.units_name_counter[0])
            quotient_units = self.xml_create_element(
                u'units', NSS[u'cml'], attributes={u'name': uname})
            quotient_units._cml_generated = True
            quotient_units.xml_parent = self._best_parent(other_units)
            # Invent references to the constituent units
            u = self.xml_create_element(u'unit', NSS[u'cml'],
                                        attributes={u'units': self.name})
            u._set_units_element(self)
            quotient_units.xml_append(u)
            u = self.xml_create_element(u'unit', NSS[u'cml'],
                                        attributes={u'units': self.name,
                                                    u'exponent': u'-1'})
            u._set_units_element(other_units)
            quotient_units.xml_append(u)
            # Cache
            if self.model._is_new_units_obj(quotient_units):
                quotient_units.xml_parent.add_units(uname, quotient_units)
            quotient_units = self.model._get_units_obj(quotient_units)
            self._cml_quotients[other_units] = quotient_units
        return self._cml_quotients[other_units]

    def dimensionally_equivalent(self, other_units):
        """Return True iff other_units is dimensionally equivalent to self.
        
        As per appendix C.2.2, two units definitions have dimensional
        equivalence if, when each is recursively expanded and
        simplified into a product of base units, each has the same set
        of base units with the same exponent on corresponding base units.
        """
        u1 = self.expand().simplify()
        if isinstance(other_units, UnitsSet):
            other_units = other_units.extract()
        u2 = other_units.expand().simplify()
        # Build dictionaries mapping base_unit_obj to exponent.
        d1, d2 = {}, {}
        if u1.is_base_unit():
            d1[u1] = 1
        else:
            for u in u1.unit:
                d1[u.get_units_element()] = u.get_exponent()
        if u2.is_base_unit():
            d2[u2] = 1
        else:
            for u in u2.unit:
                d2[u.get_units_element()] = u.get_exponent()
        # Compare keys: do u1 & u2 have the same set of base units?
        sym_diff = set(d1.keys()) ^ set(d2.keys())
        if sym_diff:
            dimensionless = self.get_units_by_name(u'dimensionless')
            if not sym_diff == set([dimensionless]):
                # Symmetric difference is non-empty, but not an ignorable
                # instance of dimensionless
                return False
        # Compare corresponding exponents
        for k in d1:
            try:
                if d1[k] != d2[k]: return False
            except KeyError:
                # One may have dimensionless as a key
                pass
        # We have a match!
        return True


class cellml_unit(element_base):
    """Specialised class for <unit> elements.
    
    Maintains a reference to the object representing the units definition
    it references, provides some helpful accessor type methods, and allows
    safe, easy cloning of <unit> elements.
    """
    
    def __init__(self):
        element_base.__init__(self)
        self._cml_units_obj = None
        return

    def _hash_tup(self):
        """Create a tuple to be used for hashing/equality tests."""
        return (self.get_units_element(), getattr(self, u'prefix_', ''),
                self.get_multiplier(), self.get_exponent(),
                self.get_offset())

    def __eq__(self, other):
        """Compare two <unit> elements.

        Two <unit> elements are equal if they reference the same <units>
        element, and have the same prefix, multiplier, etc.
        """
        eq = False
        if isinstance(other, cellml_unit):
            eq = self._hash_tup() == other._hash_tup()
        return eq
    def __ne__(self, other):
        """The inverse of self.__eq__(other)."""
        return not self.__eq__(other)
    def __hash__(self):
        """Richer hashing function than the default based on object id.

        Returns the hash of a tuple of relevant attributes."""
        return hash(self._hash_tup())

    def get_units_element(self):
        """
        Return the object representing the <units> element that this <unit> element references.
        """
        if self._cml_units_obj is None:
            # Chase the reference and cache it
            self._cml_units_obj = self.xml_parent.get_units_by_name(self.units)
        return self._cml_units_obj
    def _set_units_element(self, obj, override=False):
        """
        Set the object representing the <units> element that this <unit> element references.
        
        Don't use unless you know what you're doing.
        """
        assert override or self._cml_units_obj is None
        self._cml_units_obj = obj
        return

    SI_PREFIXES = {"yotta": 24, "zetta": 21, "exa": 18, "peta": 15,
                   "tera": 12, "giga": 9, "mega": 6, "kilo": 3,
                   "hecto": 2, "deka": 1,
                   "deci": -1, "centi": -2,
                   "milli": -3, "micro": -6, "nano": -9, "pico": -12,
                   "femto": -15, "atto": -18, "zepto": -21, "yocto": -24}

    def get_multiplicative_factor(self):
        """Return the factor this units reference is multiplied by.
        
        Return the quantity m.p^e as a floating point number, where:
          m is the multiplier (default value 1.0)
          p is the multiplicative factor due to the prefix (default 10^0=1)
          e is the exponent (default 1)
        """
        m = self.get_multiplier()
        e = self.get_exponent()
        p = getattr(self, u'prefix_', 0) # Since prefix is a method :(
        if p in self.SI_PREFIXES:
            p = self.SI_PREFIXES[p]
        else:
            p = int(p) # RELAX NG schema says it's an integer
        p = 10**p
        return m * (p**e)

    def get_multiplier(self):
        """Return the multiplier of this units reference, as a float."""
        return float(getattr(self, u'multiplier', 1))

    def get_exponent(self):
        """Return the exponent of this units reference, as a float."""
        return float(getattr(self, u'exponent', 1))

    def get_offset(self):
        """Return the offset of this units reference, as a float."""
        return float(getattr(self, u'offset', 0))

    def clone(self):
        """Clone this object.
        
        Return a new <unit> element that has the same attributes as this
        one.
        """
        attrs = {}
        for apyname, aname in self.xml_attributes.iteritems():
            attrs[aname] = getattr(self, apyname)
        new = self.xml_create_element(u'unit', NSS[u'cml'], attributes=attrs)
        if self._cml_units_obj:
            new._set_units_element(self._cml_units_obj)
        new.xml_parent = self.xml_parent
        return new


class EvaluationError(Exception):
    """
    Exception class for errors raised trying to evaluate a MathML expression.
    """
    pass

class MathsError(Exception):
    """
    Exception class for validation errors raised while checking mathematics.
    """
    def __init__(self, context_obj, message, warn=False, level=None):
        """Create a mathematics validation error.
        
        context_class should be the object that is reporting the error.
        message gives more details on what went wrong.
        If warn is set to true then produce a warning message, not an error.
        level gives the level of the message logged.
        """
        self.message = message
        self.warn = warn
        self.level = level or (logging.ERROR,logging.WARNING)[warn]
        self.show_xml_context = False
        self._show_xml_context_only = False

        self.cname = context_obj.component.name
        self.ename = context_obj.localName
        self.context = context_obj
        
        # Nicer context explanation
        if isinstance(context_obj, cellml_variable):
            self.expr_index = self.math_index = 0
            self.reaction_spec = ''
            return
        expr_root = context_obj.xml_xpath(u'ancestor-or-self::*[local-name(..)="math"]')[0]
        self.expr_index = self.math_index = 1
        math = expr_root.xml_parent
        for elt in math.xml_element_children():
            if elt is expr_root: break
            if elt.localName in [u'apply', u'semantics']:
                self.expr_index += 1
        for elt in math.xml_element_children(math.xml_parent):
            if elt is math: break
            if elt.localName == u'math':
                self.math_index += 1
        # Is this in a reaction?
        vref = context_obj.xml_xpath(u'ancestor::cml:variable_ref')
        if vref:
            self.reaction_spec = ' in mathematics for variable "%s" in a reaction' % vref[0].variable
        else:
            self.reaction_spec = ''
        return

    def show_xml_context_only(self):
        """Only show the XML where the error occurred."""
        self.show_xml_context = True
        self._show_xml_context_only = True

    def __str__(self):
        msg = unicode(self)
        return msg.encode('UTF-8')

    def ordinal(self, i):
        "Convert integer i to an ordinal string."
        if i // 10 == 1:  suf = 'th'
        elif i % 10 == 1: suf = 'st'
        elif i % 10 == 2: suf = 'nd'
        elif i % 10 == 3: suf = 'rd'
        else:             suf = 'th'
        return "%d%s" % (i, suf)

    def _generate_message(self, where):
        if self.warn: type = 'Warning'
        else: type = 'Error'
        msg = [type, ' ', where, ': ', self.message]
        if not self._show_xml_context_only:
            msg.extend(['\n  Context: ', self.ordinal(self.expr_index),
                         ' expression in the ', self.ordinal(self.math_index),
                         ' math element', self.reaction_spec,
                         ' in component ', self.cname, '\n  XPath: ',
                        element_xpath(self.context)])
        msg = u''.join(msg)
        if self.show_xml_context:
            # Return the context XML tree as well.
            xml = self.context.xml(indent = u'yes',
                                   omitXmlDeclaration = u'yes')
            msg = msg + u'\n' + unicode(xml, encoding='UTF-8')
        return msg

    def __unicode__(self):
        return self._generate_message('checking mathematics')

    
class UnitsError(MathsError):
    """
    Exception class for validation errors raised during units checking
    of mathematics.
    """
    def __init__(self, context_obj, message, warn=False, level=None):
        """Create a units validation error.
        
        context_class should be the object that is reporting the error.
        message gives more details on what went wrong.
        If warn is set to true then produce a warning message, not an error.
        level gives the level of the message logged.
        """
        MathsError.__init__(self, context_obj, message, warn=warn, level=level)
        return

    def __unicode__(self):
        return self._generate_message('checking units')
        

def child_i(elt, i):
    "Return the i'th child element of elt.  Indexing starts at 1."
    j = 0
    for e in elt.xml_children:
        if getattr(e, 'nodeType', None) == Node.ELEMENT_NODE:
            j += 1
            if j == i: return e
    else:
        # Not found :(
        raise ValueError("<"+elt.localName+"> does not have "+str(i)+
                         " child element(s).")
def _child1(elt):
    "Get the first child element of elt."
    return child_i(elt, 1)


######################################################################
#                          MathML elements                           #
######################################################################

class mathml_units_mixin(object):
    """Base class for units mixin classes."""
    def _add_units_conversion(self, expr, defn_units, to_units, no_act=False):
        """Add mathematics for an explicit units conversion.

        Wraps expr in the expression
        m[to_units/defn_units]*(expr-o1[defn_units]) + o2[to_units].
        """
#         print '_add_units_conv for', element_xpath(expr), 'from', defn_units.description(), 'to', to_units.description()
        if hasattr(expr.model, '_cml_special_units_converter') and not defn_units.dimensionally_equivalent(to_units):
            # This may be a special conversion case defined by a functional curation protocol
            if no_act:
                model._cml_conversions_needed = True
                return
            else:
                expr = expr.model._cml_special_units_converter(expr, defn_units, to_units)
#             print 'post special, expr units=', expr.get_units().description(), 'm=', expr.get_units().extract().expand().simplify().get_multiplicative_factor()
#             print 'orig defn_units=', defn_units.description(), 'm=', defn_units.expand().simplify().get_multiplicative_factor()
#             print 'orig to_units=', to_units.description(), 'm=', to_units.expand().simplify().get_multiplicative_factor()
            try:
                defn_units = expr.get_units().extract(check_equality=True)
            except:
                print 'ouch', expr.xml()
                for u in expr.get_units():
                    print u.description(), u.get_multiplier(), expr.get_units()._get_sources(u)
                raise
        defn_units_exp = defn_units.expand().simplify()
        to_units_exp = to_units.expand().simplify()
        # Conversion factor
        m = (defn_units_exp.get_multiplicative_factor() / to_units_exp.get_multiplicative_factor())
        # Replace expr by  m[to_units/defn_units]*(expr-o1[defn_units]) + o2[to_units]
        orig_expr, parent = expr, expr.xml_parent
        dummy = expr.xml_create_element(u'dummy', NSS[u'm'])
        model = expr.model # So we still have a reference after the next line
        parent.replace_child(expr, dummy) # Mark where to put the new elt
        if defn_units_exp.get_offset() != 0:
            # Create expr-o1 expression
            uattr = orig_expr._ensure_units_exist(defn_units, no_act=no_act)
            new_expr = mathml_apply.create_new(expr, u'minus',
                                               [expr, (unicode(defn_units_exp.get_offset()), uattr)])
            new_expr._cml_units = defn_units
            expr = new_expr
        if m != 1:
            quotient_units = to_units.quotient(defn_units)
            # Add units element to model if needed
            uattr = orig_expr._ensure_units_exist(quotient_units, no_act=no_act)
            # Create m*expr expression
            new_expr = mathml_apply.create_new(expr, u'times', [(unicode(m), uattr), expr])
            new_expr._cml_units = to_units
            expr = new_expr
        if to_units_exp.get_offset() != 0:
            # Create expr+o2 expression
            uattr = orig_expr._ensure_units_exist(to_units, no_act=no_act)
            new_expr = mathml_apply.create_new(expr, u'plus',
                                               [expr, (unicode(to_units_exp.get_offset()), uattr)])
            new_expr._cml_units = to_units
            expr = new_expr
        # Note that the model needed conversions
        if expr is not orig_expr:
            model._cml_conversions_needed = True
        if no_act:
            expr = orig_expr
        parent.replace_child(dummy, expr)
        return
        
    def _set_element_in_units(self, elt, units, no_act=False):
        """Try to set the units of the given element.

        Generates a debug message if this isn't possible.
        """
        if hasattr(elt, '_set_in_units') and callable(elt._set_in_units):
            elt._set_in_units(units, no_act)
        elif elt.localName in [u'false', u'true']:
            boolean = self.component.get_units_by_name('cellml:boolean')
            if boolean is not units:
                # TODO: *blink* this should never happen
                self._add_units_conversion(elt, boolean, units, no_act)
        elif elt.localName in [u'notanumber', u'pi', u'infinity', u'exponentiale']:
            dimensionless = self.component.get_units_by_name('dimensionless')
            if dimensionless is not units:
                # TODO: *blink* this should never happen
                self._add_units_conversion(elt, dimensionless, units, no_act)
        else:
            DEBUG('validator',
                  'Cannot set units (to', units.description(), ') for element', elt.localName)
        return

class mathml_units_mixin_tokens(mathml_units_mixin):
    """Contains the _set_in_units method for ci, cn, etc."""
    def _set_in_units(self, units, no_act=False):
        """Set the units this element should be expressed in.

        Where these aren't the units it's defined in, replace self by suitable units conversion mathematics.
        """
        defn_units = self.get_units(return_set=False)
        if defn_units != units:
            self._add_units_conversion(self, defn_units, units, no_act)
        # Store the units
        if not no_act:
            self._cml_units = units
        return

class mathml_units_mixin_set_operands(mathml_units_mixin):
    def _set_in_units(self, units, no_act=False):
        """Set the units of the application of this operator.

        The default behaviour for many operators is to simply set all operands to have the given units.
        """
        # TODO: Do the conversion at this level sometimes rather than pushing it down the tree?
        app = self.xml_parent
        # We need to convert the operands to a list, because the tree
        # may be modified if a conversion is thought to be needed
        for operand in list(app.operands()):
            self._set_element_in_units(operand, units, no_act)
        # And set our units to what they now are
        if not no_act:
            app._cml_units = units
        return

class mathml_units_mixin_equalise_operands(mathml_units_mixin):
    def _set_in_units(self, units, no_act=False):
        """Set the units of the application of this operator.

        This method is used for the relational operators.  It ignores the
        given units, and instead ensures that all operands have the same
        units.  The units it chooses are those that are 'least' amongst the
        possibilities for the operand units, i.e. that have the smallest
        multiplicative factor when expanded.
        """
        app = self.xml_parent
        min_factor, operand_units = None, None
        for us in app._get_operand_units():
            if isinstance(us, cellml_units):
                us = [us]
            for u in us:
                f = u.expand().get_multiplicative_factor()
                if f < min_factor or min_factor is None:
                    min_factor = f
                    operand_units = u
        # We need to convert the operands to a list, because the tree
        # may be modified if a conversion is thought to be needed
        for operand in list(app.operands()):
            # TODO: Modify tree to collect same conversions together?
            self._set_element_in_units(operand, operand_units, no_act)
        # And set our units to what they now are
        if not no_act:
            app._cml_units = units
        return

class mathml_units_mixin_choose_nearest(mathml_units_mixin):
    def _set_in_units(self, desired_units, no_act=False):
        """Set the units of the application of this operator.

        This mixin is used for the <times> and <divide> operators.
        There are 2 possible strategies here.  One is to pick one of the
        operands and convert it so that the overall units match those
        required.  The other is to pick units from the set of those
        possible for this application, and convert the result to the
        desired units.  We go with the latter option, picking the units
        that are least in the sense that they have the least multipicative
        factor, but where possible that factor is no less than that on the
        desired units.
        """
        app = self.xml_parent
        min_factor, best_factor = None, None
        least_units, best_units = None, None
        desired_factor = desired_units.expand().get_multiplicative_factor()
        DEBUG('validator', '>',self.localName,':',desired_factor,
              desired_units.description())
        units_set = app.get_units().get_consistent_set(desired_units)
        for possible_units in units_set:
            f = possible_units.expand().get_multiplicative_factor()
            if min_factor is None or f<min_factor:
                least_units, min_factor = possible_units, f
            if f >= desired_factor and (best_factor is None or f < best_factor):
                best_units, best_factor = possible_units, f
        if best_units is None:
            # All factors were less than that desired, so just go with the least
            units = least_units
        else:
            units = best_units
        DEBUG('validator', '\t<-',
              units.expand().get_multiplicative_factor(),
              units.description())
        # Add units conversion code
        app._add_units_conversion(app, units, desired_units, no_act)
        # Set the operand units
        for src_units_set, src_units in app.get_units()._get_sources(units):
            expr = src_units_set.get_expression()
            DEBUG('validator', '#',self.localName,':',
                  src_units.description(),expr.localName)
            self._set_element_in_units(expr, src_units, no_act)
        # Record which units we used
        if not no_act:
            app._cml_units = units
        return
    
class mathml_units_mixin_container(mathml_units_mixin):
    def _set_in_units(self, units, no_act=False):
        """Set the units of this element.

        For container elements, we set the units of the child(ren).
        """
        # We need to copy the children list, because the tree
        # may be modified if a conversion is thought to be needed
        for elt in self.xml_children[:]:
            if getattr(elt, 'nodeType', None) == Node.ELEMENT_NODE:
                self._set_element_in_units(elt, units, no_act)
        # And set our units to what they now are
        if not no_act:
            self._cml_units = units
        return

class mathml(element_base):
    """Base class for MathML elements."""
    def __init__(self):
        super(mathml, self).__init__()
        self._cml_component = None
        self._cml_model = None
        self._cml_source_expr_for_clone = None
        return
    
    def __repr__(self):
        return '<%s (%s) at 0x%x>' % (self.__class__.__name__, str(self), id(self))

    def __deepcopy__(self, memo):
        """Customise copying of MathML expressions.

        When copying an expression tree, we only want to deepcopy the
        XML, not the additional info that we have attached to it -
        these should be copied as references to the originals.
        """
        new_elt = copy.copy(self)
        # Children may refer to us, so need to update memo before copying children
        assert id(self) not in memo
        memo[id(self)] = new_elt
        new_dict = {}
        for name, value in self.__dict__.iteritems():
            name_copy = copy.deepcopy(name, memo)
            if not name.startswith('_cml'):
                new_dict[name_copy] = copy.deepcopy(value, memo)
            else:
                new_dict[name_copy] = value
        new_elt.__dict__.update(new_dict)
        return new_elt

    @staticmethod
    def clone(expr):
        """Properly clone a MathML sub-expression.

        Makes sure siblings and parent don't get copied too.
        """
#        print "Cloning MathML", prid(expr, True)
        next_elem, par = expr.next_elem, getattr(expr, 'xml_parent', None)
        expr.next_elem = None  # Make sure we don't copy siblings...
        expr.xml_parent = None # ...and parent
        new_expr = copy.deepcopy(expr) # Do the clone
        expr.next_elem = next_elem # Restore siblings...
        expr.xml_parent = par      # ...and parent to original
        return new_expr
    
    def clone_self(self, register=False):
        """Properly clone this expression.
        
        If register is True, then keep a link to this expression in the clone.
        """
        clone = mathml.clone(self)
        if register:
            clone._cml_source_expr_for_clone = self
        return clone
    
    def get_original_of_clone(self):
        """
        If this is a clone with a registered original expression, return it.
        Otherwise returns None.
        """
        return self._cml_source_expr_for_clone

    def get_component(self):
        "Cache & return the enclosing component element."
        if self._cml_component is None:
            def get_ancestor(elt, name):
                while elt and elt.localName != name:
                    elt = elt.parentNode
                return elt
            comp = get_ancestor(self, u'component')
            if comp:
                self._cml_component = comp
            else:
                # It may be in the solver_info section, in which case fake a component
                solver_info = get_ancestor(self, u'solver_info')
                if solver_info:
                    self._cml_component = self.model.get_component_by_name(u'')
                else:
                    raise ValueError("MathML element " + str(self) + " does not occur in a model component!")
        return self._cml_component
    component = property(get_component)

    def _unset_cached_links(self, elt=None):
        """Forget cached component and variable references in this MathML tree.
        
        Used by partial evaluator when moving maths to a new component, and by
        simulation protocols.
        """
        if elt is None:
            elt = self
        if isinstance(elt, mathml):
            elt._cml_component = None
        for child in self.xml_element_children(elt):
            if hasattr(child, '_unset_cached_links'):
                child._unset_cached_links()
            else:
                self._unset_cached_links(child)
        return

    @property
    def model(self):
        """Cache & return the enclosing model element."""
        if self._cml_model is None:
            self._cml_model = self.rootNode.model
        return self._cml_model
    
    def eval(self, elt):
        """Evaluate the given element.

        Tries to evaluate the given element, and raises an EvaluationError
        if this is not possible.
        """
        if hasattr(elt, 'evaluate') and callable(elt.evaluate):
            return elt.evaluate()
        elif elt.localName == u'pi':
            return math.pi
        else:
            # No evaluate() method on elt
            raise EvaluationError("Don't know how to evaluate element " +
                                  elt.localName)
        return

    def _ensure_units_exist(self, units=None, no_act=False):
        """Ensure that there is an element in the XML tree giving this
        expression's units.

        Add a new <units> element if this expression has generated units.

        If units is not None, use the given units rather than those of
        this expression.

        Return an attribute dictionary with the appropriate units
        attribute."""
        if no_act:
            # Doesn't matter what we return, as it wont be used
            return {(u'cml:units', NSS[u'cml']): u'#UNUSED#'}
        try:
            if units is None:
                units = self.get_units().extract()
##            _u = units
            units = self.model._get_units_obj(units)
            if units._cml_generated and units.name[:3] == "___":
##                print "Adding",units.name, hash(units), units.description(),
##                print "(was",id(_u),"now",id(units),")"
                # Ensure referenced units exist
                for unit in getattr(units, u'unit', []):
                    self._ensure_units_exist(unit.get_units_element())
                    unit.units = unit.get_units_element().name
                # Rename units and add to XML tree
                msg = "Adding units " + units.name + " as "
                units.name = units.description(cellml=True)
                msg = msg + units.name
                DEBUG('partial-evaluator', msg.encode('UTF-8'))
                if units.name == units.unit.units:
                    # Uh-oh
                    DEBUG('partial-evaluator', 'Generated units',
                          units.name, 'identical to referenced units; ignoring.')
                    assert units.get_multiplicative_factor() == 1
                    assert units.get_offset() == 0
                else:
                    units.xml_parent.add_units(units.name, units)
                    units.xml_parent.xml_append(units)
            attrs = {(u'cml:units', NSS[u'cml']): units.name}
        except UnitsError:
            # Hack to allow PE on broken (wrt units) models
            attrs = {(u'cml:units', NSS[u'cml']): u'#FUDGE#'}
        return attrs

    def varobj(self, ci_elt):
        """Return the variable object for the given ci element.

        This method is more general than ci_elt.variable, working even
        for ci elements outside of a component.  Such elements *must*
        contain a fully qualified variable name, i.e. including the
        name of the component the variable occurs in.  This method
        handles a variety of encodings of variable names that contain
        the component name.
        """
        try:
            var = ci_elt.variable
        except:
            varname = unicode(ci_elt).strip()
            var = cellml_variable.get_variable_object(self.model, varname)
        return var

    def vars_in(self, expr):
        """Return a list of 'variable' objects used in the given expression.

        This method doesn't make use of the dependency information
        generated when validating the model, but parses the
        mathematics afresh.  It is used to refresh the dependency
        lists after partial evaluation, and to determine dependencies
        in mathematics added outside the normal model structure
        (e.g. Jacobian calculation).

        If an ODE appears, includes the mathml_apply instance defining
        the ODE.  Otherwise all returned objects will be
        cellml_variable instances.
        """
        res = set()
        if isinstance(expr, mathml_ci):
            res.add(self.varobj(expr))
        elif isinstance(expr, mathml_apply) and \
                 expr.operator().localName == u'diff':
            dep_var = self.varobj(expr.ci)
            indep_var = self.varobj(expr.bvar.ci)
            res.add(dep_var.get_ode_dependency(indep_var))
        elif hasattr(expr, 'xml_children'):
            for child in expr.xml_children:
                res.update(self.vars_in(child))
        return res
    
    def same_tree(self, other, this=None):
        """Check whether this element represents the same tree as a given element."""
        if this is None: this = self
        equal = (this.localName == other.localName
                 and len(getattr(this, 'xml_children', [])) == len(getattr(other, 'xml_children', [])))
        if equal and this.localName in [u'cn', u'ci']:
            equal = unicode(this) == unicode(other)
        if equal and hasattr(this, 'xml_children'):
            for tc, oc in zip(self.xml_element_children(this), self.xml_element_children(other)):
                if not self.same_tree(oc, tc):
                    equal = False
                    break
        return equal

    def _xfer_complexity(self, new_elt):
        """Transfer our complexity to the new element.
        
        PE is replacing us by a new element.  If we are annotated with
        a complexity - the complexity of this expression prior to PE -
        then transfer the annotation to the new element.
        """
        try:
            new_elt._cml_complexity = self._cml_complexity
        except AttributeError:
            pass
        return
    def _adjust_complexity(self, old_elt, new_elt):
        """Adjust ancestor complexity because old_elt changed to new_elt.

        The purpose of this method is to allow us to keep track of
        what the complexity of each expression node *was* prior to PE
        being performed.  Thus we cannot just re-compute complexities,
        but must update values using the original complexities.  If a
        variable definition is instantiated, then the complexity of
        the expression containing the lookup must be adjusted to
        reflect the additional expense of evaluating the defining
        expression.

        When this method is called, only new_elt is a child of self.
        """
        #print "Adjusting", element_xpath(self), "due to", element_xpath(old_elt),
        #if isinstance(old_elt, mathml_ci):
        #    print unicode(old_elt)
        #else:
        #    print
        try:
            new = new_elt._cml_complexity
            old = old_elt._cml_complexity
        except AttributeError:
            return
        #print "  by", new-old
        elt = self
        while elt:
            if isinstance(elt, mathml_piecewise):
                # Piecewise is tricky to adjust!  So we 're-compute' instead.
                ac, piece_ac = 0, []
                for piece in getattr(elt, u'piece', []):
                    ac += child_i(piece, 2)._cml_complexity
                    piece_ac.append(child_i(piece, 1)._cml_complexity)
                if hasattr(elt, u'otherwise'):
                    piece_ac.append(child_i(elt.otherwise, 1)._cml_complexity)
                ac += max(piece_ac)
                elt._cml_complexity = ac
            elif hasattr(elt, '_cml_complexity'):
                elt._cml_complexity += (new - old)
            elt = getattr(elt, 'xml_parent', None)
        return
    
    def classify_child_variables(self, elt, **kwargs):
        """Classify variables in the given expression according to how they are used.
        
        In the process, compute and return a set of variables on which that expression depends.

        If dependencies_only then the variable classification will not be
        done, only dependencies will be analysed.  This is useful for doing
        a 'light' re-analysis if the dependency set has been reduced; if the
        set has increased then the topological sort of equations may need to
        be redone.
        
        The function needs_special_treatment may be supplied to override the
        default recursion into sub-trees.  It takes a single sub-tree as
        argument, and should either return the dependency set for that
        sub-tree, or None to use the default recursion.  This is used when
        re-analysing dependencies after applying lookup tables, since table
        lookups only depend on the keying variable.
        """
        if hasattr(elt, 'classify_variables'):
            dependencies = elt.classify_variables(**kwargs)
        else:
            dependencies = set()
            needs_special_treatment = kwargs.get('needs_special_treatment', lambda e: None)
            for e in elt.xml_element_children():
                child_deps = needs_special_treatment(e)
                if child_deps is None:
                    child_deps = self.classify_child_variables(e, **kwargs)
                dependencies.update(child_deps)
        return dependencies

class mathml_math(mathml):
    def __init__(self):
        super(mathml_math, self).__init__()
        return

class mathml_constructor(mathml):
    """
    Base class for MathML constructor elements, e.g. apply and piecewise.
    """
    def __init__(self):
        super(mathml_constructor, self).__init__()
        return

    def _tree_complexity(self, elt, **kw):
        """
        Calculate a rough estimate of the computation time for
        evaluating the given element.
        
        If lookup_tables is True, then assume we're using lookup tables
        where possible.
        If store_result is True, the complexity is saved to the
        _cml_complexity attribute.
        If algebraic is True, the complexity is calculated as a dictionary,
        mapping node types to the number of occurences of that type.
        """
        kw['lookup_tables'] = kw.get('lookup_tables', False)
        kw['store_result'] = kw.get('store_result', False)
        kw['algebraic'] = kw.get('algebraic', False)
        if kw['algebraic']:
            ac = {}
        if kw['lookup_tables'] and \
             elt.getAttributeNS(NSS['lut'], u'possible', u'no') == u'yes':
            # This elt will be replaced by a lookup table
            if hasattr(self.rootNode, 'num_lookup_tables'):
                # Count the number of used lookup tables
                self.rootNode.num_lookup_tables += 1
            # Cost per table: 2 lookups, 2 +, -, *, 3 ci
            if kw['algebraic']:
                ac['lookup'] = 2
                ac['op'] = 3
                ac['times'] = 1
                ac['variable'] = 3
            else:
                ac = 2*1 + 2*1 + 1 + 1 + 3*0.7
        elif hasattr(elt, 'tree_complexity') \
           and callable(elt.tree_complexity):
            ac = elt.tree_complexity(**kw)
        elif elt.localName in ['true', 'false', 'cn', 'exponentiale', 'pi']:
            if kw['algebraic']: ac['constant'] = 1
            else: ac = 0.5
        elif elt.localName == 'ci':
            if kw['algebraic']: ac['variable'] = 1
            else: ac = 0.7
        elif elt.localName in ['degree', 'logbase']:
            ac = self._tree_complexity(child_i(elt, 1), **kw)
        else:
            raise EvaluationError("Don't know complexity of element " +
                                  elt.localName)
        if kw['store_result'] and isinstance(elt, mathml):
            elt._cml_complexity = ac
        return ac

    def _get_element_binding_time(self, elt):
        """Helper method to get the binding time of a MathML element."""
        if hasattr(elt, '_get_binding_time'):
            # Call the element's method
            return elt._get_binding_time()
        elif elt.localName in [u'true', u'false', u'exponentiale', u'pi']:
            return BINDING_TIMES.static
        else:
            raise EvaluationError("Don't know how to compute binding time"
                                  " of element " + elt.localName)

    def _get_element_units(self, elt, return_set=True):
        """Helper method to get the units of a MathML element."""
        if hasattr(elt, 'get_units'):
            # Element has a method to get its units, so call it
            u = elt.get_units()
        elif hasattr(elt, '_cml_units') and elt._cml_units:
            # Units have been cached
            u = elt._cml_units
        else:
            # Let's figure it out ourselves...
            if elt.localName in [u'false', u'true']:
                u = UnitsSet(
                    [self.component.get_units_by_name('cellml:boolean')],
                    expression=elt)
            elif elt.localName in [u'notanumber', u'pi', u'infinity',
                                   u'exponentiale']:
                u = UnitsSet(
                    [self.component.get_units_by_name('dimensionless')],
                    expression=elt)
            else:
                # Unknown or unexpected element
                raise UnitsError(self, u''.join([
                    u'Unsupported element "', elt.localName, '".']))
        if return_set:
            assert isinstance(u, UnitsSet)
        else:
            u = u.extract()
        return u

    def _reduce_elt(self, elt):
        """Try to reduce the given element.

        Call the _reduce method on elt, if it has one.
        If not, do nothing (we assume elt cannot be reduced).
        """
        if hasattr(elt, '_reduce') and callable(elt._reduce):
            elt._reduce()
        else:
            DEBUG('partial-evaluator', "Don't know how to reduce",
                  elt.localName)
        return

    def _eval_self(self):
        """Evaluate self and return <cn>, <true> or <false>, as appropriate."""
        value = self.evaluate()
        if value is True:
            new_elt = self.xml_create_element(u'true', NSS[u'm'])
        elif value is False:
            new_elt = self.xml_create_element(u'false', NSS[u'm'])
        else:
            # Add a new <units> element to the document if needed
            attrs = self._ensure_units_exist()
            new_elt = self.xml_create_element(u'cn', NSS[u'm'],
                                              content=unicode("%.17g" % value),
                                              attributes=attrs)
        return new_elt

    def _update_usage_counts(self, expr, remove=False):
        """Update usage counts of variables used in the given expression.

        By default, increment the usage count of any variable occuring
        in a <ci> element within expr.  If remove is set to True,
        then decrement the usage counts instead.
        """
        if isinstance(expr, mathml_ci):
            if remove:
                expr.variable._decrement_usage_count()
            else:
                raise NotImplementedError("_update_usage_counts currently only reliable for remove=True")
                expr.variable._used()
        elif isinstance(expr, mathml_apply) and isinstance(expr.operator(),
                                                           mathml_diff):
            # TODO: Check if this is a suitable handling of ODEs on a RHS
            # It matches the current behaviour in apply.classify_variables.
            pass
        else:
            for e in self.xml_element_children(expr):
                self._update_usage_counts(e, remove=remove)

class mathml_cn(mathml, mathml_units_mixin_tokens):
    def __init__(self):
        super(mathml_cn, self).__init__()
        self._cml_units = None
        return

    def evaluate(self):
        """
        Convert the text content of this element to a floating point
        value and return it.  Will handle the type attribute and, if
        relevant to the type, the sep child element, but does not yet
        handle the base attribute.
        """
        if hasattr(self, u'base'):
            raise ValueError('pycml does not yet support the base attribute on cn elements')
        if hasattr(self, u'type'):
            if self.type == u'real':
                val = float(unicode(self))
            elif self.type == u'integer':
                val = int(unicode(self))
            elif self.type == u'e-notation':
                assert len(self.xml_children) == 3
                assert self.xml_children[1] is self.sep
                mantissa = unicode(self.xml_children[0]).strip()
                exponent = unicode(self.xml_children[2]).strip()
                val = float(mantissa + 'e' + exponent)
            elif self.type == u'rational':
                assert len(self.xml_children) == 3
                assert self.xml_children[1] is self.sep
                numer = int(unicode(self.xml_children[0]))
                denom = int(unicode(self.xml_children[2]))
                val = numer / denom
            else:
                raise ValueError('Unsupported type attribute for cn element: '
                                 + self.type)
        else:
            val = float(unicode(self))
        return val

    def _get_binding_time(self):
        """Return the binding time of this expression.

        The binding time of a <cn> element is always static,
        unless the CellML is annotated otherwise.
        """
        bt = self.getAttributeNS(NSS['pe'], u'binding_time', u'static')
        return getattr(BINDING_TIMES, bt)
    
    def _reduce(self):
        """Reduce this expression by evaluating its static parts.
        
        Is actually a no-op; we must have been annotated explicitly as dynamic.
        """
        return

    def get_units(self, return_set=True):
        """Return the units this number is expressed in."""
        if not self._cml_units:
            self._cml_units = UnitsSet([self.component.get_units_by_name(self.units)], expression=self)
        if not return_set:
            u = self._cml_units.extract()
        else:
            if not isinstance(self._cml_units, UnitsSet):
                self._cml_units = UnitsSet([self._cml_units], expression=self)
            u = self._cml_units
        return u

    @staticmethod
    def create_new(elt, value, units):
        """Create a new <cn> element with the given value and units."""
        attrs = {(u'cml:units', NSS[u'cml']): unicode(units)}
        new_elt = elt.xml_create_element(u'cn', NSS[u'm'],
                                         attributes=attrs,
                                         content=unicode(value))
        return new_elt

class mathml_ci(mathml, mathml_units_mixin_tokens):
    def __init__(self):
        super(mathml_ci, self).__init__()
        self._cml_variable = None
        self._cml_units = None
        return
    
    def _unset_cached_links(self, elt=None):
        """Forget cached component and variable references in this MathML tree.
        
        Used by partial evaluator when moving maths to a new component, and by
        simulation protocols.
        """
        self._cml_variable = None
        super(mathml_ci, self)._unset_cached_links()

    @property
    def variable(self):
        """Cache & return the variable object refered to by this element."""
        if self._cml_variable is None:
            vname = unicode(self).strip()
            self._rename(vname) # Remove excess whitespace from our text content
            self._cml_variable = self.component.get_variable_by_name(vname)
        return self._cml_variable
    
    def _set_variable_obj(self, var):
        """Set the variable object referred to by this element."""
        self._cml_variable = var

    def get_units(self, return_set=True):
        """Return the units of the variable represented by this element."""
        if not self._cml_units:
            self._cml_units = UnitsSet([self.component.get_units_by_name(self.variable.units)],
                                       expression=self)
        if not return_set:
            u = self._cml_units.extract()
        else:
            if not isinstance(self._cml_units, UnitsSet):
                self._cml_units = UnitsSet([self._cml_units], expression=self)
            u = self._cml_units
        return u
    
    def evaluate(self):
        """
        Evaluate this expression by returning the value of the
        variable it represents.
        """
        return self.variable.get_value()

    def _get_binding_time(self):
        """Return the binding time of this expression.

        The binding time of a <ci> element is that of the variable it
        represents.
        """
        return self.variable._get_binding_time()

    def _rename(self, new_name=None):
        """Update the variable reference to use a canonical name."""
        self.xml_remove_child(unicode(self))
        if new_name is None:
            new_name = self.variable.fullname(cellml=True)
        self.xml_append(unicode(new_name))
        return

    def _reduce(self):
        """Reduce this expression by evaluating its static parts.

        If this is a static variable, replace by its value (as a <cn> element).

        Otherwise the behaviour depends on the number of uses of this
        variable.  If there is only one, instantiate the definition of
        this variable here in place of the <ci> element, otherwise
        leave the element unchanged to avoid code duplication.
        """
        bt = self._get_binding_time()
        DEBUG('partial-evaluator', "Reducing", self.variable.fullname(),
              "which is", bt)
        if bt is BINDING_TIMES.static:
            value = self.evaluate()
            attrs = {(u'cml:units', NSS[u'cml']): self.variable.units}
            cn = self.xml_create_element(u'cn', NSS[u'm'],
                                         content=unicode("%.17g" % value),
                                         attributes=attrs)
            DEBUG('partial-evaluator', "  value =", unicode(cn))
            self._xfer_complexity(cn)
            self.xml_parent.xml_insert_after(self, cn)
            self.xml_parent.xml_remove_child(self)
            self.variable._decrement_usage_count()
        else:
            defns = self.variable.get_dependencies()
            if defns:
                defn = defns[0]
            else:
                # Just need to update the name to be canonical - done in later pass
                defn = None
            if isinstance(defn, cellml_variable):
                if self.variable.pe_keep:
                    # Don't remove this var, just reduce its source
                    DEBUG('partial-evaluator', "Keeping",
                          self.variable.fullname())
                    self.variable._reduce(update_usage=True)
                else:
                    # Create a new <ci> element
                    ci = self.xml_create_element(
                        u'ci', NSS[u'm'], content=defn.fullname(cellml=True))
                    self._xfer_complexity(ci)
                    ci._set_variable_obj(defn)
                    DEBUG('partial-evaluator', "  to", defn.fullname())
                    self.xml_parent.xml_insert_after(self, ci)
                    self.xml_parent.xml_remove_child(self)
                    # Decrement the usage count of just us, not source vars
                    self.variable._decrement_usage_count(follow_maps=False)
                    # May need to recurse down maps
                    ci._reduce()
            elif isinstance(defn, mathml_apply):
                if (not self.variable.pe_keep and
                    (self.variable.get_usage_count() == 1 or
                     self.rootNode.partial_evaluator.is_instantiable(defn))):
                    # defn is defining expression, so will be a MathML element already,
                    # and should be reduced already as well due to topological sort.
                    # Clone the RHS and instantiate it here.
                    rhs = mathml.clone(list(defn.operands())[1])
                    DEBUG('partial-evaluator', "  to", rhs)
                    parent = self.xml_parent
                    parent.xml_insert_after(self, rhs)
                    parent.xml_remove_child(self)
                    parent._adjust_complexity(self, rhs)
                    # Flag the defining expression for removal
                    defn._pe_process = u'remove'
                    self.variable._decrement_usage_count()
            elif defn is not None:
                raise ValueError("Unexpected variable definition: " + defn.xml())
        return
    
    def classify_variables(self, dependencies_only=False,
                           needs_special_treatment=lambda n: None):
        """Classify variables in this expression according to how they are used.
        
        For ci elements we just return a set containing the referenced variable
        as the single dependency.  If dependencies_only is False, we also mark
        the variable as used.
        """
        var = self.variable
        if not dependencies_only:
            var._used()
        return set([var])

    @staticmethod
    def create_new(elt, variable_name):
        """Create a new <ci> element with the given variable name."""
        new_elt = elt.xml_create_element(u'ci', NSS[u'm'],
                                         content=unicode(variable_name))
        return new_elt

class mathml_apply(Colourable, mathml_constructor, mathml_units_mixin):

    QUALIFIERS = frozenset(('degree', 'bvar', 'logbase',
                            'lowlimit', 'uplimit', 'interval', 'condition',
                            'domainofapplication', 'momentabout'))

    class OPS:
        """Classifications of operators."""
        absRound = frozenset(('abs', 'floor', 'ceiling', 'rem'))
        timesDivide = frozenset(('times', 'divide'))
        plusMinus = frozenset(('plus', 'minus'))
        trig = frozenset(('sin', 'cos', 'tan', 'sec', 'csc', 'cot',
                          'sinh', 'cosh', 'tanh', 'sech', 'csch', 'coth',
                          'arcsin', 'arccos', 'arctan',
                          'arcsec', 'arccsc', 'arccot',
                          'arcsinh', 'arccosh', 'arctanh',
                          'arcsech', 'arccsch', 'arccoth'))
        elementary = frozenset(('exp', 'log', 'ln')).union(trig)
        relations = frozenset(('eq', 'neq', 'gt', 'lt', 'geq', 'leq'))
        logical = frozenset(('and', 'or', 'xor', 'not'))
    
    def __init__(self):
        super(mathml_apply, self).__init__()
        self._cml_units = None
        self.clear_dependency_info()
    
    def clear_dependency_info(self):
        """Clear the type, dependency, etc. information for this equation.
        
        This allows us to re-run the type & dependency analysis for the model."""
        self._cml_binding_time = None
        # Dependency graph edges
        self._cml_depends_on = []
        self._cml_assigns_to = None
        self.clear_colour()

    def get_dependencies(self):
        """Return the list of variables this expression depends on."""
        return self._cml_depends_on

    def operator(self):
        """Return the element representing the operator being applied."""
        return _child1(self)

    def _is_qualifier(self, element):
        """Return True iff element is a qualifier element."""
        return element.localName in self.QUALIFIERS

    def qualifiers(self):
        """
        Return an iterable over the elements representing the qualifiers for
        this application.
        """
        quals = self.xml_element_children()
        return filter(self._is_qualifier, quals)

    def operands(self):
        """
        Return an iterable over the elements representing the operands for
        this application.
        """
        # Get all element children and strip the first (the operator)
        operands = self.xml_element_children()
        operands.next()
        # Now strip qualifiers from the front
        return itertools.dropwhile(self._is_qualifier, operands)

    def _get_operand_units(self):
        """
        Return an iterable containing a <units> element for each operand
        of this expression.
        """
        for o in self.operands():
            yield self._get_element_units(o)

    def evaluate(self):
        """
        Evaluate this expression, and return its value, if possible.
        """
        # Result depends on the operator
        op = self.operator()
        if hasattr(op, 'evaluate') and callable(op.evaluate):
            return op.evaluate()
        else:
            raise EvaluationError("Don't know how to evaluate the operator " + op.localName)

    def tree_complexity(self, **kw):
        """
        Calculate a rough estimate of the computation time for
        evaluating this <apply> element.

        Operates recursively, so the complexity of a function call is
        given by summing the complexity of the arguments and the time
        for evaluating the function itself.
        
        If lookup_tables is True, then assume we're using lookup tables
        where possible.
        If algebraic is True, the complexity is calculated as a dictionary,
        mapping node types to the number of occurences of that type.
        """
        kw['algebraic'] = kw.get('algebraic', False)
        if kw['algebraic']: ac = {}
        else: ac = 0

        # Complexity of this function
        op_name = self.operator().localName
        OPS = self.OPS
        if op_name in OPS.plusMinus or op_name in OPS.logical or op_name in OPS.relations:
            if kw['algebraic']: ac['op'] = (len(list(self.operands())) - 1)
            else: ac += 1 * (len(list(self.operands())) - 1)
        elif op_name in OPS.absRound:
            if op_name == 'abs':
                if kw['algebraic']: ac['abs'] = 1
                else: ac += 5
            else:
                if kw['algebraic']: ac['round'] = 1
                else: ac += 20
        elif op_name in OPS.elementary:
            if kw['algebraic']: ac['elementary'] = 1
            else: ac += 70
        elif op_name == 'times':
            if kw['algebraic']: ac['times'] = (len(list(self.operands())) - 1)
            else: ac += 1 * (len(list(self.operands())) - 1)
        elif op_name == 'divide':
            if kw['algebraic']: ac['divide'] = 1
            else: ac += 15
        elif op_name == 'power':
            # This will vary depending on the exponent - gcc can optimise for 2 and 3, it seems.
            exponent = list(self.operands())[1]
            if exponent.localName == u'cn' and unicode(exponent).strip() in [u'2', u'3']:
                if kw['algebraic']: ac['power2'] = 1
                else: ac += 5
            else:
                if kw['algebraic']: ac['power'] = 1
                else: ac += 30
        elif op_name == 'root':
            if kw['algebraic']: ac['root'] = 1
            else: ac += 30
        elif op_name == 'diff':
            if kw['algebraic']: ac['variable'] = 1
            else: ac += 0.7
        else:
            raise EvaluationError("Don't know complexity of operator " + op_name)

        # Complexity of operands
        for elt in self.operands():
            if kw['algebraic']:
                add_dicts(ac, self._tree_complexity(elt, **kw))
            else: ac += self._tree_complexity(elt, **kw)
        return ac

    def _set_in_units(self, units, no_act=False):
        """Set the units this expression should be given in.

        If these aren't our natural units (as given by an initial
        get_units) then we need to add units conversion code.
        """
        # First check we have something to do
        current_units = self.get_units()
        if units is current_units:
            return
        # Next, check if the required units can be achieved by suitable choices for operand units
        done = False
        if units in current_units:
            # They can!
            if not no_act:
                self._cml_units = units
            for src_units_set, src_units in current_units._get_sources(units):
                expr = src_units_set.get_expression()
                self._set_element_in_units(expr, src_units, no_act)
                done = True
            if not done and not no_act:
                # Some operators need this to be a UnitsSet
                self._cml_units = UnitsSet([units], self)
        if not done:
            # The behaviour now depends on the operator
            op = self.operator()
            if hasattr(op, '_set_in_units') and callable(op._set_in_units):
                op._set_in_units(units, no_act)
            else:
                raise UnitsError(self, u' '.join([
                    "Don't know how to select units for operands of operator",
                    op.localName, "when its units are", units.description()]))

    def get_units(self):
        """Recursively check this expression for dimensional consistency.

        Checks that the operands have suitable units.
        What constitutes 'suitable' depends on the operator; see appendix
        C.3.2 of the CellML 1.0 spec.

        If yes, returns a <units> element for the whole expression, based
        on the rules in appendix C.3.3.

        Throws a UnitsError if the units are inconsistent.
        """
        if self._cml_units:
            return self._cml_units
        our_units = None
        op = self.operator().localName
        operand_units = self._get_operand_units()
        # Tuples where second item is an index
        operand_units_idx = itertools.izip(operand_units, itertools.count(1))
        # Standard units objects
        dimensionless = self.model.get_units_by_name(u'dimensionless')
        boolean = self.model.get_units_by_name(u'cellml:boolean')

        if op in self.OPS.relations | self.OPS.plusMinus:
            our_units = operand_units.next().copy()
            # Operands mustn't be booleans
            if boolean in our_units:
                raise UnitsError(self, u' '.join([
                    u'Operator',op,u'has boolean operands, which does not make sense.']))
            # Operand units must be 'equivalent' (perhaps dimensionally)
            for u in operand_units:
                if not hasattr(self.model, '_cml_special_units_converter') and not our_units.dimensionally_equivalent(u):
                    raise UnitsError(self, u' '.join([
                        u'Operator',op,u'requires its operands to have',
                        u'dimensionally equivalent units;',u.description(),
                        u'and',our_units.description(),u'differ']))
                our_units.update(u)
            if op in self.OPS.relations:
                # Result has cellml:boolean units
                our_units = UnitsSet([boolean])
        elif op in self.OPS.logical:
            # Operand units must be cellml:boolean
            for u, i in operand_units_idx:
                if not boolean in u:
                    raise UnitsError(self, u' '.join([
                        u'Operator',op,u'requires operands to be booleans;',
                        u'operand',str(i),u'has units',u.description()]))
            # Result has cellml:boolean units
            our_units = UnitsSet([boolean])
        elif op in self.OPS.elementary:
            # Operand units must be dimensionless
            for u, i in operand_units_idx:
                if not u.dimensionally_equivalent(dimensionless):
                    raise UnitsError(self, u' '.join([
                        u'Operator',op,u'requires operands to be dimensionless;',
                        u'operand',str(i),u'has units',u.description()]))
            if op == 'log':
                # <logbase> qualifier must have units dimensionless
                if hasattr(self, u'logbase'):
                    base = _child1(self.logbase)
                    u = self._get_element_units(base)
                    if not u.dimensionally_equivalent(dimensionless):
                        raise UnitsError(self, u' '.join([u'The logbase qualifier must have dimensionless',
                                                          u'units, not',u.description()]))
            # Result has units of dimensionless
            our_units = UnitsSet([dimensionless])
        elif op == 'power':
            # Arg1 : any, Arg2 : dimensionless
            arg_units = operand_units.next()
            if boolean in arg_units:
                raise UnitsError(self, u'The argument of <power> should not be boolean')
            exponent_units = operand_units.next()
            if not exponent_units.dimensionally_equivalent(dimensionless):
                raise UnitsError(self, u' '.join([u'The second operand to power must have dimensionless',
                                                  u'units, not',exponent_units.description()]))
            # Result has units that are the units on the (first)
            # operand raised to the power of the second operand.  If
            # units on the first operand are dimensionless, then so is
            # the result.
            # TODO: Check how we could allow equiv. to d'less, instead of equal.
            # Need to consider any multiplicative factor...
            if arg_units.equals(dimensionless):
                our_units = UnitsSet([dimensionless])
            else:
                opers = self.operands()
                opers.next()
                # Make sure exponent is static
                expt = opers.next()
                if self._get_element_binding_time(expt) != BINDING_TIMES.static:
                    raise UnitsError(self, 'Unable to units check power with an exponent that can vary at run-time',
                                     warn=True,
                                     level=logging.WARNING_TRANSLATE_ERROR)
                # Try to evaluate the exponent
                try:
                    expt = self.eval(expt)
                except EvaluationError, e:
                    raise UnitsError(self, u' '.join([u'Unable to evaluate the exponent of a power element:', unicode(e)]),
                                     warn=True,
                                     level=logging.WARNING_TRANSLATE_ERROR)
                our_units = dimensionless.simplify(arg_units, expt)
        elif op == 'root':
            # Arg : any, <degree> : dimensionless
            arg_units = operand_units.next()
            if boolean in arg_units:
                raise UnitsError(self, u'The argument of <root> should not be boolean')
            if hasattr(self, u'degree'):
                degree = _child1(self.degree)
                u = self._get_element_units(degree)
                if not u.dimensionally_equivalent(dimensionless):
                    raise UnitsError(self, u' '.join([
                        u'The degree qualifier must have dimensionless units, not',u.description()]))
            else:
                degree = 2.0 # Default is square root
            # Result has units that are the units on the (first) operand
            # raised to the power of the reciprocal of the value of the
            # degree qualifier.
            # TODO: If units on the first operand are dimensionless,
            # then so is the result.
            if not type(degree) is float:
                # Make sure degree is static
                if self._get_element_binding_time(degree) != BINDING_TIMES.static:
                    raise UnitsError(self, 'Unable to units check root with a degree that can vary at run-time',
                                     warn=True,
                                     level=logging.WARNING_TRANSLATE_ERROR)
                try:
                    degree = self.eval(degree)
                except EvaluationError, e:
                    raise UnitsError(self, u' '.join([u'Unable to evaluate the degree of a root element:', unicode(e)]),
                                     warn=True,
                                     level=logging.WARNING_TRANSLATE_ERROR)
            our_units = dimensionless.simplify(arg_units, 1/degree)
        elif op == 'diff':
            # Arg : any, <bvar> : any, <degree> : dimensionless
            arg_units = operand_units.next()
            if boolean in arg_units:
                raise UnitsError(self, u'The argument of <diff> should not be boolean')
            if hasattr(self, u'bvar'):
                if hasattr(self.bvar, u'degree'):
                    degree = _child1(self.bvar.degree)
                    u = self._get_element_units(degree)
                    if not u.dimensionally_equivalent(dimensionless):
                        raise UnitsError(self, u' '.join([
                            u'The degree qualifier must have dimensionless units, not',u.description()]))
                else:
                    degree = 1.0 # Default is first derivative
            else:
                raise UnitsError(self, u'A diff operator must have a bvar qualifier')
            # Result has units that are the quotient of the units of the
            # operand, over the units of the term in the bvar qualifier
            # raised to the value of the degree qualifier
            if not type(degree) is float:
                # Make sure exponent is static
                if self._get_element_binding_time(degree) != BINDING_TIMES.static:
                    raise UnitsError(self, 'Unable to units check derivative with a degree that can vary at run-time',
                                     warn=True,
                                     level=logging.WARNING_TRANSLATE_ERROR)
                try:
                    degree = self.eval(degree)
                except EvaluationError, e:
                    raise UnitsError(self, u' '.join([u'Unable to evaluate the degree of a diff element:', unicode(e)]),
                                     warn=True,
                                     level=logging.WARNING_TRANSLATE_ERROR)
            for e in self.xml_element_children(self.bvar):
                if not e.localName == u'degree':
                    bvar_units = self._get_element_units(e)
                    break
            else:
                raise UnitsError(self, u'diff element does not have a valid bvar')
            our_units = arg_units.simplify(bvar_units, -degree)
        elif op in self.OPS.absRound | self.OPS.timesDivide:
            # No restrictions on operand units, except that they shouldn't be boolean
            for u in self._get_operand_units():
                if boolean in u:
                    raise UnitsError(self, u' '.join([
                        u'Operator',op,u'has boolean operands, which does not make sense.']))
            if op == 'times':
                # Result has units that are the product of the operand units
                our_units = operand_units.next().copy()
                for u in operand_units:
                    our_units = our_units.simplify(other_units=u)
            elif op == 'divide':
                # Result has units that are the quotient of the units
                # on the first and second operands
                our_units = operand_units.next()
                our_units = our_units.simplify(other_units=operand_units.next(), other_exponent=-1)
            else:
                # Result has same units as operands
                our_units = operand_units.next().copy()
        else:
            # Warning: unsupported operator!
            raise UnitsError(self, u' '.join([u'Unsupported operator for units checking:', op]),
                             warn=True,
                             level=logging.WARNING_TRANSLATE_ERROR)

        # Cache & return result
        if isinstance(our_units, cellml_units):
            # Units conversion has been done, then PE
            our_units = UnitsSet([our_units])
        self._cml_units = our_units
        our_units.set_expression(self)
        return self._cml_units

    def check_assigned_var(self):
        """Check the current component owns the variable being assigned to.
        
        Should only be called if this object represents an assignment
        expression.  Checks that the variable being assigned to doesn't
        have an interface value of 'in'.  If this isn't a simple assignment
        (i.e. the LHS isn't a plain ci element) then the check succeeds
        automatically.
        
        Adds to the model's error list if the check fails.  This method
        always returns None.
        """
        # Check this is an application of eq
        if self.operator().localName != u'eq':
            raise MathsError(self, u'Top-level mathematics expressions should be assigment expressions.')
        first_operand = self.operands().next()
        if first_operand.localName == u'ci':
            # We've already checked that the variable exists
            var = first_operand.variable
            for iface in [u'public_interface', u'private_interface']:
                if getattr(var, iface, u'none') == u'in':
                    raise MathsError(self, u' '.join([
                        u'Variable', var.fullname(),
                        u'is assigned to in a math element, but has its',
                        iface, u'set to "in".']))

    def classify_variables(self, root=False,
                           dependencies_only=False,
                           needs_special_treatment=lambda n: None):
        """
        Classify variables in this expression according to how they are
        used.
        
        In the process, compute and return a set of variables on which
        this expression depends.  If root is True, store this set as a
        list, to represent edges of a dependency graph.
        Also, if root is True then this node is the root of an expression
        (so it will be an application of eq); treat the LHS differently.
        
        If dependencies_only then the variable classification will not be
        done, only dependencies will be analysed.  This is useful for doing
        a 'light' re-analysis if the dependency set has been reduced; if the
        set has increased then the topological sort of equations may need to
        be redone.
        
        The function needs_special_treatment may be supplied to override the
        default recursion into sub-trees.  It takes a single sub-tree as
        argument, and should either return the dependency set for that
        sub-tree, or None to use the default recursion.  This is used when
        re-analysing dependencies after applying lookup tables, since table
        lookups only depend on the keying variable.
        """
        dependencies = set()
        ode_indep_var = None
        op = self.operator()
        if op.localName == u'diff':
            # This is a derivative dy/dx on the RHS of an assignment.
            # Store the dependency as a pair (y,x)
            dependencies.add((op.dependent_variable, op.independent_variable))
            if not dependencies_only:
                # Set variable types
                op._set_var_types()
        else:
            opers = self.operands()
            if root:
                # Treat the LHS of the assignment
                lhs = opers.next()
                if lhs.localName == u'ci':
                    # Direct assignment to variable
                    var = lhs.variable
                    var._add_dependency(self)
                    self._cml_assigns_to = var
                    if not dependencies_only:
                        # Check for possibly conflicting types
                        t = var.get_type()
                        if t == VarTypes.Constant or t == VarTypes.MaybeConstant:
                            self.model.validation_warning(
                                u' '.join([
                                u'Variable',var.fullname(),u'is assigned to',
                                u'and has an initial value set.']),
                                level=logging.WARNING_TRANSLATE_ERROR)
                        elif t == VarTypes.State or t == VarTypes.Free:
                            self.model.validation_warning(
                                u' '.join([
                                u'Variable',var.fullname(),u'is assigned to',
                                u'and appears on the LHS of an ODE.']),
                                level=logging.WARNING_TRANSLATE_ERROR)
                        var._set_type(VarTypes.Computed)
                elif lhs.localName == u'apply':
                    # This could be an ODE
                    diff = lhs.operator()
                    if diff.localName == u'diff':
                        # It is an ODE. TODO: Record it somewhere?
                        if not dependencies_only:
                            diff._set_var_types()
                        dep = diff.dependent_variable
                        indep = diff.independent_variable
                        dep._add_ode_dependency(indep, self)
                        # An ODE should depend on its independent variable
                        ode_indep_var = indep
                        if not dependencies_only:
                            indep._used()
                        # TODO: Hack; may remove.
                        self._cml_assigns_to = (dep, indep)
                    else:
                        raise MathsError(self, u'Assignment statements are expected to be an ODE or assign to a variable.',
                                         warn=True,
                                         level=logging.WARNING_TRANSLATE_ERROR)
                else:
                    raise MathsError(self, u'Assignment statements are expected to be an ODE or assign to a variable.',
                                     warn=True,
                                     level=logging.WARNING_TRANSLATE_ERROR)

            # Consider operands other than the LHS of an assignment
            for oper in opers:
                dependencies.update(self.classify_child_variables(oper, dependencies_only=dependencies_only,
                                                                  needs_special_treatment=needs_special_treatment))

        if ode_indep_var:
            # ODEs should depend on their independent variable.
            # However, for code generation we wish to distinguish
            # whether the independent variable appears on the RHS or
            # not.
            if ode_indep_var in dependencies:
                self._cml_ode_has_free_var_on_rhs = True
            else:
                self._cml_ode_has_free_var_on_rhs = False
                dependencies.add(ode_indep_var)

        if root:
            # Store dependencies
            self._cml_depends_on = list(dependencies)
        return dependencies

    def is_top_level(self):
        """Test whether this is a top-level assignment expression."""
        return self._cml_assigns_to is not None

    def is_ode(self):
        """Return True iff this is the assignment of an ODE.

        Only makes sense if called on a top-level assignment
        expression, and checks if it represents an ODE, i.e. if the
        LHS is a derivative.
        """
        if not self.is_top_level():
            return False
        return type(self._cml_assigns_to) == types.TupleType

    def is_assignment(self):
        """Return True iff this is a straightforward assignment expression.

        Only makes sense if called on a top-level assignment expression.
        Checks that this is *not* an ODE, but assigns to a single variable.
        """
        if not self.is_top_level():
            return False
        return isinstance(self._cml_assigns_to, cellml_variable)

    def assigned_variable(self):
        """Return the variable assigned to by this assignment.

        Should only be called on a top-level assignment expression.
        
        If it's a straightforward assignment (so self.is_assignment()
        returns True) then return the cellml_variable object
        representing the variable assigned to.

        If it's an ODE, return a pair
        (dependent variable, independent variable).
        """
        if not self.is_top_level():
            raise TypeError("not a top-level apply element")
        else:
            return self._cml_assigns_to

    def _get_binding_time(self, check_operator=True):
        """Return the binding time of this expression.

        The binding time will be computed recursively and cached.
        It will also be made available as an attribute in the XML.

        It is computed by taking the least upper bound of the binding
        times of our operands, unless the operator possesses an
        alternative method.
        """
        if self._cml_binding_time is not None:
            return self._cml_binding_time

        # Do we have an annotation?
        if hasattr(self, u'binding_time'):
            self._cml_binding_time = getattr(BINDING_TIMES, self.binding_time)
            return self._cml_binding_time

        # Does operator have a specialised method for this?
        op = self.operator()
        if check_operator and hasattr(op, '_get_binding_time'):
            self._cml_binding_time = op._get_binding_time()
        else:
            # Compute operand binding times
            bts = [BINDING_TIMES.static]
            for operand in self.operands():
                bts.append(self._get_element_binding_time(operand))
            # Take l.u.b.
            self._cml_binding_time = max(bts)

        # Annotate the element with the binding time
        self.xml_set_attribute((u'pe:binding_time', NSS[u'pe']), unicode(self._cml_binding_time))
        return self._cml_binding_time

    def _reduce(self, check_operator=True):
        """Reduce this expression by evaluating its static parts."""
        # Check to see if this operator requires a special
        # reduction strategy
        op = self.operator()
        DEBUG('partial-evaluator', "Reducing", op.localName, getattr(self, u'id', ''))
        if check_operator and hasattr(op, '_reduce'):
            op._reduce()
        else:
            bt = self._get_binding_time()
            if bt == BINDING_TIMES.static:
                # Evaluate self and replace by a <cn>, <true> or <false>
                new_elt = self._eval_self()
                self._xfer_complexity(new_elt)
                self.replace_child(self, new_elt, self.xml_parent)
                # Update usage counts
                self._update_usage_counts(self, remove=True)
            elif bt == BINDING_TIMES.dynamic:
                # Recurse into operands and reduce those
                for op in self.operands():
                    self._reduce_elt(op)
        return

    @staticmethod
    def create_new(elt, operator, operands=[], qualifiers=[]):
        """Create a new MathML apply element, with given content.

        elt should be any element in the document.

        operator is used as the name of the first, empty, child.
        
        operands is a list, possibly empty, of operand elements.  If
        any member is a unicode object, it is considered to be the
        name of a variable.  If a tuple, then it should be a pair of
        unicode objects: (number, units).  (Although units can be an
        attribute dictionary.)

        qualifiers specifies a list of qualifier elements.
        """
        app = elt.xml_create_element(u'apply', NSS[u'm'])
        app.xml_append(app.xml_create_element(unicode(operator), NSS[u'm']))
        for qual in qualifiers:
            check_append_safety(qual)
            app.xml_append(qual)
        for op in operands:
            if isinstance(op, basestring):
                # Variable name
                op = app.xml_create_element(u'ci', NSS[u'm'],
                                            content=unicode(op))
            elif isinstance(op, tuple):
                # Constant with units
                if isinstance(op[1], dict):
                    attrs = op[1]
                else:
                    attrs = {(u'cml:units', NSS[u'cml']): unicode(op[1])}
                op = app.xml_create_element(u'cn', NSS[u'm'],
                                            attributes=attrs,
                                            content=unicode(op[0]))
            else:
                # Should already be an element
                check_append_safety(op)
            app.xml_append(op)
        return app

class mathml_piecewise(mathml_constructor, mathml_units_mixin):
    def __init__(self):
        super(mathml_piecewise, self).__init__()
        self._cml_units = None
        self._cml_binding_time = None

    def tree_complexity(self, **kw):
        """
        Calculate a rough estimate of the computation time for
        evaluating this <piecewise> element.

        The real evaluation time will generally depend on run time
        data, which makes things tricky.  Here we estimate by taking
        the sum of the complexities of the conditions and the maximum
        of the complexity of the cases, in order to obtain an upper
        bound.
        
        If lookup_tables is True, then assume we're using lookup tables
        where possible.
        If algebraic is True, the complexity is calculated as a dictionary,
        mapping node types to the number of occurences of that type.

        If self.rootNode.num_lookup_tables exists, this method will
        update the count of lookup tables based on this expression,
        unless the argument 'count_tables' is False or algebraic is True.
        """
        kw['algebraic'] = kw.get('algebraic', False)
        alg = kw['algebraic']
        if alg: ac, piece_dicts = {}, []
        else: ac = 0
        piece_acs = []
        count_lts = hasattr(self.rootNode, 'num_lookup_tables') and kw.get('count_tables', True) and not alg
        if count_lts:
            # Alternative method of counting number of lookup tables;
            # handles the Zhang model better!
            num_lts = self.rootNode.num_lookup_tables
            piece_num_lts = []

        for piece in getattr(self, u'piece', []):
            test_ac = self._tree_complexity(child_i(piece, 2), **kw)
            if alg: add_dicts(ac, test_ac)
            else: ac += test_ac
            if count_lts:
                nlts = self.rootNode.num_lookup_tables
            piece_ac = self._tree_complexity(child_i(piece, 1), **kw)
            if alg:
                piece_dicts.append(piece_ac)
                piece_acs.append(self._tree_complexity(child_i(piece, 1), count_tables=False))
            else:
                piece_acs.append(piece_ac)
            if count_lts:
                piece_num_lts.append(self.rootNode.num_lookup_tables - nlts)
        if hasattr(self, u'otherwise'):
            if count_lts:
                nlts = self.rootNode.num_lookup_tables
            ow_ac = self._tree_complexity(child_i(self.otherwise, 1), **kw)
            if alg:
                piece_dicts.append(ow_ac)
                piece_acs.append(
                    self._tree_complexity(child_i(self.otherwise, 1), count_tables=False))
            else:
                piece_acs.append(ow_ac)
            if count_lts:
                piece_num_lts.append(self.rootNode.num_lookup_tables - nlts)
        max_idx, max_piece_ac = max_i(piece_acs)
        if alg:
            add_dicts(ac, piece_dicts[max_idx])
        else:
            ac += max_piece_ac
        if count_lts:
            self.rootNode.num_lookup_tables -= sum(piece_num_lts)
            self.rootNode.num_lookup_tables += piece_num_lts[max_idx]
        return ac

    def _set_in_units(self, units, no_act=False):
        """Set the units this expression should be given in.

        This is done recursively by setting the units for each option.
        
        We also set the units on each condition to be boolean, since
        subexpressions of the conditions may need units conversions added.
        """
        # First, record our units
        if not no_act:
            self._cml_units = units
        # Now process our children
        boolean = self.model.get_units_by_name(u'cellml:boolean')
        for piece in getattr(self, u'piece', []):
            self._set_element_in_units(child_i(piece, 1), units, no_act)
            self._set_element_in_units(child_i(piece, 2), boolean, no_act)
        if hasattr(self, u'otherwise'):
            self._set_element_in_units(child_i(self.otherwise, 1), units, no_act)

    def get_units(self):
        """Recursively check this expression for dimensional consistency.

        The first child elements of each <piece> and <otherwise> element
        should have dimensionally equivalent units (the resulting <units>
        element will be dimensionally equivalent to these).  The second child
        elements of each <piece> should have units of cellml:boolean.

        If consistent, returns a <units> element for the whole expression.
        Throws a UnitsError if the units are inconsistent.
        """
        if self._cml_units:
            return self._cml_units
        # Check the second child of each <piece> element
        boolean = self.model.get_units_by_name(u'cellml:boolean')
        for piece in getattr(self, u'piece', []):
            cond_elt = child_i(piece, 2)
            units = self._get_element_units(cond_elt)
            if not boolean in units:
                raise UnitsError(self, u' '.join([
                    u'The second child element of a <piece> element must have'
                    u'units of cellml:boolean, not',units.description()]))
        # Compare the first child element of each <piece> and the <otherwise>,
        # if present.
        our_units = None
        if hasattr(self, u'otherwise'):
            value_elt = child_i(self.otherwise, 1)
            our_units = self._get_element_units(value_elt).copy()
        for piece in getattr(self, u'piece', []):
            value_elt = child_i(piece, 1)
            if our_units is None:
                our_units = self._get_element_units(value_elt).copy()
            else:
                units = self._get_element_units(value_elt)
                if not our_units.dimensionally_equivalent(units):
                    raise UnitsError(self, u' '.join([
                        u'The first child elements of children of a piecewise',
                        u'element must have dimensionally equivalent units;',
                        units.description(),'and',our_units.description(),
                        u'differ']))
                our_units.update(units)
        # Check that we have some units for this element
        if our_units is None:
            raise UnitsError(self, u' '.join([
                u'A piecewise element must have at least one piece or',
                u'otherwise child in order to have defined units.']))
        # Cache & return units
        self._cml_units = our_units
        our_units.set_expression(self)
        return self._cml_units

    def classify_variables(self, dependencies_only=False,
                           needs_special_treatment=lambda n: None):
        """Classify variables in this expression according to how they are used.
        
        In the process, compute and return a set of variables on which
        this expression depends.

        If dependencies_only then the variable classification will not be
        done, only dependencies will be analysed.  This is useful for doing
        a 'light' re-analysis if the dependency set has been reduced; if the
        set has increased then the topological sort of equations may need to
        be redone.
        
        The function needs_special_treatment may be supplied to override the
        default recursion into sub-trees.  It takes a single sub-tree as
        argument, and should either return the dependency set for that
        sub-tree, or None to use the default recursion.  This is used when
        re-analysing dependencies after applying lookup tables, since table
        lookups only depend on the keying variable.
        """
        dependencies = set()
        pieces = list(getattr(self, u'piece', []))
        if hasattr(self, u'otherwise'):
            pieces.append(self.otherwise)
        for piece in pieces:
            dependencies.update(self.classify_child_variables(piece, dependencies_only=dependencies_only,
                                                              needs_special_treatment=needs_special_treatment))
        return dependencies

    def evaluate(self):
        """Evaluate this piecewise expression.

        Tests choices in the order they occur in the file.
        Only evaluates a choice if its condition evaluates to True.
        """
        for piece in getattr(self, u'piece', []):
            condition = child_i(piece, 2)
            cond_value = self.eval(condition)
            if cond_value is True:
                # This is the option to take
                value = self.eval(child_i(piece, 1))
                break
        else:
            # Evaluate the <otherwise>
            if hasattr(self, u'otherwise'):
                value = self.eval(_child1(self.otherwise))
            else:
                raise EvaluationError(u' '.join([
                    "A piecewise element where the pieces aren't mutually",
                    "exhaustive requires an otherwise element."]))
        return value

    def _get_binding_time(self):
        """Return the binding time of this expression.

        The binding time will be computed recursively and cached.
        It will also be made available as an attribute in the XML.

        It is computed by taking the least upper bound of the binding
        times of (some of) the conditions and cases.
        
        Condition & case binding times are computed in the order given
        in the file.  If a condition is static with value False, its
        associated case is not considered.  If a condition is static
        with value True, subsequent conditions & cases and the
        otherwise (if present) will not be considered.
        """
        if self._cml_binding_time is not None:
            return self._cml_binding_time

        # Do we have an annotation?
        if hasattr(self, u'binding_time'):
            self._cml_binding_time = getattr(BINDING_TIMES,
                                             self.binding_time)
            return self._cml_binding_time

        # Compute condition binding times
        bts = [BINDING_TIMES.static]
        for piece in getattr(self, u'piece', []):
            condition = child_i(piece, 2)
            bt = self._get_element_binding_time(condition)
            if bt is BINDING_TIMES.static:
                cond_value = self.eval(condition)
                if cond_value is True:
                    # Compute BT for associated case
                    bts.append(self._get_element_binding_time(
                        child_i(piece, 1)))
                    # Skip remaining conditions & otherwise
                    break
            else:
                # Don't need to append extra statics, since bts
                # contains at least one member that is static
                bts.append(bt)
                # Compute BT for associated case
                bts.append(self._get_element_binding_time(
                    child_i(piece, 1)))
        else:
            # Consider the <otherwise> element
            if hasattr(self, u'otherwise'):
                bts.append(self._get_element_binding_time(
                    child_i(self.otherwise, 1)))
        # Take least upper bound of appropriate binding times
        self._cml_binding_time = max(bts)

        # Annotate the element with the binding time
        self.xml_set_attribute((u'pe:binding_time', NSS[u'pe']),
                               unicode(self._cml_binding_time))
        return self._cml_binding_time

    def _reduce(self, check_operator=True):
        """Reduce this expression by evaluating its static parts.

        Even in a dynamic conditional, where a condition is static and
        evaluates to False, the associated case is discarded.
        """
        # Element to replace this <piecewise> with, if any
        new_elt = None
        if self._get_binding_time() == BINDING_TIMES.static:
            # Evaluate self and replace by a <cn>, <true> or <false>
            new_elt = self._eval_self()
        elif self._get_binding_time() == BINDING_TIMES.dynamic:
            # Go through pieces and reduce where appropriate
            deletable_pieces = []
            found_dynamic_piece = False
            for piece in getattr(self, u'piece', []):
                condition = child_i(piece, 2)
                bt = self._get_element_binding_time(condition)
                if bt is BINDING_TIMES.static:
                    cond_value = self.eval(condition)
                    if cond_value is True:
                        if not found_dynamic_piece:
                            # Replace the entire piecewise element by our case
                            # We don't replace if a previous piece had a
                            # dynamic condition, since this would change
                            # execution order, which could alter the semantics.
                            new_elt = child_i(piece, 1)
                            break
                    else:
                        # Discard this condition & case
                        deletable_pieces.append(piece)
                elif bt is BINDING_TIMES.dynamic:
                    found_dynamic_piece = True
                    # Reduce the condition & case
                    self._reduce_elt(condition)
                    self._reduce_elt(child_i(piece, 1))
            else:
                # Didn't replace entire conditional
                # Remove pieces with False conditions
                for piece in deletable_pieces:
                    self._update_usage_counts(piece, remove=True)
                    self.xml_remove_child(piece)
                # Consider the <otherwise> element
                if hasattr(self, u'otherwise'):
                    if not found_dynamic_piece:
                        # All the <piece> elements were removed, so replace
                        # the entire conditional by this <otherwise>
                        new_elt = child_i(self.otherwise, 1)
                    else:
                        # Just reduce the <otherwise>
                        self._reduce_elt(child_i(self.otherwise, 1))
        # Replace this element, if required
        if new_elt is not None:
            # Update usage counts for removed expressions
            for piece in getattr(self, u'piece', []):
                if not new_elt is child_i(piece, 1):
                    self._update_usage_counts(piece, remove=True)
                else:
                    # Condition is being removed
                    self._update_usage_counts(child_i(piece, 2), remove=True)
                    piece.xml_remove_child(child_i(piece, 2))
                    piece.xml_remove_child(new_elt)
            if hasattr(self, u'otherwise') and \
               not new_elt is child_i(self.otherwise, 1):
                self._update_usage_counts(child_i(self.otherwise, 1),
                                          remove=True)
            # Do the replace
            self._xfer_complexity(new_elt)
            self.replace_child(self, new_elt, self.xml_parent)
            # May need to reduce our replacement
            self._reduce_elt(new_elt)

    @staticmethod
    def create_new(elt, pieces, otherwise=None):
        """Create a new piecewise element.

        elt is any element in the current document.

        pieces is a list of pairs of expressions: (case, condition).

        otherwise, if given, is the default case.
        """
        piecewise = elt.xml_create_element(u'piecewise', NSS[u'm'])
        for piece in pieces:
            case, cond = piece
            check_append_safety(case)
            check_append_safety(cond)
            piece_elt = elt.xml_create_element(u'piece', NSS[u'm'])
            piece_elt.xml_append(case)
            piece_elt.xml_append(cond)
            piecewise.xml_append(piece_elt)
        if otherwise:
            check_append_safety(otherwise)
            otherwise_elt = elt.xml_create_element(u'otherwise', NSS[u'm'])
            otherwise_elt.xml_append(otherwise)
            piecewise.xml_append(otherwise_elt)
        return piecewise


class mathml_lambda(mathml_constructor):
    """Class representing the MathML lambda construct.
    
    Note that we don't support lambda occuring in CellML models.  However, it is used
    for defining special units conversion rules using the protocol syntax.
    """
    @staticmethod
    def create_new(elt, bound_var_names, body_expr):
        """Create a new lambda from the sequence of bvar names and expression."""
        lambda_ = elt.xml_create_element(u'lambda', NSS[u'm'])
        for bvar_name in bound_var_names:
            bvar = elt.xml_create_element(u'bvar', NSS[u'm'])
            bvar.xml_append(mathml_ci.create_new(elt, bvar_name))
            lambda_.xml_append(bvar)
        check_append_safety(body_expr)
        lambda_.xml_append(body_expr)
        return lambda_


class mathml_operator(mathml):
    """Base class for MathML operator elements."""
    def wrong_number_of_operands(self, found, wanted):
        """Raise an EvaluationError due to wrong operand count.

        found is the number of operands found; wanted is a list of suitable
        numbers of operands.
        """
        raise EvaluationError(u''.join([
            "Wrong number of operands for <", self.localName, "> ",
            "(found ", str(found), "; wanted", ' or '.join(map(str, wanted)), ")"]))


class mathml_diff(mathml_operator):
    """Class representing the diff element, containing some useful methods."""
    @property
    def independent_variable(self):
        """
        Return the variable object w.r.t which we are differentiating.
        """
        # Note that after units checking the <bvar> qualifier can be
        # assumed to exist.
        apply_elt = self.xml_parent
        if not hasattr(apply_elt.bvar, u'ci'):
            raise MathsError(apply_elt, u'Differential operator does not have a variable element as bound variable.')
        return apply_elt.bvar.ci.variable

    @property
    def dependent_variable(self):
        """
        Return the variable object being differentiated.
        """
        apply_elt = self.xml_parent
        operand = apply_elt.operands().next()
        if not operand.localName == u'ci':
            raise MathsError(apply_elt, u'Derivatives of non-variables are not supported.')
        return operand.variable

    def _set_var_types(self):
        """
        Set the types of the dependent & independent variables: State for
        the dependent variable and Free for the independent variable.
        Gives a validation warning if they already have 'incompatible'
        types.
        """
        dep, indep = self.dependent_variable, self.independent_variable
        model = self.xml_parent.model
        
        # The dependent variable should have an initial value
        if not dep.get_type() == VarTypes.Mapped and \
               not hasattr(dep, u'initial_value'):
            model.validation_warning(u' '.join([
                u'The state variable',dep.fullname(),
                u'does not have an initial value given.']),
                                     level=logging.WARNING_TRANSLATE_ERROR)
        # It doesn't make sense to compute a state variable
        if dep.get_type(follow_maps=True) == VarTypes.Computed:
            model.validation_warning(u' '.join([
                u'The state variable',dep.fullname(),
                u'is also assigned to directly.']),
                                     level=logging.WARNING_TRANSLATE_ERROR)
        dep._set_type(VarTypes.State)
        
        t = indep.get_type(follow_maps=True)
        if t != VarTypes.Free:
            if t != VarTypes.Unknown:
                if t == VarTypes.Computed:
                    reason = u'is computed in an expression.'
                elif t == VarTypes.State:
                    reason = u'is a state variable itself.'
                else:
                    reason = u'has an initial value specified.'
                model.validation_warning(u' '.join([
                    u'The derivative of',dep.fullname(),
                    u'is taken with respect to',indep.fullname(),
                    u'but the latter', reason]),
                                         level=logging.WARNING_TRANSLATE_ERROR)
            # TODO: Add to list of independent vars?
            indep._set_type(VarTypes.Free)

    def _get_binding_time(self):
        """Return the binding time of the enclosing <apply> element.

        This is the binding time of the expression defining this ODE.
        """
        expr = self.dependent_variable.get_ode_dependency(
            self.independent_variable)
        return expr._get_binding_time()

    def _reduce(self):
        """Reduce this expression by evaluating its static parts.

        If the whole expression is static, proceed as normal for an
        <apply>.  Otherwise just rename the variable references.  We
        can't instantiate the definition, because there will always be
        another user - external code.
        
        This operator is special cased because we need to alter its
        qualifier, but mathml_apply only considers operands.  MathML
        data binding can be annoying at times!
        """
        app = self.xml_parent
        bt = app._get_binding_time()
        if bt == BINDING_TIMES.static:
            # Evaluate this expression as normal
            app._reduce(check_operator=False)
        else:
            # Just update names to be canonical.
            for ci in [app.ci, app.bvar.ci]:
                ci._set_variable_obj(ci.variable.get_source_variable(recurse=True))
                ci._rename()

    @staticmethod
    def create_new(elt, bvar, state_var, rhs):
        """Construct an ODE expression: d(state_var)/d(bvar) = rhs."""
        bvar_elt = elt.xml_create_element(u'bvar', NSS[u'm'])
        bvar_elt.xml_append(mathml_ci.create_new(elt, bvar))
        diff = mathml_apply.create_new(elt, u'diff', [state_var], [bvar_elt])
        ode = mathml_apply.create_new(elt, u'eq', [diff, rhs])
        return ode

class reduce_commutative_nary(object):
    def _reduce(self):
        """Reduce this expression by evaluating its static parts.

        If the whole expression is static, proceed as normal for an
        <apply>.  Otherwise check if we have more than one static
        operand.  If we do, we can combine them in a new static
        expression and evaluate that as a whole.
        """
        app = self.xml_parent
        bt = app._get_binding_time()
        if bt == BINDING_TIMES.static or not self.model.get_option('partial_pe_commutative'):
            # Evaluate this expression as normal
            app._reduce(check_operator=False)
        else:
            # Get binding times of operands
            static_opers = filter(lambda e: app._get_element_binding_time(e) == BINDING_TIMES.static,
                                  app.operands())
            if len(static_opers) > 1:
                # Remove them from app
                for oper in static_opers:
                    app.safe_remove_child(oper)
                # Create the new expression
                new_expr = mathml_apply.create_new(self, self.localName, static_opers)
                # Put it in the model and reduce it to evaluate it properly
                app.xml_append(new_expr)
                new_expr._reduce()
            # Now reduce all our operands as normal
            app._reduce(check_operator=False)

class mathml_plus(mathml_operator, mathml_units_mixin_set_operands, reduce_commutative_nary):
    """Class representing the MathML <plus> operator."""
    def __init__(self):
        super(mathml_plus, self).__init__()
        return

    def evaluate(self):
        """Evaluate by summing the operands of the enclosing <apply>."""
        app = self.xml_parent
        ops = app.operands()
        value = 0
        for operand in ops:
            value += self.eval(operand)
        return value

class mathml_minus(mathml_operator, mathml_units_mixin_set_operands):
    """Class representing the MathML <minus> operator."""
    def evaluate(self):
        """Evaluate the enclosing <apply> element.

        Behaviour depends on the number of operands: we perform
        either a unary or binary minus.
        """
        app = self.xml_parent
        ops = list(app.operands())
        if len(ops) == 1:
            value = -self.eval(ops[0])
        elif len(ops) == 2:
            value = self.eval(ops[0]) - self.eval(ops[1])
        else:
            self.wrong_number_of_operands(len(ops), [1, 2])
        return value

class mathml_times(mathml_operator, mathml_units_mixin_choose_nearest, reduce_commutative_nary):
    """Class representing the MathML <times> operator."""
    def evaluate(self):
        """
        Evaluate by taking the product of the operands of the enclosing <apply>.
        """
        app = self.xml_parent
        ops = app.operands()
        value = 1
        for operand in ops:
            value *= self.eval(operand)
        return value

class mathml_divide(mathml_operator, mathml_units_mixin_choose_nearest):
    """Class representing the MathML <divide> operator."""
    def evaluate(self):
        """Evaluate by dividing the 2 operands of the enclosing <apply>."""
        app = self.xml_parent
        ops = list(app.operands())
        if len(ops) != 2:
            self.wrong_number_of_operands(len(ops), [2])
        numer = self.eval(ops[0])
        denom = self.eval(ops[1])
        return numer/denom

    def _reduce(self):
        """Reduce this expression by evaluating its static parts.

        If the whole expression is static, proceed as normal for an
        <apply>.  If just the denominator is static, transform the
        expression into a multiplication.
        """
        app = self.xml_parent
        bt = app._get_binding_time()
        if bt == BINDING_TIMES.static:
            # Evaluate this expression as normal
            app._reduce(check_operator=False)
        else:
            # Check binding time of the denominator
            ops = list(app.operands())
            if len(ops) != 2:
                self.wrong_number_of_operands(len(ops), [2])
            bt = app._get_element_binding_time(ops[1])
            if bt == BINDING_TIMES.static:
                # Create inverse expression and evaluate it
                dummy = self.xml_create_element(u'dummy', NSS[u'm'])
                app.replace_child(ops[1], dummy)
                new_expr = mathml_apply.create_new(
                    self, u'divide', [(u'1', u'dimensionless'),
                                      ops[1]])
                app.replace_child(dummy, new_expr)
                app._reduce_elt(new_expr)
                # Change this expression to a <times>
                times = self.xml_create_element(u'times', NSS[u'm'])
                app.replace_child(self, times)
                # And finally reduce it as normal
                app._reduce(check_operator=False)
            else:
                # Evaluate this expression as normal
                app._reduce(check_operator=False)

class mathml_exp(mathml_operator, mathml_units_mixin_set_operands):
    """Class representing the MathML <exp> operator."""
    def evaluate(self):
        """Return e to the power of the single operand."""
        app = self.xml_parent
        ops = list(app.operands())
        if len(ops) != 1:
            self.wrong_number_of_operands(len(ops), [1])
        return math.exp(self.eval(ops[0]))

class mathml_ln(mathml_operator, mathml_units_mixin_set_operands):
    """Class representing the MathML <ln> operator."""
    def evaluate(self):
        """Return the natural logarithm of the single operand."""
        app = self.xml_parent
        ops = list(app.operands())
        if len(ops) != 1:
            self.wrong_number_of_operands(len(ops), [1])
        return math.log(self.eval(ops[0]))

class mathml_log(mathml_operator, mathml_units_mixin_set_operands):
    """Class representing the MathML <log> operator."""
    def evaluate(self):
        """Return the logarithm of the single operand."""
        app = self.xml_parent
        ops = list(app.operands())
        if len(ops) != 1:
            self.wrong_number_of_operands(len(ops), [1])
        if hasattr(app, u'logbase'):
            base = self.eval(app.logbase)
        else:
            base = 10
        return math.log(self.eval(ops[0]), base)

class mathml_abs(mathml_operator, mathml_units_mixin_set_operands):
    """Class representing the MathML <abs> operator."""
    def evaluate(self):
        """Return the absolute value of the single operand."""
        app = self.xml_parent
        ops = list(app.operands())
        if len(ops) != 1:
            self.wrong_number_of_operands(len(ops), [1])
        return abs(self.eval(ops[0]))

class mathml_power(mathml_operator, mathml_units_mixin):
    """Class representing the MathML <power> operator."""
    def _set_in_units(self, units, no_act=False):
        """Set the units of the application of this operator.

        Set the exponent to have units of dimensionless, and the operand to
        have an arbitrary member of its possible units set.

        Where these mean the <apply> doesn't have the given units, wrap it
        in suitable units conversion mathematics.
        """
        app = self.xml_parent
        defn_units_set = app.get_units()
        defn_units = defn_units_set.extract()
        app._add_units_conversion(app, defn_units, units, no_act)
        # Record which member of the set we used
        if not no_act:
            app._cml_units = defn_units
        # Set exponent units
        dimensionless = app.model.get_units_by_name('dimensionless')
        ops = list(app.operands())
        self._set_element_in_units(ops[1], dimensionless, no_act)
        # Set operand units
        for src_units_set, src_units in defn_units_set._get_sources(defn_units):
            expr = src_units_set.get_expression()
            self._set_element_in_units(expr, src_units, no_act)
        return

    def evaluate(self):
        """Return the first operand to the power of the second."""
        app = self.xml_parent
        ops = list(app.operands())
        if len(ops) != 2:
            self.wrong_number_of_operands(len(ops), [2])
        return self.eval(ops[0]) ** self.eval(ops[1])

    def _reduce(self):
        """Reduce this expression by evaluating its static parts.

        If the whole expression is static, proceed as normal for an <apply>.
        Otherwise check if the exponent is static and the expression being exponentiated
        is a <ci>.  If so, and the exponent is equal to 2, 3, or 4, convert the expression
        to a multiplication.
        """
        app = self.xml_parent
        bt = app._get_binding_time()
        converted = False
        if bt != BINDING_TIMES.static and self.model.get_option('pe_convert_power'):
            base, expt = list(app.operands())
            expt_bt = app._get_element_binding_time(expt)
            if expt_bt == BINDING_TIMES.static and isinstance(base, mathml_ci):
                expt_val = self.eval(expt)
                if expt_val in [2,3,4]:
                    # Do the conversion
                    app.safe_remove_child(base)
                    operands = [base]
                    for _ in range(1, expt_val):
                        operands.append(base.clone_self())
                        base.variable._used()
                    new_app = mathml_apply.create_new(app, u'times', operands)
                    app.replace_child(app, new_app, app.xml_parent)
                    # Finally, reduce the new expression, just to be safe!
                    new_app._reduce()
                    converted = True
        if not converted:
            # Evaluate this expression as normal
            app._reduce(check_operator=False)

class mathml_root(mathml_operator, mathml_units_mixin):
    """Class representing the MathML <root> operator."""
    def _set_in_units(self, units, no_act=False):
        """Set the units of the application of this operator.

        Set the degree to have units of dimensionless, and the operand to
        have an arbitrary member of its possible units set.

        Where these mean the <apply> doesn't have the given units, wrap it
        in suitable units conversion mathematics.
        """
        app = self.xml_parent
        defn_units_set = app.get_units()
        defn_units = defn_units_set.extract()
        app._add_units_conversion(app, defn_units, units, no_act)
        # Record which member of the set we used
        if not no_act:
            app._cml_units = defn_units
        # Set degree units
        if hasattr(app, u'degree'):
            dimensionless = app.model.get_units_by_name('dimensionless')
            self._set_element_in_units(_child1(app.degree), dimensionless, no_act)
        # Set operand units
        for src_units_set, src_units in defn_units_set._get_sources(defn_units):
            expr = src_units_set.get_expression()
            self._set_element_in_units(expr, src_units, no_act)

    def evaluate(self):
        """
        Return the operand to the given degree, if present.
        Otherwise return the square root of the operand.
        """
        app = self.xml_parent
        ops = list(app.operands())
        if len(ops) != 1:
            self.wrong_number_of_operands(len(ops), [1])
        if hasattr(app, u'degree'):
            degree = self.eval(app.degree)
        else:
            degree = 2
        return self.eval(ops[0]) ** (1/degree)

class mathml_and(mathml_operator, mathml_units_mixin_equalise_operands):
    """Class representing the MathML <and> operator."""
    def evaluate(self):
        """Return the logical conjunction of the operands.

        Evaluates operands in the order given in the file, and will
        short-circuit at the first which evaluates to False.
        """
        app = self.xml_parent
        ops = app.operands()
        value = True
        for operand in ops:
            value = value and self.eval(operand)
            if not value: break
        return value

    def _get_binding_time(self):
        """Return the binding time of the enclosing <apply> element.

        Short-circuit if a static False operand occurs before any dynamic
        operands, returning static.  Otherwise return the least upper bound
        of operand binding times, as usual.
        """
        app = self.xml_parent
        bts = [BINDING_TIMES.static]
        for operand in app.operands():
            bt = app._get_element_binding_time(operand)
            if bt is BINDING_TIMES.static:
                value = self.eval(operand)
                if not value and len(bts) == 1:
                    # Short-circuit
                    break
            else:
                bts.append(bt)
        # Take least upper bound
        return max(bts)

    # TODO: Write a _reduce method

class mathml_or(mathml_operator, mathml_units_mixin_equalise_operands):
    """Class representing the MathML <or> operator."""
    def evaluate(self):
        """Return the logical disjunction of the operands.

        Evaluates operands in the order given in the file, and will
        short-circuit at the first which evaluates to True.
        """
        app = self.xml_parent
        ops = app.operands()
        value = False
        for operand in ops:
            value = value or self.eval(operand)
            if value: break
        return value

    def _get_binding_time(self):
        """Return the binding time of the enclosing <apply> element.

        Short-circuit if a static True operand occurs before any dynamic
        operands, returning static.  Otherwise return the least upper bound
        of operand binding times, as usual.
        """
        app = self.xml_parent
        bts = [BINDING_TIMES.static]
        for operand in app.operands():
            bt = app._get_element_binding_time(operand)
            if bt is BINDING_TIMES.static:
                value = self.eval(operand)
                if value and len(bts) == 1:
                    # Short-circuit
                    break
            else:
                bts.append(bt)
        # Take least upper bound
        return max(bts)

    # TODO: Write a _reduce method

class mathml_leq(mathml_operator, mathml_units_mixin_equalise_operands):
    """Class representing the MathML <leq> operator."""
    def evaluate(self):
        """
        Return True iff the value of the first operand is
        less than or equal to the value of the second.
        """
        app = self.xml_parent
        ops = list(app.operands())
        if len(ops) != 2:
            self.wrong_number_of_operands(len(ops), [2])
        return self.eval(ops[0]) <= self.eval(ops[1])

class mathml_lt(mathml_operator, mathml_units_mixin_equalise_operands):
    """Class representing the MathML <lt> operator."""
    def evaluate(self):
        """
        Return True iff the value of the first operand is
        less than the value of the second.
        """
        app = self.xml_parent
        ops = list(app.operands())
        if len(ops) != 2:
            self.wrong_number_of_operands(len(ops), [2])
        return self.eval(ops[0]) < self.eval(ops[1])

class mathml_geq(mathml_operator, mathml_units_mixin_equalise_operands):
    """Class representing the MathML <geq> operator."""
    def evaluate(self):
        """
        Return True iff the value of the first operand is
        greater than or equal to the value of the second.
        """
        app = self.xml_parent
        ops = list(app.operands())
        if len(ops) != 2:
            self.wrong_number_of_operands(len(ops), [2])
        return self.eval(ops[0]) >= self.eval(ops[1])

class mathml_gt(mathml_operator, mathml_units_mixin_equalise_operands):
    """Class representing the MathML <gt> operator."""
    def evaluate(self):
        """
        Return True iff the value of the first operand is
        greater than the value of the second.
        """
        app = self.xml_parent
        ops = list(app.operands())
        if len(ops) != 2:
            self.wrong_number_of_operands(len(ops), [2])
        return self.eval(ops[0]) > self.eval(ops[1])

class mathml_neq(mathml_operator, mathml_units_mixin_equalise_operands):
    """Class representing the MathML <neq> operator."""
    def evaluate(self):
        """Evaluate the enclosing <apply> element.

        Return True iff the 2 operands are not equal.
        """
        app = self.xml_parent
        ops = list(app.operands())
        if len(ops) != 2:
            self.wrong_number_of_operands(len(ops), [2])

        return (self.eval(ops[0]) != self.eval(ops[1]))

class mathml_eq(mathml_operator, mathml_units_mixin_equalise_operands):
    """Class representing the MathML <eq> operator."""
    def _is_top_level(self):
        """Return True iff the enclosing <apply> is a top-level expression."""
        return self.xml_parent.is_top_level()

    def _set_in_units(self, units, no_act=False):
        """Set the units of the application of this operator.

        If this is a top-level <eq/>, then force the RHS to take the units
        of the LHS.  Otherwise, behave as for other relational operators.
        """
        if self._is_top_level():
            ops = self.xml_parent.operands()
            lhs = ops.next()
            lhs_units = lhs.get_units().extract()
            self._set_element_in_units(lhs, lhs_units, no_act)
            self._set_element_in_units(ops.next(), lhs_units, no_act)
            if not no_act:
                self.xml_parent._cml_units = units
        else:
            super(mathml_eq, self)._set_in_units(units, no_act)

    def evaluate(self):
        """Evaluate the enclosing <apply> element.

        The behaviour depends on whether the enclosing <apply> is a
        top-level expression or not, i.e. whether this is an
        assignment or a comparison.

        If an assignment, evaluate the RHS, assign the value to
        the variable on the LHS, and return it.

        If a comparison, return True iff the 2 operands are equal.
        """
        app = self.xml_parent
        ops = list(app.operands())
        if len(ops) != 2:
            self.wrong_number_of_operands(len(ops), [2])
        
        if self._is_top_level():
            # This is a top-level assignment or ODE
            value = self.eval(ops[1])
            var = app.assigned_variable()
            if app.is_assignment():
                var.set_value(value)
            elif app.is_ode():
                indepvar = var[1].get_source_variable(recurse=True)
                var[0].set_value(value, ode=indepvar)
            else:
                raise EvaluationError("Weird sort of assignment expression.")
        else:
            # This is a comparison
            value = (self.eval(ops[0]) == self.eval(ops[1]))
        return value

    def _get_binding_time(self):
        """Return the binding time of the enclosing <apply> element.

        If this is a top-level expression, then only recurse into the RHS,
        otherwise proceed as normal for an apply.

        There is one further special case: if this is a top-level
        expression and the variable assigned to is annotated to be
        kept in the specialised model, then the expression is dynamic.
        """
        app = self.xml_parent
        if self._is_top_level():
            annotated_as_kept = False
            if app.is_ode():
                DEBUG('partial-evaluator', "BT ODE",
                      map(lambda v: v.fullname(), app.assigned_variable()))
            else:
                DEBUG('partial-evaluator', "BT expr",
                      app.assigned_variable().fullname())
                if app.assigned_variable().pe_keep:
                    annotated_as_kept = True
            ops = list(app.operands())
            if len(ops) != 2:
                self.wrong_number_of_operands(len(ops), [2])
            rhs = ops[1]
            bt = app._get_element_binding_time(rhs)
            if annotated_as_kept:
                bt = BINDING_TIMES.dynamic
        else:
            bt = app._get_binding_time(check_operator=False)
        return bt

    def _reduce(self):
        """Reduce this expression by evaluating its static parts.

        If this is a top-level assignment, then just reduce the RHS.
        Otherwise proceed as normal for an <apply>.
        """
        app = self.xml_parent
        if self._is_top_level():
            ops = list(app.operands())
            if len(ops) != 2:
                self.wrong_number_of_operands(len(ops), [2])
            rhs = ops[1]
            app._reduce_elt(rhs)
        else:
            app._reduce(check_operator=False)

    @property
    def rhs(self):
        """Return the right hand side of this expression.

        Should only be called if we're actually an assignment.
        """
        if self._is_top_level():
            ops = self.xml_parent.operands()
            ops.next()
            return ops.next()
        else:
            raise ValueError("Not an assignment expression.")
    @property
    def lhs(self):
        """Return the left hand side of this expression.

        Should only be called if we're actually an assignment.
        """
        if self._is_top_level():
            ops = self.xml_parent.operands()
            return ops.next()
        else:
            raise ValueError("Not an assignment expression.")

class mathml_rem(mathml_operator, mathml_units_mixin_set_operands):
    """Class representing the MathML <rem> operator."""
    def evaluate(self):
        """Return the remainder when the first operand is divided by the second."""
        app = self.xml_parent
        ops = list(app.operands())
        if len(ops) != 2:
            self.wrong_number_of_operands(len(ops), [2])
        return self.eval(ops[0]) % self.eval(ops[1])

class mathml_logbase(mathml, mathml_units_mixin_container):
    """Class representing the MathML <logbase> element."""
    def evaluate(self):
        """Evaluate this element, by evaluating its child."""
        return self.eval(_child1(self))

class mathml_degree(mathml, mathml_units_mixin_container):
    """Class representing the MathML <degree> element."""
    def evaluate(self):
        """Evaluate this element, by evaluating its child."""
        return self.eval(_child1(self))

class mathml_otherwise(mathml):
    """Class representing the MathML <otherwise> element.
    
    Only defined to make it inherit from mathml.
    """
    pass

class mathml_piece(mathml):
    """Class representing the MathML <piece> element.
    
    Only defined to make it inherit from mathml.
    """
    pass

class mathml_sin(mathml_operator, mathml_units_mixin_set_operands):
    """Class representing the MathML <sin> operator."""
    def evaluate(self):
        """Return the sine of the single operand."""
        app = self.xml_parent
        ops = list(app.operands())
        if len(ops) != 1:
            self.wrong_number_of_operands(len(ops), [1])
        return math.sin(self.eval(ops[0]))

class mathml_cos(mathml_operator, mathml_units_mixin_set_operands):
    """Class representing the MathML <cos> operator."""
    def evaluate(self):
        """Return the cosine of the single operand."""
        app = self.xml_parent
        ops = list(app.operands())
        if len(ops) != 1:
            self.wrong_number_of_operands(len(ops), [1])
        return math.cos(self.eval(ops[0]))

class mathml_tan(mathml_operator, mathml_units_mixin_set_operands):
    """Class representing the MathML <tan> operator."""
    def evaluate(self):
        """Return the tangent of the single operand."""
        app = self.xml_parent
        ops = list(app.operands())
        if len(ops) != 1:
            self.wrong_number_of_operands(len(ops), [1])
        return math.tan(self.eval(ops[0]))

class mathml_arcsin(mathml_operator, mathml_units_mixin_set_operands):
    """Class representing the MathML <arcsin> operator."""
    def evaluate(self):
        """Return the arc sine of the single operand, in radians."""
        app = self.xml_parent
        ops = list(app.operands())
        if len(ops) != 1:
            self.wrong_number_of_operands(len(ops), [1])
        return math.asin(self.eval(ops[0]))

class mathml_arccos(mathml_operator, mathml_units_mixin_set_operands):
    """Class representing the MathML <arccos> operator."""
    def evaluate(self):
        """Return the arc cosine of the single operand, in radians."""
        app = self.xml_parent
        ops = list(app.operands())
        if len(ops) != 1:
            self.wrong_number_of_operands(len(ops), [1])
        return math.acos(self.eval(ops[0]))

class mathml_arctan(mathml_operator, mathml_units_mixin_set_operands):
    """Class representing the MathML <arctan> operator."""
    def evaluate(self):
        """Return the arc tangent of the single operand, in radians."""
        app = self.xml_parent
        ops = list(app.operands())
        if len(ops) != 1:
            self.wrong_number_of_operands(len(ops), [1])
        return math.atan(self.eval(ops[0]))


class mathml_csymbol(mathml_operator, mathml_units_mixin):
    """Class representing the MathML <csymbol> operator.
    
    This is used to represent special-case operations that are treated uniquely by the code generation phase.
    """
    def evaluate(self):
        raise NotImplementedError

    def _set_in_units(self, units, no_act=False):
        """Set the units of the application of this operator, adding a conversion if needed."""
        app = self.xml_parent
        defn_units = app.get_units().extract()
        if defn_units != units:
            self._add_units_conversion(app, defn_units, units, no_act)
        # Store the units
        if not no_act:
            app._cml_units = units
        return


## Don't export module imports to people who do "from pycml import *"
#__all__ = filter(lambda n: type(globals()[n]) is not types.ModuleType, globals().keys())
