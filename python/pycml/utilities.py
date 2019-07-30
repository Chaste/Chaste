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
Various helpful utility functions/classes used by PyCml.
"""

import logging
import sys
from xml.dom import Node # For nodeType values

import amara

from _enum import Enum

__all__ = ['OnlyWarningsFilter', 'OnlyDebugFilter', 'OnlyTheseSourcesFilter', 'NotifyHandler',
           'DEBUG', 'LOG',
           'Colourable', 'DFS', 'Sentinel', 'unitary_iterator',
           'amara_parse', 'element_path', 'element_path_cmp', 'element_xpath', 'brief_xml',
           'call_if', 'max_i', 'prid', 'add_dicts',
           'open_output_stream', 'close_output_stream']

################################################################################
#                                                                              #
#                              Logging                                         #
#                                                                              #
################################################################################

class OnlyWarningsFilter(logging.Filter):
    """A filter that only passes warning messages."""
    def filter(self, rec):
        return (logging.WARNING <= rec.levelno < logging.ERROR)


class OnlyDebugFilter(logging.Filter):
    """A filter that only passes debug messages."""
    def filter(self, rec):
        return (logging.DEBUG <= rec.levelno < logging.INFO)


class OnlyTheseSourcesFilter(logging.Filter):
    """A filter that only emits messages from the given sources."""
    def __init__(self, sources):
        logging.Filter.__init__(self)
        self.__sources = sources
    def filter(self, rec):
        return rec.name in self.__sources


class NotifyHandler(logging.Handler):
    """
    A logging handler that just notes if any messages are logged.
    """
    def __init__(self, level=logging.NOTSET):
        logging.Handler.__init__(self, level=level)
        self.reset()

    def emit(self, record):
        self.messages = True

    def reset(self):
        """Reset the handler, as if no messages have occurred."""
        self.messages = False


def DEBUG(facility, *args):
    """Log a debug message to facility.

    Arguments are treated as for the print statement.
    """
    logger = logging.getLogger(facility)
    if logger.isEnabledFor(logging.DEBUG):
        logger.debug(' '.join(map(str, args)))


def LOG(facility, level, *args):
    """Log a message to facility with the given level.

    Arguments are treated as for the print statement.
    """
    logger = logging.getLogger(facility)
    if logger.isEnabledFor(level):
        logger.log(level, ' '.join(map(str, args)))

################################################################################
#                                                                              #
#                        Miscellaneous classes                                 #
#                                                                              #
################################################################################

# Useful constants for depth-first search
DFS = Enum('White', 'Gray', 'Black')

class Colourable(object):
    """
    A mixin class for objects that have a colour attribute, and so support
    a depth-first search.
    """
    def __init__(self, *args, **kwargs):
        super(Colourable, self).__init__(*args, **kwargs)
        self.clear_colour()
    
    def set_colour(self, colour):
        self._cml_colour = colour
    
    def get_colour(self):
        return self._cml_colour
    
    def clear_colour(self):
        self._cml_colour = DFS.White


class Sentinel(object):
    """A simple true/false store that can be used as a default parameter value in recursive calls.
    
    This provides a slightly nicer looking alternative to code such as:
    
    def f(done=[False]):
       if not done[0]:
           # do some stuff
           if xyz: done[0] = True
           f(done)

    The list can be replaced with a Sentinel instance.
    """
    def __init__(self, tf=False):
        """Create a new, unset sentinel."""
        self._tf = tf
    
    def __nonzero__(self):
        """Test whether the sentinel has been set."""
        return self._tf
    
    def set(self, tf=True):
        """Set the sentinel."""
        self._tf = tf


class unitary_iterator(object):
    """An iterator over a single item."""
    def __init__(self, start):
        self.curr = start
        return

    def __iter__(self):
        return self

    def next(self):
        if not self.curr:
            raise StopIteration()
        result = self.curr
        self.curr = None
        return result


################################################################################
#                                                                              #
#                        XML-related functions                                 #
#                                                                              #
################################################################################

def amara_parse(source, uri=None, rules=None, binderobj=None,
                prefixes=None):
    """Convenience function for parsing XML.
    
    Works just as amara.parse, except that if source is '-' then
    it reads from standard input.
    """
    if source == '-':
        return amara.parse(sys.stdin, uri=uri, rules=rules,
                           binderobj=binderobj, prefixes=prefixes)
    else:
        return amara.parse(source, uri=uri, rules=rules,
                           binderobj=binderobj, prefixes=prefixes)


def element_path(elt):
    """Find the path from the root element to this element."""
    if hasattr(elt, 'xml_parent'):
        idx = 0
        for child in elt.xml_parent.xml_children:
            if getattr(child, 'nodeType', None) == Node.ELEMENT_NODE:
                idx += 1
                if child is elt:
                    break
        return element_path(elt.xml_parent) + [idx]
    else:
        return []
    

def element_path_cmp(e1, e2):
    """Compare 2 elements by comparing their paths from the root element."""
    return cmp(element_path(e1), element_path(e2))


def element_xpath(elt):
    """Return an xpath expression that will select this element."""
    indices = element_path(elt)
    xpath = u'/*[' + u']/*['.join(map(str, indices)) + u']'
    return xpath


def brief_xml(elt):
    """Print a more concise version of elt.xml() that omits all attributes."""
    s = ''
    if getattr(elt, 'nodeType', None) == Node.ELEMENT_NODE:
        children = getattr(elt, 'xml_children', [])
        if children:
            s += '<' + elt.localName + '>'
            for child in children:
                s += brief_xml(child)
            s += '</' + elt.localName + '>'
        else:
            s += '<' + elt.localName + '/>'
    else:
        s += str(elt)
    return s


################################################################################
#                                                                              #
#                       Miscellaneous functions                                #
#                                                                              #
################################################################################

def call_if(tf, callable, *args, **kwargs):
    """Call the given callable with the given arguments iff tf is True."""
    if tf:
        return callable(*args, **kwargs)


def max_i(it):
    """Find the maximum entry of an iterable, and return it with its index.

    Returns (i, m) where i is the index of m, the maximum entry of `it`.
    """
    idx, m = None, None
    for i, val in enumerate(it):
        if m is None or val > m:
            m, idx = val, i
    return idx, m


def prid(obj, show_cls=False):
    """Get the id of an object as a hex string, optionally with its class/type."""
    if obj is None:
        s = 'None'
    else:
        s = hex(id(obj))
        if show_cls:
            s += str(type(obj))
    return s


def add_dicts(r, *ds):
    """Add multiple dictionaries together.

    Updates the first input dictionary by adding values from
    subsequent inputs.  Produces a dictionary with keys taken from the
    input dictionaries, and values being the sum of the corresponding
    values in all inputs.
    Assumes values are numeric.
    """
    for d in ds:
        for k, v in d.iteritems():
            r[k] = r.get(k, 0) + v


def open_output_stream(fname):
    """Open fname for output.
    
    fname should be a local filename, or '-' for standard output.
    Additionally, the names 'stdout' and 'stderr' have the usual
    special meanings.
    """
    if fname == '-' or fname == 'stdout':
        stream = sys.stdout
    elif fname == 'stderr':
        stream = sys.stderr
    else:
        stream = open(fname, 'w')
    return stream


def close_output_stream(stream):
    """
    Close the given output stream, unless it's one of the standard streams
    (i.e. sys.stdout or sys.stderr).
    Note that closing a stream multiple times is safe.
    """
    if not stream is sys.stdout and not stream is sys.stderr:
        stream.close()
