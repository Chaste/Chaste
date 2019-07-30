#!/usr/bin/env python

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

# mod_python interface to the CellML 1.0 validator.
# Author: Jonathan Cooper

# CGI library
import cgi
import cgitb
# Send exceptions to the browser
cgitb.enable()

from cStringIO import StringIO
import os, sys

# So we can find rvp
os.environ['PATH'] = os.environ['PATH'] + ':/usr/local/bin'

# Helpful output functions

def _footer():
    """
    Print the page footer.
    """
    print """
  <hr />
  <a href="./">CellML Tools index</a>
</body>

</html>
"""
    return

def _error(msg):
    """
    Print the error message, and exit.
    """
    print '<p class="error">%s</p>' % msg
    _footer()
    sys.exit()
    return

def _flush():
    sys.stdout.flush()
    sys.stderr.flush()


def _htmlescape(string):
    """
    Convert problematic characters in string to use HTML entities.
    Handles angle brackets, double quotes and ampersand.
    """
    for char, rep in [('&', 'amp'), ('<', 'lt'), ('>', 'gt'), ('"', 'quot')]:
        string = string.replace(char, '&%s;' % rep)
    return string

def _nl2br(string):
    "Convert newline characters to <br /> tags."
    return string.replace('\n', '<br />\n')

import re
_sp = re.compile(r'^(  +)', re.M)
def _respect_spaces(string):
    "Convert runs of spaces to use &nbsp;"
    def repl(match):
        return '&nbsp;' * len(match.group(0))
    return _sp.sub(repl, string)

def _makehtml(string):
    "Take an error/warning message and return an HTML version."
    return _respect_spaces(_nl2br(_htmlescape(string)))


# Get form data
form = cgi.FieldStorage()

# Determine user's prefered style
style_pref = form.getfirst('style', 'dark')
if style_pref == 'plain':
    stylesheet = """
    .errors   { font-family: monospace; }
    .warnings { font-family: monospace; }"""
else:
    stylesheet = """
    body { background-color: #222222; color: #bbbbbb; }
    .errors   { color: red; font-family: monospace; }
    .warnings { color: orange; font-family: monospace; }"""
stylesheet = stylesheet + """
    .valid   { color: green;
               text-align: center; font-size: larger;
             }
    .invalid { color: red;
               text-align: center; font-size: larger;
             }
    .error   { color: red; }"""


# Start page output
print """Content-Type: text/html

<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE html 
     PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
          "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">

<head>
  <title>CellML Tools - Validator</title>
  <style type="text/css">%s
  </style>
</head>

<body>
  <h1>Validation results against CellML 1.0.</h1>
""" % stylesheet


# Process remaining user form input
if not form.has_key('cellml'):
    # Error - no file
    _error('You must upload a CellML file.')

cellml = form['cellml']
if not cellml.file:
    # Error - not a file
    _error('A CellML file must be uploaded.')

show_errors = bool(form.getfirst('show_errors', False))
show_warnings = bool(form.getfirst('show_warnings', False))
xml_context = bool(form.getfirst('xml_context', False))
check_units_conversions = bool(form.getfirst('units_conv', False))
units_warn = bool(form.getfirst('units_warn', False))


## Debugging
#print """<p>
#show_warnings=%s<br />
#show_errors=%s
#</p>""" % (show_warnings, show_errors)

#print "<code>%s</code>" % str(dir(cellml.file))
#print "<p>%s</p>" % cellml.file.readline().replace('<', '&lt;')

#cgi.print_directory()
#cgi.print_environ()



# Load validator module and create validator
import validator as val_mod
validator = val_mod.CellMLValidator()

# Create StringIO objects to collect any warning/error messages,
# if required
if show_warnings:
    warnings = StringIO()
else:
    warnings = None
if show_errors:
    errors = StringIO()
else:
    errors = None

print """
  <h2>Validating file '%s'</h2>
""" % cellml.filename

# Validate the file.
_flush()
valid = validator.validate(
    cellml.file,
    show_errors    = show_errors,
    show_warnings  = show_warnings,
    error_stream   = errors,
    warning_stream = warnings,
    space_errors   = True,
    xml_context    = xml_context,
    warn_on_units_errors=units_warn,
    check_for_units_conversions=check_units_conversions)
validator.quit()

# Show validation result
if valid:
    style = 'valid'
    msg   = 'This file is valid CellML 1.0.'
else:
    style = 'invalid'
    msg   = 'This file is NOT valid CellML 1.0.'
        
print """
  <h3>Validation Result</h3>
  <p class="%s">%s</p>
""" % (style, msg)

# Display any error messages
if show_errors and not valid:
    print """
  <h3>Errors:</h3>
  <div class="errors">
%s  </div>
""" % _makehtml(errors.getvalue())

# Add warnings to the HTML output
if show_warnings and warnings.getvalue():
    print """
  <h3>Warnings:</h3>
  <div class="warnings">
%s  </div>
""" % _makehtml(warnings.getvalue())

# Page footer
_footer()



