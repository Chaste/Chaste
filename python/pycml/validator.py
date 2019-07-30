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

# Validator for CellML 1.0
# Author: Jonathan Cooper

# We want 1/2==0.5
from __future__ import division

import codecs
import optparse
import os
import re
import subprocess
import sys

# Make sure PyCml is on sys.path
pycml_path = os.path.dirname(os.path.realpath(__file__))
sys.path[0:0] = [pycml_path]

# Common CellML processing stuff
import pycml
from pycml import *  # Put contents in the local namespace as well

__version__ = "$Revision$"[11:-2]



class ValidatorError(Exception):
    """Base class for errors trying to run validator stages."""
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)



class CellMLValidator(object):
    def __init__(self, create_relaxng_validator=True):
        """Initialise a validator for CellML files."""
        # Create validator from RELAX NG schema
        self.relaxng_validator = None
        if create_relaxng_validator:
            schema_base = os.path.join(pycml_path, 'cellml1.0')
            run_errors = []
            for klass in [LxmlRelaxngValidator, RvpRelaxngValidator]:
                try:
                    self.relaxng_validator = klass(schema_base)
                except ValidatorError, e:
                    run_errors.append(str(e))
                else:
                    break
            if not self.relaxng_validator:
                msg = '\n\t'.join(["Unable to run a RELAX NG validator.  Please install lxml or rvp."]
                                  + run_errors)
                raise ValidatorError(msg)

    def quit(self):
        """
        Since using __del__ is precarious, we provide this method to allow
        the RVP process to be killed cleanly.  Call it when the validator
        is finished with, or you'll get an interesting error when the program
        terminates (if using RVP).
        """
        if self.relaxng_validator:
            self.relaxng_validator.quit()
        return
    
    @staticmethod
    def setup_logging(show_errors=True, error_stream=sys.stderr,
                      show_warnings=True, warning_stream=sys.stderr,
                      space_errors=False, loglevel=logging.WARNING,
                      **kwargs):
        """Set up loggers for validation errors/warnings.
        
        Set show_errors or show_warnings to False to suppress the output of
        validation error or warning messages, respectively.  When not
        suppressed, the messages will be output to the streams given by
        error_stream and warning_stream; these should be file-like objects.

        If space_errors is True, a blank line will be inserted between
        each message.
        """
        logger = logging.getLogger('validator')
        logger.setLevel(loglevel)
        if space_errors:
            formatter = logging.Formatter(fmt="%(message)s\n")
        else:
            formatter = logging.Formatter(fmt="%(message)s")
        error_handler = logging.StreamHandler(error_stream)
        error_handler.setLevel(logging.ERROR)
        error_handler.setFormatter(formatter)
        warning_handler = logging.StreamHandler(warning_stream)
        warning_handler.addFilter(OnlyWarningsFilter())
        warning_handler.setFormatter(formatter)
        if not show_errors:
            error_handler.setLevel(logging.CRITICAL)
        if not show_warnings:
            warning_handler.setLevel(logging.CRITICAL)
        logger.addHandler(error_handler)
        logger.addHandler(warning_handler)
        return error_handler, warning_handler
    
    @staticmethod
    def cleanup_logging(handlers):
        """Flush logger & remove handlers."""
        logger = logging.getLogger('validator')
        for handler in handlers:
            handler.flush()
            logger.removeHandler(handler)

    def validate(self, source, return_doc=False, assume_valid=False, **kw):
        """Validate the given document.

        source should be a file-like object, URI, local file name,
        or '-' for standard input.  If a file-like object, it must support
        the seek method to reset it.

        If return_doc is True then the result is a tuple (valid, document),
        where if valid==True then document is an Amara binding of the CellML
        document.
        Otherwise just return True iff the document is valid.

        If xml_context is True, then the failing XML tree will be displayed
        with every units error.

        The assume_valid option allows you to skip RELAX NG validation, along
        with many of the checks in the Python code.  This is useful for speeding
        transformation of models that are known to pass these checks.
        
        See cellml_model.validate and setup_logging for other keyword arguments.
        """
        logging_info = CellMLValidator.setup_logging(**kw)
        
        # Validate against RELAX NG schema
        DEBUG('validator', 'CellML Validator version', __version__)
        res = True
        # Get stream of CellML document
        if source == '-':
            stream = sys.stdin
        elif hasattr(source, 'read'):
            stream = source
        elif bt.Uri.IsAbsolute(source):
            stream = bt.Url.UrlOpen(source)
        else:
            if not os.path.isfile(source):
                res = False
                logging.getLogger('validator').error('File ' + source + ' does not exist.')
            else:
                stream = file(source, 'r')
        # Parse & validate
        if res and not assume_valid:
            DEBUG('validator', 'Starting RELAX NG validation')
            res = self.relaxng_validator.validate(stream)
            if stream == source:
                source.seek(0)
            elif not stream == sys.stdin:
                stream.close()
            DEBUG('validator', 'Finished RELAX NG:', res)

        # Check further rules that can't be expressed by a (RELAX NG) schema.
        # We use our own Python code for this.
        if res:
            DEBUG('validator', 'Loading model with Amara')
            doc = amara_parse_cellml(source)
            DEBUG('validator', 'Validating loaded model')
            res = doc.model.validate(assume_valid=assume_valid, **kw)
            DEBUG('validator', 'Validation complete:', res)
        else:
            doc = None

        # Flush logger & remove handlers
        CellMLValidator.cleanup_logging(logging_info)

        # Return result
        if return_doc:
            return (res, doc)
        else:
            if doc:
                doc.model.clean_up()
            return res


class LxmlRelaxngValidator(object):
    """
    A RELAX NG validator built on top of lxml (http://lxml.de/validation.html#relaxng).
    Can validate against schemas written in the XML syntax.
    """
    def __init__(self, schemaBase):
        """Initialise the RELAX NG validator.
        
        Parses the schema into memory, and constructs lxml's validator object.
        We are passed the path to the schema with no extension.
        """
        try:
            from lxml import etree
        except ImportError, e:
            raise ValidatorError("Unable to import lxml: " + str(e))
        fp = open(schemaBase + '.rng', 'r')
        schema_doc = etree.parse(fp)
        self._validator = etree.RelaxNG(schema_doc)
    
    def validate(self, stream):
        """Validate an XML document, returning a boolean.
        
        stream should be a file-like object containing the document to be validated.
        Returns True iff the document was valid."""
        from lxml import etree
        doc = etree.parse(stream)
        res = self._validator.validate(doc)
        # Report error details via the logger
        logger = logging.getLogger('validator')
        for e in self._validator.error_log:
            logger.error(e)
        return res

    def quit(self):
        """Providing for compatibility with RvpRelaxngValidator; does nothing."""
        pass


class RvpRelaxngValidator(object):
    """
    A RELAX NG validator built on top of RVP (http://www.davidashen.net/rnv.html).
    Can validate against schemas written in the compact syntax.
    """
    class RvpProtocolError(ValidatorError):
        """Raised if the response from RVP is not understood."""
        pass

    def __init__(self, schemaBase):
        """Initialise the RELAX NG validator.
        Launches RVP as a parallel process.
        schema_filename should be the name of a file containing the RELAX NG schema, in compact syntax.
        """
        self._ws = re.compile('[^\t\n\r ]')
        schema_filename = schemaBase + '.rnc'
        # Launch RVP
        try:
            self._rvp_pipe = subprocess.Popen(['rvp', schema_filename],
                                              stdin=subprocess.PIPE, stdout=subprocess.PIPE, close_fds=True)
        except OSError, e:
            raise self.RvpProtocolError("Failed to run rvp for CellML syntax validation: " + str(e))
        # We use os.read & os.write, so store file descriptors
        self._rvpin, self._rvpout = self._rvp_pipe.stdin.fileno(), self._rvp_pipe.stdout.fileno()
        # Import Expat parser module
        import xml.parsers.expat
        self.expat = xml.parsers.expat

##    def __del__(self):
##        """
##        Tell our RVP process to quit.
##        This doesn't work well, since __del__ isn't necessarily called at
##        program exit.  Hence, there is a manual quit() method.
##        """
##        self.quit()
    def quit(self):
        """Call this method when the validator is finished with.
        It terminates the associated RVP process cleanly.
        Failure to do so will probably result in an error when your program exits.
        """
        self._send('quit')
        return self._resp()

    def validate(self, stream):
        """Validate an XML document, returning a boolean.
        
        stream should be a file-like object containing the document to
        be validated.
        Returns True iff the document was valid.
        """
        # Initialise
        self._text = ''
        self._errors, self._error_messages = False, []
        self._prevline, self._prevcol = -1, -1
        self._pat = self._start()
        # Create parser
        parser = self.expat.ParserCreate(namespace_separator=':')
        parser.StartElementHandler = self._start_element
        parser.EndElementHandler = self._end_element
        parser.CharacterDataHandler = self._characters
        self._parser = parser
        # Parse & validate
        try:
            self._parser.ParseFile(stream)
        except self.expat.ExpatError, e:
            self._errors = True
            self._error_found(e.lineno, e.offset, str(e))
        # Any errors
        return not self._errors

    def get_validation_errors(self):
        """
        Return the list of all errors found while validating
        the current file.
        """
        return self._error_messages

    def _error_found(self, line, col, msg):
        report = "%d,%d:%s" % (line, col, msg.strip())
        logging.getLogger('validator').error(report)
        self._error_messages.append(report)

    # RVP protocol methods
    def _start_tag_open(self, cur, name):
        self._send('start-tag-open '+cur+' '+name)
        return self._resp()
    def _attribute(self, cur, name, val):
        self._send('attribute '+cur+' '+name+' '+val)
        return self._resp()
    def _start_tag_close(self, cur, name):
        self._send('start-tag-close '+cur+' '+name)
        return self._resp()
    def _end_tag(self, cur, name):
        self._send('end-tag '+cur+' '+name)
        return self._resp()
    def _textonly(self, cur, text):
        self._send('text '+cur+' '+text)
        return self._resp()
    def _mixed(self, cur, text):
        """
        In mixed content, whitespace is discarded, and any non-whitespace
        is counted as equal.
        """
        if self._ws.search(text):
            self._send('mixed '+cur+' .')
            return self._resp()
        else:
            return cur
    def _start(self, grammar = '0'):
        self._send('start '+grammar)
        return self._resp()

    # Low-level communication with RVP
    def _send(self, s):
        """
        Terminate string with zero, encode in UTF-8 and send to RVP.
        """
        os.write(self._rvpin, s.encode('UTF-8') + '\0')
    def _recv(self):
        """
        Receive a zero-terminated response from RVP; drop zero byte.
        """
        s = ''
        while True:
            # 16 is a good buffer length for ok responses; errors
            # should be rare
            s = s + os.read(self._rvpout, 16)
            if s[-1] == '\0': break
        return s[:-1]
    def _resp(self):
        """
        Get a reply from RVP.
        If an error occurs, log the message.
        Return the current pattern value.
        """
        r = self._recv().split(' ', 3)
        if r[0] == 'ok': return r[1]
        if r[0] == 'error':
            self._errors = True
            if r[3] != '': # Only report if we have a message
                line = self._parser.CurrentLineNumber
                col  = self._parser.CurrentColumnNumber
                if line != self._prevline or col != self._prevcol:
                    # One report per file position
                    self._error_found(line, col, r[3])
                    self._prevline, self._prevcol = line, col
            return r[1]
        if r[0] == 'er':
            self._errors = True
            return r[1]
        # Unknown response
        raise self.RvpProtocolError, "unexpected response '"+r[0]+"'"
                    

    # Expat handlers
    def _flush_text(self):
        """
        Apparently Expat doesn't concatenate text nodes, so we do it
        manually; the CharDataHandler collects the text, and this
        method passes it to the validator.
        """
        if self._text:
            if self._ismixed:
                self._pat = self._mixed(self._pat, self._text)
            else:
                self._pat = self._textonly(self._pat, self._text)
            self._text = ''
    def _start_element(self, name, attrs):
        self._ismixed = True
        self._flush_text()
        self._pat = self._start_tag_open(self._pat, name)
        self._ismixed = False
        for n, v in attrs.items():
            self._pat = self._attribute(self._pat, n, v)
        self._pat = self._start_tag_close(self._pat, name)
    def _end_element(self, name):
        self._flush_text()
        self._pat = self._end_tag(self._pat, name)
        self._ismixed = True
    def _characters(self, data):
        self._text = self._text + data






######################################################################
#                        Convenience functions                       #
######################################################################

def check_repo(repo_dir = '../../models/all_from_repository',
               model_suffix = 'xml',
               invalid_if_warnings = False,
               compare=True):
    """
    Validate every model in the CellML repository, and return a list
    of invalid models.
    
    If compare is False, writes errors & warnings to log files in the
    same folder as the models, otherwise compares the output to log
    files already present, and notes differences.
    
    Displays total run time.
    """
    def close_log_file(stream, filename):
        stream.close()
        try:
            size = os.path.getsize(filename)
            if size == 0:
                os.remove(filename)
        except OSError:
            pass
    import glob, time, gc
    start_time = time.time()
    v = CellMLValidator()
    invalid = []
    files = glob.glob(repo_dir + '/*.' + model_suffix)
    files.sort()
    for filename in files:
        model = os.path.basename(filename)[:-4]
        fn = os.path.splitext(filename)[0]
        warnfn, errfn = fn + '_warnings.log', fn + '_errors.log'
        if compare:
            warn_stream = StringIO()
            err_stream = StringIO()
        else:
            warn_stream = open(warnfn, 'w')
            err_stream = open(errfn, 'w')
        print "Checking model",model,"at",time.strftime("%X %x"),
        sys.stdout.flush()
        valid = v.validate(filename, error_stream=err_stream,
                           warning_stream=warn_stream,
                           invalid_if_warnings=invalid_if_warnings)
        if not valid:
            print max(1,4 - (5+len(model))//8) * '\t', "X"
        else:
            print
        if compare:
            compare_output_files(warn_stream, warnfn)
            compare_output_files(err_stream, errfn)
        else:
            close_log_file(err_stream, errfn)
            close_log_file(warn_stream, warnfn)
        if not valid:
            invalid.append(model)
        gc.collect()
    elapsed_time = time.time() - start_time
    mins,secs = int(elapsed_time//60), int(elapsed_time%60)
    print len(files),"models checked in",mins,"minutes",secs,"seconds."
    print len(invalid),"invalid"
    v.quit()
    return invalid

def compare_output_files(new_stream, old_filename):
    def save_new_output():
        nfp = open(old_filename + '-new', 'w')
        nfp.write(new_stream.getvalue())
        nfp.close()
    new_stream.seek(0, 2)
    new_len = new_stream.tell()
    try:
        fp = open(old_filename, 'r')
    except IOError:
        if new_len > 0:
            print "Log file", old_filename, "doesn't exist,", \
                  "but we have new output"
            try:
                ofp = open(os.path.join(os.path.dirname(old_filename),
                                        'new'), 'a')
                print >>ofp, "Log file", old_filename, "doesn't exist,", \
                      "but we have new output"
                ofp.close()
            except IOError:
                pass
            save_new_output()
        return
    new_stream.seek(0)
    new_lines = set(new_stream.readlines())
    old_lines = set(fp.readlines())
    if old_lines != new_lines:
        print "Output set differs from log file", old_filename
        print "Lines added:", new_lines - old_lines
        print "Lines removed:", old_lines - new_lines
        try:
            ofp = open(os.path.join(os.path.dirname(old_filename), 'new'), 'a')
            print >>ofp, "Output set differs from log file", old_filename
            print >>ofp, "Lines added:", new_lines - old_lines
            print >>ofp, "Lines removed:", old_lines - new_lines, "\n"
            ofp.close()
        except IOError:
            pass
        save_new_output()
    new_stream.close()
    fp.close()
    return

######################################################################
#                    For running as an executable                    #
######################################################################

def get_options(args):
    """get_options(args):
    Process our command-line options.

    args is a list of options & positional arguments.
    """
    usage = 'usage: %prog [options] <cellml file or URI> ...'
    parser = optparse.OptionParser(version="%%prog %s" % __version__,
                                   usage=usage)
    parser.add_option('-o', dest='outfilename', metavar='OUTFILE',
                      help='write *all* output to OUTFILE (overrides -e and -w)')
    parser.add_option('-e', '--error-file',
                      dest='errfilename', metavar='ERRFILE',
                      default='stderr',
                      help='write errors to ERRFILE [default %default]')
    parser.add_option('-w', '--warning-file',
                      dest='warnfilename', metavar='WARNFILE',
                      default='stderr',
                      help='write warnings to WARNFILE [default %default]')
    parser.add_option('-q', '--quiet',
                      dest='quiet', action='store_true', default=False,
                      help="don't display any output, just set exit status")
    parser.add_option('--no-warnings', '--only-errors',
                      dest='show_warnings', action='store_false',
                      default=True,
                      help="don't display warning messages")
    parser.add_option('-s', '--space-messages',
                      dest='space_errors', action='store_true',
                      default=False,
                      help="print a blank line after every warning/error message")
    parser.add_option('-x', '--xml-context',
                      dest='xml_context', action='store_true', default=False,
                      help="display the MathML tree of any expression with invalid units")
    parser.add_option('-u', '--warn-on-unit-conversions',
                      action='store_true', default=False,
                      help="generate a warning if unit conversions are required")
    parser.add_option('--Wu', '--warn-on-units-errors',
                      action='store_true', default=False,
                      dest='warn_on_units_errors',
                      help="give a warning instead of an error for"
                      " dimensional inconsistencies")
    parser.add_option('-i', '--interactive',
                      action='store_true', default=False,
                      help="use with python -i to enter an interactive session"
                      " after validation")
    parser.add_option('-d', '--debug', action='store_true', default=False,
                      help="output debug info to stderr")
    parser.add_option('--profile', action='store_true', default=False,
                      help="turn on profiling of PyCml")

    options, args = parser.parse_args(args)
    if len(args) < 1:
        parser.error("an input CellML file must be specified")
    return options, args


def run():
    # Validate all files specified on the command line
    options, files = get_options(sys.argv[1:])
    validator = CellMLValidator()

    # Open output streams
    if not options.outfilename is None:
        out_s = open_output_stream(options.outfilename)
        err_s = warn_s = out_s
    else:
        out_s  = sys.stdout
        err_s  = open_output_stream(options.errfilename)
        warn_s = open_output_stream(options.warnfilename)

    # Keyword arguments for validator
    kwargs = {'show_warnings': options.show_warnings,
              'warning_stream': warn_s,
              'show_errors': not options.quiet,
              'error_stream': err_s,
              'space_errors': options.space_errors,
              'xml_context': options.xml_context,
              'warn_on_units_errors': options.warn_on_units_errors,
              'check_for_units_conversions': options.warn_on_unit_conversions}
    
    if options.debug:
        formatter = logging.Formatter(fmt="%(name)s %(asctime)s: %(message)s")
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(formatter)
        handler.addFilter(OnlyDebugFilter())
        logging.getLogger().addHandler(handler)
        logging.getLogger().setLevel(logging.DEBUG)
        kwargs['loglevel'] = logging.DEBUG
    
    # Validate files
    result = True
    for f in files:
        if not options.quiet:
            print >>out_s, "Validating file", f, "against CellML 1.0"
        res = validator.validate(f, **kwargs)
        if not options.quiet:
            print >>out_s, "File is%s valid CellML 1.0" % ((' NOT','')[res])
        result = result and res
    
    # Close output streams
    close_output_stream(out_s)
    close_output_stream(err_s)
    close_output_stream(warn_s)

    # Tidy up validator
    validator.quit()

    # Set exit status
    if not options.interactive:
        sys.exit(not result)


if __name__ == '__main__':
    if '--profile' in sys.argv:
        import time, cProfile
        profile_name = '/tmp/pycml-profile-%f-%d' % (time.time(), os.getpid()) 
        cProfile.run('run()', profile_name)
    else:
        run()
