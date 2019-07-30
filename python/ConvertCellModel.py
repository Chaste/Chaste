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

"""
A little helper script to facilitate calling PyCml for common Chaste usage scenarios.
"""

import optparse
import os
import subprocess
import sys
import tempfile


# Parse command-line options
class OptionParser(optparse.OptionParser):
    def __init__(self, *args, **kwargs):
        if 'our_short_options' in kwargs:
            self.__our_short_options = kwargs['our_short_options']
            del kwargs['our_short_options']
        else:
            self.__our_short_options = []
        self.__our_short_options.append('-h')
        optparse.OptionParser.__init__(self, *args, **kwargs)
        
    def _process_args(self, largs, rargs, values):
        """
        Override to catch BadOption errors and pass these options on to PyCml.
        """
        while rargs:
            arg = rargs[0]
            if arg == "--":
                del rargs[0]
                return
            elif arg[0:2] == "--":
                try:
                    self._process_long_opt(rargs, values)
                except optparse.BadOptionError:
                    if '=' in arg:
                        # _process_long_opt puts the value as the next option
                        del rargs[0]
                    largs.append(arg)
            elif arg[0] == "-":
                if arg in self.__our_short_options:
                    self._process_short_opts(rargs, values)
                else:
                    largs.append(rargs.pop(0))
            elif self.allow_interspersed_args:
                largs.append(arg)
                del rargs[0]
            else:
                return                  # stop now, leave this arg in rargs
usage = 'usage: %prog [options] <file1.cellml> ...'
parser = OptionParser(usage=usage, our_short_options=['-y'])
parser.add_option('--opt', action='store_true', default=False,
                  help="apply default optimisations to all generated code")
parser.add_option('--normal', action='store_true', default=False,
                  help="generate standard cell model code, suitable for"
                  " use with Chaste's ODE solvers")
parser.add_option('--cvode', action='store_true', default=False,
                  help="generate a subclass of AbstractCvodeCell,"
                  " i.e. that can be simulated with CVODE")
parser.add_option('--cvode-data-clamp', action='store_true', default=False,
                  help="generate a subclass of AbstractCvodeCellWithDataClamp,"
                  " i.e. that can be simulated with CVODE, and clamped to experimental data")
parser.add_option('--backward-euler', action='store_true', default=False,
                  help="generate a version of the cell model that can be"
                  " solved using a backward Euler method.  Requires the"
                  " presence of a .out file accompanying the CellML.")
parser.add_option('--rush-larsen', action='store_true', default=False,
                  help="generate a version of the cell model that can be"
                  " solved using the Rush-Larsen method.")
parser.add_option('--grl1', action='store_true', default=False,
                  help="generate a version of the cell model that can be"
                  " solved using the GRL1 method.")     
parser.add_option('--grl2', action='store_true', default=False,
                  help="generate a version of the cell model that can be"
                  " solved using the GRL2 method.")           
parser.add_option('--output-dir', action='store',
                  help="directory to place output files in")
parser.add_option('--show-outputs', action='store_true', default=False,
                  help="don't actually run PyCml, just show what files would"
                  " be generated, one per line")
parser.add_option('--config-file',
                  action='store',
                  help="pathname of configuration file.  This can be used to"
                  " override options supplied on the command line, and will"
                  " also be passed on to PyCml itself, as will any options"
                  " not listed above.")
parser.add_option('-y', '--dll', '--dynamically-loadable',
                  dest='dynamically_loadable',
                  action='store_true', default=False,
                  help="add code to allow the model to be compiled to a "
                  "shared library and dynamically loaded.  If this option is "
                  "given, only one type of output will be generated (so you "
                  "can't combine, e.g. --cvode --normal).")
parser.add_option('--assume-valid',
                  action='store_true', default=False,
                  help="skip some of the model validation checks")
parser.add_option('-q', '--quiet', action='store_true', default=False,
                  help="don't print the command(s) being run, and pass the option through to pycml")
options, args = parser.parse_args()

option_names = ['opt', 'normal', 'cvode', 'cvode_data_clamp', 'backward_euler', 'rush_larsen', 'grl1', 'grl2']
def arg2name(arg):
    return str(arg)[2:].replace('-', '_')


# Use external PyCml if requested
if 'PYCML_DIR' in os.environ and os.path.isdir(os.environ['PYCML_DIR']):
    if not options.show_outputs:
        print 'Using external PyCml from PYCML_DIR =', os.environ['PYCML_DIR']
    pycml_dir = os.environ['PYCML_DIR']
else:
    pycml_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'pycml')

# Options that we will supply to PyCml anyway
essential_options = ['--conf=' + os.path.join(pycml_dir, 'config.xml'),
                     '--use-chaste-stimulus',
                     '--convert-interfaces',
                     '--warn-on-unit-conversions',
                     '--Wu']
#essential_options.append('--profile')
# Options supplied if the user doesn't give a config file
default_options = []


# Read further arguments from config file?
if options.config_file:
    # Parse the config file and extract any options
    import xml.dom.minidom
    def getText(nodelist):
        rc = []
        for node in nodelist:
            if node.nodeType == node.TEXT_NODE:
                rc.append(node.data)
        return ''.join(rc)
    config_doc = xml.dom.minidom.parse(options.config_file)
    config_modified = False
    cl_args_elt = config_doc.getElementsByTagName('command_line_args')
    if cl_args_elt:
        cl_args_elt = cl_args_elt[0]
        arg_elts = cl_args_elt.getElementsByTagName('arg')
        config_args = []
        # Strip from the file any arguments only understood by this script
        for arg_elt in arg_elts:
            arg = getText(arg_elt.childNodes)
            config_args.append(arg)
            if arg2name(arg) in option_names:
                cl_args_elt.removeChild(arg_elt)
                config_modified = True
        if config_modified:
            # If the config file supplied such arguments, then pretend there weren't
            # any on the command line, since we can't turn options off.
            for option in option_names:
                setattr(options, option, False)
            # An empty command_line_args isn't allowed
            for node in cl_args_elt.childNodes:
                if node.nodeType == node.ELEMENT_NODE:
                    break
            else:
                for node in config_doc.childNodes:
                    if node.nodeType == node.ELEMENT_NODE:
                        node.removeChild(cl_args_elt)
        # Parse additional options
        options, extra_args = parser.parse_args(config_args, options)
        args.extend(extra_args)
    if config_modified:
        # Write a new config file
        fp, config_path = tempfile.mkstemp(suffix='.xml', text=True)
        config_file = os.fdopen(fp, 'w')
        config_doc.writexml(config_file)
        config_file.close()
        essential_options.append('--conf=' + config_path)
    else:
        # Just pass on the one we were given
        essential_options.append('--conf=' + options.config_file)
    del config_doc


# If no options supplied, default to --normal --opt
number_of_options = len(filter(None, [getattr(options, opt_name) for opt_name in option_names]))
if number_of_options == 0:
    options.normal = True
    if not options.dynamically_loadable:
        options.opt = True
        number_of_options = 2
    else:
        number_of_options = 1

# Special case for -y --[something] --opt
if number_of_options == 2 and options.opt:
    dyn_opt = True
else:
    dyn_opt = False

# Check for .so creation
if options.dynamically_loadable:
    if number_of_options > 1 and not dyn_opt and not options.cvode_data_clamp:
        print 'You asked for ', options
        parser.error("Only one output type may be specified if creating a dynamic library")
        
    essential_options.append('-y')

# What options should be passed on to PyCml?
pycml_options = filter(lambda a: a.startswith('-'), args)
if not pycml_options:
    pycml_options = default_options
pycml_options.extend(essential_options)
if options.quiet:
    pycml_options.append('--quiet')

# Models to process
models = []
for model in filter(lambda a: not a.startswith('-'), args):
    if os.path.exists(model):
        models.append(os.path.abspath(model))
    else:
        print >>sys.stderr, "Skipping", model, "because it does not exist"

def tidy_up():
    """Clean up temporary file, if created."""
    if options.config_file and config_modified:
        os.remove(config_path)

def do_cmd(cmd, outputs):
    """Print and execute a command.
    
    If the command fails, remove any generated outputs and exit.
    """
    if options.show_outputs:
        for output in outputs:
            print output
    else:
        if not options.quiet:
            print ' '.join(cmd)
        rc = subprocess.call(cmd)
        if rc:
            for output in outputs:
                try:
                    os.remove(output)
                except OSError:
                    pass
            tidy_up()
            sys.exit(rc)

# We should only do validation once
done_validation = False

def add_out_opts(base_options, output_dir, classname, file_base, file_extra=''):
    """Add options specifying output path and class name.
    
    Returns extended options list and list of output file paths.
    """
    global done_validation
    if options.dynamically_loadable:
        filename = file_base + '.cpp'
    else:
        filename = file_base + file_extra + '.cpp'
    cpp_path = os.path.join(output_dir, filename)
    full_options = base_options + ['-c', classname, '-o', cpp_path]
    outputs = [cpp_path, cpp_path[:-3] + 'hpp']
    if done_validation or options.assume_valid:
        full_options.append('--assume-valid')
    done_validation = True # Don't do validation next time
    return (full_options, outputs)

def convert(model, output_dir):
    """The main workhorse function."""
    model_dir = os.path.dirname(model)
    model_base = os.path.basename(model)
    model_base = os.path.splitext(model_base)[0]
    if options.dynamically_loadable:
        # If you try to use both normal and dynamic cells in the same test,
        # strange things happen otherwise!
        name_prefix = 'Dynamic'
    else:
        name_prefix = 'Cell'
    class_name = name_prefix + model_base.replace('-', '_') + "FromCellML"
    if not output_dir:
        output_dir = model_dir

    command_base = [os.path.join(pycml_dir, 'translate.py'), model] + pycml_options

    if options.normal:
        # Basic class
        cmd, outputs = add_out_opts(command_base, output_dir, class_name, model_base)
        do_cmd(cmd, outputs)

    if options.opt and (options.normal or number_of_options == 1):
        # Normal with optimisation
        cmd, outputs = add_out_opts(command_base + ['-p', '-l'], output_dir,
                                    class_name + 'Opt', model_base, 'Opt')
        do_cmd(cmd, outputs)
    
    maple_output = os.path.splitext(model)[0] + '.out'
    maple_options = []
    if os.path.exists(maple_output):
        maple_options.extend(['-j', maple_output])

    if options.cvode:
        # For use with CVODE
        if not dyn_opt:
            cmd, outputs = add_out_opts(command_base + ['-t', 'CVODE'] + maple_options,
                                        output_dir,
                                        class_name + 'Cvode', model_base, 'Cvode')
            do_cmd(cmd, outputs)

        if options.opt:
            # With optimisation
            cmd, outputs = add_out_opts(command_base + ['-p', '-l', '-t', 'CVODE'] + maple_options,
                                        output_dir,
                                        class_name + 'CvodeOpt',
                                        model_base, 'CvodeOpt')
            do_cmd(cmd, outputs)
            
    if options.cvode_data_clamp:
        # For use with CVODE and a data clamp
        if not dyn_opt:
            cmd, outputs = add_out_opts(command_base + ['--use-data-clamp', '-t', 'CVODE'] + maple_options,
                                        output_dir,
                                        class_name + 'CvodeDataClamp', model_base, 'CvodeDataClamp')
            do_cmd(cmd, outputs)

        if options.opt:
            # With optimisation
            cmd, outputs = add_out_opts(command_base + ['--use-data-clamp', '-p', '-l', '-t', 'CVODE'] + maple_options,
                                        output_dir,
                                        class_name + 'CvodeDataClampOpt',
                                        model_base, 'CvodeDataClampOpt')
            do_cmd(cmd, outputs)
    
    if options.backward_euler:
        assert maple_options, "A Maple output file is required for Backward Euler code generation"
        be_opts = maple_options + ['-p', '-l', '--backward-euler']
        cmd, outputs = add_out_opts(command_base + be_opts,
                                    output_dir,
                                    class_name + 'BackwardEuler',
                                    model_base, 'BackwardEuler')
        do_cmd(cmd, outputs)
    
    rush_larsen_variants = {'--rush-larsen': 'RushLarsen', '--grl1': 'GRL1', '--grl2': 'GRL2'}
    for rl_opt, rl_name in rush_larsen_variants.iteritems():
        if getattr(options, arg2name(rl_opt), False):
            opts = [rl_opt]
            if rl_opt.startswith('--grl'):
                opts.extend(maple_options)
            if not dyn_opt:
                cmd, outputs = add_out_opts(command_base + opts, output_dir,
                                            class_name + rl_name, model_base, rl_name)
                do_cmd(cmd, outputs)
                
            if options.opt:
                opts.extend(['-p', '-l'])
                cmd, outputs = add_out_opts(command_base + opts, output_dir,
                                            class_name + rl_name + 'Opt',
                                            model_base, rl_name + 'Opt')
                do_cmd(cmd, outputs)


for model in models:
    convert(model, options.output_dir)

tidy_up()
