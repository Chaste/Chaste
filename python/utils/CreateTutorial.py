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


"""Script to convert a tutorial test file into wiki code

Notes: 
1. The script starts converting everything after the first #define line
2. The script stops converting after the final #endif.  It is assumed that there
   will be a #if line immediately after the first #define, and a count is kept of
   #if and #endif lines seen in order to allow matching #if/endif blocks within
   the converted section.
3. All C-style block comments '/*' to '*/' are converted to wiki text.
4. All other lines, including C++ style comments '//', are kept as code lines
5. In C-block comments (ie wiki text), whitespace is removed. Bulleted lists 
   will work but nested bulleted lists won't.
6. To print an empty line in the wiki page, just leave a blank line between paragraphs
   in the comment, as for wiki pages.  The older style of writing EMPTYLINE in a
   paragraph by itself also works, but is deprecated.
7. Lines inside a block comment which start with a '*', i.e.
     /* my comment is
      * two lines long */
   are ok, the initial '*' is removed.
"""


import optparse
import os
import sys

# This had better match GenerateHowTo.py!
HOWTO_TAG = "HOW_TO_TAG"

def ConvertFileToWikiText(fileobj, filepath):
    """Convert a single tutorial source file to wiki markup, returning the page source."""
    # "State machine" state variables
    parsing = False
    ST_NONE, ST_TEXT, ST_CODE, ST_HOWTO = 0, 1, 2, 3 # Status codes
    status = ST_NONE
    in_list = False
    ifdefs_seen = 0
    # Configuration based on file type
    if filepath.endswith('.py'):
        StartsComment = lambda stripped_line: stripped_line.startswith('##')
        IsStillComment = lambda stripped_line: stripped_line.startswith('#')
        EndsComment = lambda stripped_line: not stripped_line.startswith('#')
        CleanEndComment = lambda stripped_line: stripped_line
        end_can_be_start = False # Whether EndsComment and StartsComment can apply to the same line
    else:
        StartsComment = lambda stripped_line: stripped_line.startswith('/*') and not stripped_line.startswith('/**') # Ignore Doxygen comments
        IsStillComment = lambda stripped_line: stripped_line.startswith('*')
        EndsComment = lambda stripped_line: stripped_line.endswith('*/') or stripped_line == '/'
        CleanEndComment = lambda stripped_line: stripped_line[:-2].strip()
        end_can_be_start = True
    code_block_opener = CodeBlockOpener(filepath)
    # Output
    output = []
    code_store = [] # A store of every code-line, to print at the end
    line_store = [] # A copy of the input lines, in case we don't find anything wiki-like
    last_line = ''
    
    for line in fileobj:
        line = line.rstrip() # Note: also removes '\n'
        line_store.append(line)
        # Remove all whitespace and save to a new string.
        # We don't remove the initial whitespace as it will be needed for code lines.
        stripped_line = line.strip()
        
        # We stop processing input after an #endif matching the initial include guard
        if stripped_line.startswith('#endif'):
            assert ifdefs_seen > 0, "#endif seen before #if"
            ifdefs_seen -= 1
            if ifdefs_seen == 0:
                if status is ST_CODE:
                    # close code block
                    output.append('}}}\n')
                parsing = False
    
        # If in Parsing mode
        if parsing:
            if status in [ST_TEXT, ST_HOWTO]:
                # We are still in a comment line, so strip it
                line = stripped_line
            
            # Check if the line is a new text line
            comment_started = False
            if StartsComment(stripped_line):
                comment_started = True
                # remove all whitespace and the '/*'
                stripped_line = line = stripped_line[2:].strip()
                # if the last line was code, close the output code block
                if status is ST_CODE:
                    output.append('}}}\n')
                # set the status as text
                status = ST_TEXT
            elif status in [ST_TEXT, ST_HOWTO] and IsStillComment(stripped_line):
                # we are in a comment, so get rid of whitespace and the initial '*'
                stripped_line = line = line[1:].strip()
            elif status is ST_NONE and len(stripped_line) > 0:
                # Line has content and isn't a comment => it's code
                output.append(code_block_opener)
                status = ST_CODE
            
            # Check if comment ends
            if EndsComment(stripped_line) and (not comment_started or end_can_be_start):
                # If it's not a Doxygen comment, switch state to unknown
                if status in [ST_TEXT, ST_HOWTO]:
                    # get rid of whitespace and '*/'
                    stripped_line = line = CleanEndComment(stripped_line)
                    status = ST_NONE
            
            # Check for (and strip) HOWTO tagging
            if status is ST_TEXT and stripped_line.startswith(HOWTO_TAG):
                status = ST_HOWTO
            if status is ST_HOWTO:
                if not stripped_line:
                    status = ST_TEXT # Blank comment line ends tagging
                else:
                    stripped_line = line = '' # Strip tag content
            
            if status is ST_TEXT and stripped_line and stripped_line[0] == '*':
                # It's a list, so needs some indentation!
                in_list = True
            if status is ST_TEXT and in_list:
                line = ' '+line
            if in_list and (len(stripped_line) == 0 or status is not ST_TEXT):
                in_list = False
    
            # If the line is a comment just saying 'EMPTYLINE', we'll print a blank line
            if stripped_line == 'EMPTYLINE':
                stripped_line = line = ''
            # We print the line unless we'd get 2 empty lines
            if len(stripped_line) > 0 or len(last_line) > 0:
                output.append(line+'\n')
            
            # If the line is a code line we store it,
            # unless there would be two consecutive empty lines
            if status is ST_CODE:
                if len(stripped_line) > 0 or len(code_store[-1].strip()) > 0:
                    code_store.append(line)
    
        # We start processing lines AFTER the first #define..
        if stripped_line.startswith('#define'):
            parsing = True
        if stripped_line.startswith('#if'):
            ifdefs_seen += 1
        last_line = stripped_line
    
    if not output:
        # It's probably not C++ or Python, so let's include it all just as raw 'code'
        code_store = line_store

    return ''.join(output), code_store

def CodeBlockOpener(file_name):
    """Return the opener string for a Trac wiki code block with syntax highlighting based on file extension."""
    ext = os.path.splitext(file_name)[1]
    highlight_code = {'.hpp': '#!cpp\n', '.cpp': '#!cpp\n', '.py': '#!python\n', '.sh': '#!sh\n'}.get(ext, '')
    return '{{{\n' + highlight_code

def AddCodeOutput(file_name, code, output):
    output.append('\n\n== File name `%s` ==\n\n' % file_name)
    output.append(CodeBlockOpener(file_name))
    output.append('\n'.join(code))
    output.append('\n}}}\n\n')

def ConvertTutorialToWikiText(test_file_path, test_file, other_files, revision=''):
    """Convert a tutorial, possibly comprised of multiple files, to wiki markup.
    
    test_file is the content of the tutorial test .hpp file, as an object which will
    return each line in turn when iterated.
    other_files is a list of subsidiary files, which may be empty.  Each entry should
    be a pair, the first item of which is the basename of the file, and the second is
    the contents as for test_file.
    """
    if revision:
        revision = str(revision.strip())

        if test_file_path.strip().startswith('projects/'):
            assert revision.isdigit(), 'Expected svn-style integer revision'
            revision = 'at revision r{}'.format(revision)
        else:
            assert len(revision) == 40,  'Expected 40-digit git commit hash'
            revision = 'at revision [changeset:{}/git_repo]'.format(revision[0:12])  # abbreviate to consistent length

    output = []

    # Header
    output.append('This tutorial is automatically generated from the file {} {}.\n'.format(test_file_path, revision))
    output.append('Note that the code is given in full at the bottom of the page.\n\n\n')

    # Convert each file in turn
    test_output, test_code = ConvertFileToWikiText(test_file, test_file_path)
    output.append(test_output)
    other_code = {}
    for other_file in other_files:
        file_output, file_code = ConvertFileToWikiText(other_file[1], other_file[0])
        if file_output:
            output.append('\n\n= Extra file %s =\n' % other_file[0])
            output.append(file_output)
        if file_code:
            other_code[other_file[0]] = file_code
    # Now output the C++ code for all files
    output.append('\n\n= Code =\nThe full code is given below\n')
    AddCodeOutput(os.path.basename(test_file_path), test_code, output)
    for filename, code in other_code.iteritems():
        AddCodeOutput(filename, code, output)
    return ''.join(output)


def ParseOptions():
    usage = "usage: %prog [options] <test_file>|- <output_file>|-"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-r', '--revision', default='',
                      help="Revision")
    parser.add_option('-f', '--real-file-path',
                      help="The real path of the test file, if it is being piped in")
    (options, args) = parser.parse_args()
    
    if len(args) != 2:
        parser.error("You must specify input and output files")
    
    return options, args


if __name__ == '__main__':
    options, args = ParseOptions()
    
    test_file = args[0]
    out_file_name = args[1]
    
    if options.real_file_path:
        real_file_path = options.real_file_path
    else:
        real_file_path = test_file

    # We're interested in the path relative to the root of the Chaste directory
    chaste_build_dir_root = os.path.abspath(os.path.join(os.getcwd(), '..', '..', '..'))
    real_file_path = os.path.relpath(test_file, chaste_build_dir_root)
    
    if test_file == '-':
        # Read from stdin (pipe mode)
        in_file = sys.stdin
    elif test_file[-3:] not in ['hpp', '.py'] or os.path.basename(test_file)[:4] != 'Test':
        print >>sys.stderr, "Syntax error:", test_file, "does not appear to be a test file"
        sys.exit(1)
    else:
        in_file = open(test_file)
    
    if out_file_name == '-':
        # Write to stdout (pipe mode)
        out_file = sys.stdout
    else:
        out_file = open(out_file_name, 'w')
    
    # Do the conversion
    out_file.write(ConvertTutorialToWikiText(real_file_path, in_file, [], options.revision))
    
    # Close files
    if in_file is not sys.stdin:
        in_file.close()
    
    if out_file is not sys.stdout:
        out_file.close()

