"""Copyright (c) 2005-2025, University of Oxford.
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

from pathlib import Path
from typing import List, Tuple
from subprocess import check_output, CalledProcessError

import re

CHASTE_SRC_DIR = Path(__file__).resolve().parent.parent.parent
HEART_DIR = CHASTE_SRC_DIR / 'heart'
USER_TUTORIAL_DIR = CHASTE_SRC_DIR / 'tutorials' / 'UserTutorials'
HOWTO_TAG = "HOW_TO_TAG"
# Status codes for use during parsing
ST_NONE, ST_TEXT, ST_CODE, ST_HOWTO = 0, 1, 2, 3


def get_list_of_user_tutorial_files() -> List[Path]:
    """
    Scans the CHASTE_SRC_DIR directory recursively to find all files with a name starting with 'Test'
    and ending with 'Tutorial.hpp', and returns a list of their resolved paths.

    :return: List of resolved Path objects of test files.
    """
    return [file.resolve() for file in CHASTE_SRC_DIR.rglob('Test*Tutorial.hpp')]


def get_last_revision_hash(abs_file_path: Path) -> str:
    """
    Retrieves the hash of the last revision in which the specified file was changed.

    :param abs_file_path: An absolute path to the file within the Git repository.
    :return: The hash of the last commit that changed the file.
    """
    try:
        rel_file_path = abs_file_path.relative_to(CHASTE_SRC_DIR)
        rev_hash = check_output(
            ['git', '-C', str(CHASTE_SRC_DIR), 'log', '-n', '1',
             '--pretty=format:%H', '--', str(rel_file_path)],
            encoding='utf-8'
        ).strip()
        return rev_hash
    except CalledProcessError as e:
        print(f'Error getting last revision for {abs_file_path}.')
        print(f'Output:\n{e.output}')
        print(f'Error:\n{e.stderr}')
        exit(e.returncode)


def get_output_file_path(abs_file_path: Path) -> Path:
    """
    Generates a markdown file Path in USER_TUTORIAL_DIR from a given file's stem,
    assuming 'Test' prefix and 'Tutorial' suffix in the filename.

    :param abs_file_path: The absolute Path of the tutorial file.
    :return: A new Path object with a modified filename and '.md' extension.
    :raises AssertionError: If the filename doesn't match the expected pattern.
    """
    file_name = abs_file_path.stem
    assert file_name.startswith('Test') and file_name.endswith(
        'Tutorial'), 'File does not appear to be a tutorial'
    return USER_TUTORIAL_DIR / f'{file_name[len("Test"):-(len("Tutorial"))]}.md'


def get_title_from_file(abs_file_path: Path) -> str:
    file_name = abs_file_path.stem
    assert file_name.startswith('Test') and file_name.endswith(
        'Tutorial'), 'File does not appear to be a tutorial'
    reduced_file_name = file_name[len("Test"):-(len("Tutorial"))]

    # Split CamelCase string into list of words
    words = re.findall(r'(\d+[A-Za-z]?|[A-Z][a-z]*|[A-Z]+(?=[A-Z]|$))', reduced_file_name)

    # Join words with space and return the result
    return ' '.join(words)


def get_summary_from_file(abs_file_path: Path) -> str:
    """
    Opens a file and searches for the first occurrence of a pattern that identifies
    a Markdown heading. Returns the text of the first Markdown heading without the '##' prefix.

    :param abs_file_path: The Path object of the file to be searched.
    :return: A string containing the text of the first Markdown heading, or an empty string if not found.
    """
    # Compile the regular expression
    pattern = re.compile(r'.*##\s+(.+)')

    # Read the file and search for the pattern
    with abs_file_path.open('r', encoding='utf-8') as file:
        for line in file:
            match = pattern.match(line)
            if match:
                return match.group(1)

    return ''


def get_revision_string_from_file(abs_file_path: Path) -> str:
    """
    Forms a string at the top of each file that looks like:

    This tutorial is automatically generated from {file_link} at revision {revision}.
    Note that the code is given in full at the bottom of the page.

    :param abs_file_path: The Path object of the file.
    :return: A string containing the text indicated above.
    """
    relative_file_path = abs_file_path.relative_to(CHASTE_SRC_DIR)
    file_link = f'[{abs_file_path.name}](https://github.com/Chaste/Chaste/blob/develop/{relative_file_path})'

    revision = get_last_revision_hash(abs_file_path)
    revision_link = f'[{revision[0:12]}](https://github.com/Chaste/Chaste/commit/{revision})'

    return f'This tutorial is automatically generated from {file_link} at revision {revision_link}. Note that the code is given in full at the bottom of the page.'


def write_tutorial_markdown(abs_file_path: Path) -> None:
    """
    Convert a tutorial to hugo markdown.

    :param abs_file_path: The Path object of the tutorial file.
    """
    revision_string = get_revision_string_from_file(abs_file_path)

    output = []

    header = f"""---
title : "{get_title_from_file(abs_file_path)}"
summary: "{get_summary_from_file(abs_file_path)}"
draft: false
images: []
toc: true
---
"""
    output.append(header)
    output.append(revision_string)

    test_output, test_code = convert_file_to_markdown(abs_file_path)
    output.append(test_output)

    # Now output the C++ code for all files
    output.append(f'\n\n## Full code\n\n{format_code_lines(test_code)}')

    full_markdown_content = ''.join(output)

    # Postprocess to remove empty lines at the end of code blocks
    full_markdown_content = full_markdown_content.replace(
        '\n\n```\n', '\n```\n')

    # Postprocess to remove more than one empty line above a code block
    full_markdown_content = full_markdown_content.replace(
        '\n\n\n```cpp\n', '\n\n```cpp\n')

    # Postprocess to remove any cases of 3 blank lines
    full_markdown_content = full_markdown_content.replace('\n\n\n\n', '\n\n')

    # Postprocess to remove any cases of 2 blank lines
    full_markdown_content = full_markdown_content.replace('\n\n\n', '\n\n')

    with open(get_output_file_path(abs_file_path), 'w', encoding='utf-8') as file:
        file.write(full_markdown_content)


def format_code_lines(code_lines: List[str]) -> str:
    """
    Formats a list of code lines into a Markdown code block with C++ syntax highlighting.

    :param code_lines: A list of strings, where each string is a line of code.
    :return: A string representing the entire code block, enclosed within Markdown code block
             fencing for C++ syntax highlighting.
    """
    joined_code = '\n'.join(code_lines)
    return f'```cpp\n{joined_code}\n```\n'


def starts_comment(_stripped_line):
    """
    Determine if the stripped line is the start of a comment block.

    :param _stripped_line: The line to be checked.
    :return: True if the line starts with '/*' but not with '/**' (ignoring Doxygen comments), False otherwise.
    """
    return _stripped_line.startswith('/*') and not _stripped_line.startswith('/**')


def is_still_comment(_stripped_line):
    """
    Check if the stripped line is part of a comment block.

    :param _stripped_line: The line to be checked.
    :return: True if the line starts with '*', indicating it's part of a comment block, False otherwise.
    """
    return _stripped_line.startswith('*')


def ends_comment(_stripped_line):
    """
    Determine if the stripped line is the end of a comment block.

    :param _stripped_line: The line to be checked.
    :return: True if the line ends with '*/' or is '/', indicating the end of a comment block, False otherwise.
    """
    return _stripped_line.endswith('*/') or _stripped_line == '/'


def clean_end_comment(_stripped_line):
    """
    Remove the closing characters of a comment block from the stripped line.

    :param _stripped_line: The line to be cleaned.
    :return: The cleaned line with the comment block end removed, if present.
    """
    return _stripped_line[:-2].strip()


def convert_file_to_markdown(abs_file_path: Path) -> Tuple[str, List[str]]:
    """
    Convert a single tutorial source file to wiki markup, returning the page source.

    :param abs_file_path: The Path object of the tutorial file.
    :return: a string representing the tutorial text, and a list of strings representing the full code of the tutorial
    """
    # "State machine" state variables
    parsing = False

    status = ST_NONE
    in_list = False
    ifdefs_seen = 0

    end_can_be_start = True
    code_block_opener = '\n```cpp\n'

    # Output
    output = []
    code_store = []  # A store of every code-line, to print at the end
    line_store = []  # A copy of the input lines, in case we don't find anything wiki-like
    last_line = ''

    with abs_file_path.open('r', encoding='utf-8') as fileobj:

        for line in fileobj:
            line = line.rstrip()  # Note: also removes '\n'
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
                        output.append('```\n\n')
                    parsing = False

            # If in Parsing mode
            if parsing:
                if status in [ST_TEXT, ST_HOWTO]:
                    # We are still in a comment line, so strip it
                    line = stripped_line

                # Check if the line is a new text line
                comment_started = False
                if starts_comment(stripped_line):
                    comment_started = True
                    # remove all whitespace and the '/*'
                    stripped_line = line = stripped_line[2:].strip()
                    # if the last line was code, close the output code block
                    if status is ST_CODE:
                        output.append('```\n\n')
                    # set the status as text
                    status = ST_TEXT
                elif status in [ST_TEXT, ST_HOWTO] and is_still_comment(stripped_line):
                    # we are in a comment, so get rid of whitespace and the initial '*'
                    stripped_line = line = line[1:].strip()
                elif status is ST_NONE and len(stripped_line) > 0:
                    # Line has content and isn't a comment => it's code
                    output.append(code_block_opener)
                    status = ST_CODE

                # Check if comment ends
                if ends_comment(stripped_line) and (not comment_started or end_can_be_start):
                    # If it's not a Doxygen comment, switch state to unknown
                    if status in [ST_TEXT, ST_HOWTO]:
                        # get rid of whitespace and '*/'
                        stripped_line = line = clean_end_comment(stripped_line)
                        status = ST_NONE

                # Check for (and strip) HOWTO tagging
                if status is ST_TEXT and stripped_line.startswith(HOWTO_TAG):
                    status = ST_HOWTO
                if status is ST_HOWTO:
                    if not stripped_line:
                        status = ST_TEXT  # Blank comment line ends tagging
                    else:
                        stripped_line = line = ''  # Strip tag content

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


if __name__ == '__main__':

    USER_TUTORIAL_DIR.mkdir(parents=True, exist_ok=True)
    if any(USER_TUTORIAL_DIR.iterdir()):
        print(f'Warning: {USER_TUTORIAL_DIR} is not empty')

    for user_tutorial_file in get_list_of_user_tutorial_files():
        write_tutorial_markdown(user_tutorial_file)
