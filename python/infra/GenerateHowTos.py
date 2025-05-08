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

from collections import Counter
from collections import defaultdict
from collections import namedtuple
from pathlib import Path


HowTo = namedtuple('HowTo', ['primary', 'secondary', 'line_number', 'test_name', 'github_link', 'how_to_text'])

CHASTE_SRC_DIR = Path(__file__).resolve().parent.parent.parent

HOW_TO_TAG = 'HOW_TO_TAG'

HOW_TO_FILE = CHASTE_SRC_DIR / 'chaste_howto_webpage.md'

PAGE_HEADER = """---
title: "How-to index"
description: "Chaste how-tos, automatically generated from the test files"
draft: false
images: []
toc: true
layout: "single"
---

As well as ensuring that Chaste functions as expected, the many Chaste tests can be effective for learning how to use Chaste.
There are often aspects of Chaste's capabilities that do not warrant a full user tutorial, but are still worth noting.
Equally, users often ask "how do I do...", and in many cases there is already a test that does something similar and can be used as a basis.

This page contains a (probably partial) index of links to useful code.
Some of the links are to code in tutorials.
The others are links to test code.
For the latter, note that since tests are written primarily for checking functionality, they are not commented to the same degree as user tutorials.
Nevertheless, they should be reasonably readable and useful.

Note that this page is generated automatically based on tags in the Chaste code.
Do not edit it manually, as your changes will be overwritten!
"""

PAGE_FOOTER = """
## Creating new how-tos

Any C-style block comment found in a test file that looks like the following will be given an entry in the index:

 * HOW_TO_TAG Section/Subsection
 * Short description (ideally one-liner)

The short description may run over multiple lines; it will be considered to end either at the end of the comment, or at a blank comment line, whichever comes first.
This index is automatically updated by a GitHub Actions workflow triggered by any new commit pushed to the `develop` branch.
"""


def get_list_of_test_files() -> list[Path]:
    """
    Scans the CHASTE_SRC_DIR directory recursively to find all files with a name starting with 'Test'
    and ending with '.hpp', and returns a list of their resolved paths.

    :return: List of resolved Path objects of test files.
    """
    return [file.resolve() for file in CHASTE_SRC_DIR.rglob('Test*.hpp')]


def stop_condition(line_to_check: str) -> bool:
    """
    Determines if a given line should signify a stop condition. A line is considered a stop condition
    if it is empty, a single asterisk, or the end of a block comment.

    :param line_to_check: The line of text to check for the stop condition.
    :return: True if the line signifies a stop condition, False otherwise.
    """
    line_to_check = line_to_check.strip()
    if line_to_check == '' or line_to_check == '*' or line_to_check == '*/':
        return True
    return False


def process_line(line_to_process: str) -> str:
    """
    Processes a line of text by removing asterisks and forward slashes and stripping whitespace from the ends.
    Adds a single space character at the beginning of the processed line.

    :param line_to_process: The line of text to be processed.
    :return: The processed line with leading space.
    """
    return ' ' + line_to_process.strip().lstrip('*/').lstrip()


def parse_how_to_line(line_to_parse: str) -> tuple[str, str]:
    """
    Parses a line of text to extract components and categories after a specified HOW_TO_TAG.
    Splits the line at the HOW_TO_TAG and then further splits the result to get the component
    and category, which are returned as a tuple.

    :param line_to_parse: The line of text containing the HOW_TO_TAG to be parsed.
    :return: A tuple containing the component and category.
    """
    _, _, after = line_to_parse.partition(HOW_TO_TAG)
    _component, _, _category = after.strip().partition('/')
    return _component, _category


def form_github_link(file_path: Path, line_num: int) -> str:
    """
    Constructs a GitHub URL to a specific line number in a file within the Chaste repository.
    The file path is made relative to the CHASTE_SRC_DIR, and the line number is appended to the URL.

    :param file_path: The path object of the file.
    :param line_num: The line number to link to.
    :return: A string containing the complete GitHub URL to the file and line number.
    """
    return f'https://github.com/Chaste/Chaste/blob/develop/{file_path.relative_to(CHASTE_SRC_DIR)}#L{line_num}'


def get_all_how_tos() -> list[HowTo]:
    """
    Retrieves a list of HowTo namedtuples from all test files within the CHASTE_SRC_DIR directory.
    Each HowTo contains a primary component, secondary category, line number, test name, GitHub link,
    and the text of the how-to. It searches for lines containing the HOW_TO_TAG, then continues to read
    the following lines until it encounters a stop condition to form the how-to text. If a how-to section
    exceeds 5 lines, a warning is printed indicating that the how-to text might be too long, which could
    mean that a stop line is missing.

    :return: A list of HowTo namedtuples containing extracted information from test files.
    """
    _list_of_how_tos = []

    for test_file in get_list_of_test_files():

        with open(test_file, 'r') as f:
            inside_how_to = False  # whether the current line is part of a how-to block
            how_to_text_lines = 0  # counter for the number of lines in the how-to text

            for line_number, line in enumerate(f, start=1):

                # If the line contains a HOW_TO_TAG then we are inside a how-to, and we parse the category
                if HOW_TO_TAG in line:
                    inside_how_to = True
                    how_to_text = ''
                    how_to_text_lines = 0  # Reset counter for each how-to block
                    primary, secondary = parse_how_to_line(line)
                    how_to_line_number = line_number
                    continue  # Skip the rest of the loop and move to the next line

                if inside_how_to:

                    # If we meet the stop condition, then we have finished parsing the how-to block and we process it
                    if stop_condition(line):
                        inside_how_to = False  # We've reached the end of the how-to block
                        if how_to_text_lines > 5:
                            print(f'Found {how_to_text_lines} lines in the how-to on line {line_number} of {test_file}:'
                                  ' is there an empty line signifying the end of the how-to?')
                            exit(1)
                        # Create the namedtuple HowTo instance
                        how_to = HowTo(
                            primary=primary, secondary=secondary, line_number=how_to_line_number,
                            test_name=test_file.name, github_link=form_github_link(test_file, how_to_line_number),
                            how_to_text=how_to_text.strip()
                        )
                        _list_of_how_tos.append(how_to)
                    else:
                        # Increment the counter as we add lines to the how-to text
                        how_to_text_lines += 1
                        how_to_text += process_line(line)

    return _list_of_how_tos


def summarise_how_tos(how_tos: list[HowTo]) -> None:
    """
    Summarises the provided how-tos by counting how many times each primary category occurs.
    Prints out the total number of how-tos and a summary of counts per primary category.

    :param how_tos: A list of HowTo namedtuples.
    :return: None. This function prints the summary directly to the console.
    """
    primary_counts = Counter([how_to.primary for how_to in how_tos])

    print(f'Found {len(how_tos)} how-tos, in the following primary categories:')
    print(primary_counts)


def howto_to_markdown_snippet(how_to: HowTo) -> str:
    """
    Converts a HowTo namedtuple into a markdown formatted string. The function formats the how-to text
    as a list item and appends a link to the specific line in the test file on GitHub.

    :param how_to: A HowTo namedtuple containing the how-to text, line number, test file name, and GitHub link.
    :return: A string formatted in markdown, with the how-to text and a link to the corresponding line in the GitHub repository.
    """
    return f'- {how_to.how_to_text}\n  - [line {how_to.line_number} of {how_to.test_name}]({how_to.github_link})'


def generate_markdown_file(how_tos: list[HowTo]) -> None:
    """
    Generates a markdown file with the provided how-tos. The how-tos are grouped by primary and secondary
    categories with headers for each section and subsection, respectively.

    :param how_tos: A list of HowTo namedtuples.
    :return: None. Writes the output directly to the file specified by HOW_TO_FILE.
    """
    print(f'Writing markdown to: {HOW_TO_FILE}')

    # Create a nested dictionary where keys are primary categories
    # and values are dictionaries with secondary categories as keys
    how_to_dict = defaultdict(lambda: defaultdict(list))
    for how_to in how_tos:
        how_to_dict[how_to.primary][how_to.secondary].append(how_to)

    with open(HOW_TO_FILE, 'w') as f:
        f.write(PAGE_HEADER)

        # Iterate over sorted primary categories
        for primary in sorted(how_to_dict):
            f.write(f'\n## {primary}\n')

            # Iterate over sorted secondary categories within each primary category
            for secondary in sorted(how_to_dict[primary]):
                if secondary != "":
                    f.write(f'\n### {secondary}\n')

                # Write each how-to markdown snippet
                for how_to in sorted(how_to_dict[primary][secondary], key=lambda x: (x.test_name, x.line_number)):
                    f.write(howto_to_markdown_snippet(how_to) + '\n')

        f.write(PAGE_FOOTER)


if __name__ == '__main__':

    list_of_how_tos = get_all_how_tos()
    summarise_how_tos(list_of_how_tos)
    generate_markdown_file(list_of_how_tos)
