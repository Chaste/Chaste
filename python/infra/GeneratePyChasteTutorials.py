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

import os
import re
import subprocess
import warnings
from pathlib import Path
from subprocess import CalledProcessError, check_output
from typing import List, Tuple

import nbformat
from nbformat import NotebookNode
from nbformat.v4 import new_code_cell, new_markdown_cell, new_notebook

CHASTE_SRC_DIR = Path(__file__).resolve().parent.parent.parent
PYCHASTE_SRC_DIR = CHASTE_SRC_DIR / "pychaste"
PYCHASTE_TUTORIAL_DIR = PYCHASTE_SRC_DIR / "src" / "py" / "doc" / "tutorial"
HOWTO_TAG = "HOW_TO_TAG"

# Status codes for use during parsing
ST_NONE, ST_TEXT, ST_CODE, ST_HOWTO = 0, 1, 2, 3

# Placeholder strings for notebook tutorials
JUPYTER_SETUP = "chaste.cell_based.SetupNotebookTest() # Set up the test"
JUPYTER_TEARDOWN = "chaste.cell_based.TearDownNotebookTest() # Tear down the test"
JUPYTER_SHOW_FIRST = "nb_manager = chaste.visualization.JupyterNotebookManager()"
JUPYTER_SHOW = "nb_manager.vtk_show(scene, height=600)"


def clean_end_comment(_stripped_line):
    """
    Remove the closing characters of a comment block from the stripped line.

    :param _stripped_line: The line to be cleaned.
    :return: The cleaned line with the comment block end removed, if present.
    """
    return _stripped_line


def ends_comment(_stripped_line):
    """
    Determine if the stripped line is the end of a comment block.

    :param _stripped_line: The line to be checked.
    :return: True if the line ends with '*/' or is '/', indicating the end of a comment block, False otherwise.
    """
    return not _stripped_line.startswith("#")


def is_still_comment(_stripped_line):
    """
    Check if the stripped line is part of a comment block.

    :param _stripped_line: The line to be checked.
    :return: True if the line starts with '*', indicating it's part of a comment block, False otherwise.
    """
    return _stripped_line.startswith("#")


def starts_comment(_stripped_line):
    """
    Determine if the stripped line is the start of a comment block.

    :param _stripped_line: The line to be checked.
    :return: True if the line starts with '/*' but not with '/**' (ignoring Doxygen comments), False otherwise.
    """
    return _stripped_line.startswith("##")


def format_code_lines(code_lines: List[str]) -> str:
    """
    Formats a list of code lines into a Markdown code block with Python syntax highlighting.

    :param code_lines: A list of strings, where each string is a line of code.
    :return: A string representing the entire code block, enclosed within Markdown code block
             fencing for Python syntax highlighting.
    """
    joined_code = "\n".join(code_lines)
    return f"```python\n{joined_code}\n```\n"


def get_last_revision_hash(abs_file_path: Path) -> str:
    """
    Retrieves the hash of the last revision in which the specified file was changed.

    :param abs_file_path: An absolute path to the file within the Git repository.
    :return: The hash of the last commit that changed the file.
    """
    try:
        rel_file_path = abs_file_path.relative_to(CHASTE_SRC_DIR)
        rev_hash = check_output(
            [
                "git",
                "-C",
                str(CHASTE_SRC_DIR),
                "log",
                "-n",
                "1",
                "--pretty=format:%H",
                "--",
                str(rel_file_path),
            ],
            encoding="utf-8",
        ).strip()
        return rev_hash

    except CalledProcessError as e:
        print(f"Error getting last revision for {abs_file_path}.")
        print(f"Output:\n{e.output}")
        print(f"Error:\n{e.stderr}")
        exit(e.returncode)


def get_list_of_tutorial_files() -> List[Path]:
    """
    Scans the PYCHASTE_SRC_DIR directory recursively to find all files with a name
    matching 'TestPy*Tutorial.py', and returns a list of their resolved paths.

    :return: List of resolved Path objects of tutorial source files.
    """
    return [file.resolve() for file in PYCHASTE_SRC_DIR.rglob("TestPy*Tutorial.py")]


def get_output_file_name(abs_file_path: Path) -> str:
    """
    Extracts the tutorial name from the source file path e.g. TestPyCellSortingTutorial.py -> CellSorting.

    :param abs_file_path: The Path object of the tutorial source file.
    :return: The tutorial name extracted from the source file path.
    """
    return abs_file_path.stem[len("TestPy") : -(len("Tutorial"))]


def get_revision_string_from_file(abs_file_path: Path) -> str:
    """
    Forms a string at the top of each file that looks like:

    This tutorial is automatically generated from {file_link} at revision {revision}.
    Note that the code is given in full at the bottom of the page.

    :param abs_file_path: The Path object of the file.
    :return: A string containing the text indicated above.
    """
    relative_file_path = abs_file_path.relative_to(CHASTE_SRC_DIR)
    file_link = f"[{abs_file_path.name}](https://github.com/Chaste/Chaste/blob/develop/{relative_file_path})"

    revision = get_last_revision_hash(abs_file_path)
    revision_link = (
        f"[{revision[0:12]}](https://github.com/Chaste/Chaste/commit/{revision})"
    )

    return f"This tutorial is automatically generated from {file_link} at revision {revision_link}."


def get_title_from_file(abs_file_path: Path) -> str:
    """
    Extracts the tutorial title from the source file path e.g. TestPyCellSortingTutorial.py -> Cell Sorting.

    :param abs_file_path: The Path object of the tutorial source file.
    :return: The tutorial title extracted from the source file path.
    """
    tutorial_name = get_output_file_name(abs_file_path)

    # Split CamelCase string into list of words
    words = re.findall(r"(\d+[A-Za-z]?|[A-Z][a-z]*|[A-Z]+(?=[A-Z]|$))", tutorial_name)

    # Join words with space and return the result
    return " ".join(words)


def convert_tutorial_file(abs_file_path: Path) -> Tuple[str, NotebookNode, List[str]]:
    """
    Convert a single tutorial source file to markdown format.

    Notes:
    1. Starts converting everything after the first #define line
    2. Stops converting after the final #endif.  It is assumed that there
       will be a #if line immediately after the first #define, and a count is kept of
       #if and #endif lines seen in order to allow matching #if/endif blocks within
       the converted section.
    3. All C-style block comments '/*' to '*/' are converted to markdown.
    4. All other lines, including C++ style comments '//', are kept as code lines
    5. In C-block comments (ie markdown), whitespace is removed. Bulleted lists
       will work but nested bulleted lists won't.
    6. To print an empty line in the markdown page, just leave a blank line between paragraphs
       in the comment, as for markdown pages.  The older style of writing EMPTYLINE in a
       paragraph by itself also works, but is deprecated.
    7. Lines inside a block comment which start with a '*', i.e.
       /* my comment is
       * two lines long */
       are ok, the initial '*' is removed.

    :param abs_file_path: The Path object of the tutorial file.
    :return: a string representing the tutorial markdown text, and a list of strings representing the full code of the tutorial.
    """
    # "State machine" state variables
    parsing = False
    status = ST_NONE
    in_list = False
    ifdefs_seen = 0

    # Whether ends_comment and starts_comment can apply to the same line
    end_can_be_start = False

    code_block_opener = "\n```python\n"

    # Output
    markdown = []

    notebook = new_notebook()  # Jupyter notebook
    notebook["cells"] = []
    notebook_cells = []  # [string, is_code]
    last_notebook_cell = ""

    code_store = []  # Code lines

    line_store = []  # Input lines
    last_line = ""

    with abs_file_path.open("r", encoding="utf-8") as fileobj:

        for line in fileobj:
            line = line.rstrip()  # Note: also removes '\n'
            line_store.append(line)
            # Remove all whitespace and save to a new string.
            # We don't remove the initial whitespace as it will be needed for code lines.
            stripped_line = line.strip()
            squashed_line = "".join(stripped_line.split())

            # We stop processing input after an #endif matching the initial include guard
            if squashed_line.startswith("#endif"):
                assert ifdefs_seen > 0, "#endif seen before #if"
                ifdefs_seen -= 1
                if ifdefs_seen == 0:
                    if status is ST_CODE:
                        # close code block
                        markdown.append("```\n\n")
                        notebook_cells.append([last_notebook_cell, True])
                        last_notebook_cell = ""
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
                        markdown.append("```\n")
                        notebook_cells.append([last_notebook_cell, True])
                        last_notebook_cell = ""
                    # set the status as text
                    status = ST_TEXT
                elif status in [ST_TEXT, ST_HOWTO] and is_still_comment(stripped_line):
                    # we are in a comment, so get rid of whitespace and the initial '*'
                    stripped_line = line = line[1:].strip()
                elif status is ST_NONE and len(stripped_line) > 0:
                    # Line has content and isn't a comment => it's code
                    markdown.append(code_block_opener)
                    notebook_cells.append([last_notebook_cell, False])
                    last_notebook_cell = ""
                    status = ST_CODE

                # Check if comment ends
                if ends_comment(stripped_line) and (
                    not comment_started or end_can_be_start
                ):
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
                        stripped_line = line = ""  # Strip tag content

                if status is ST_TEXT and stripped_line and stripped_line[0] == "*":
                    # It's a list, so needs some indentation!
                    in_list = True
                if status is ST_TEXT and in_list:
                    line = " " + line
                if in_list and (len(stripped_line) == 0 or status is not ST_TEXT):
                    in_list = False

                # If the line is a comment just saying 'EMPTYLINE', we'll print a blank line
                if stripped_line == "EMPTYLINE":
                    stripped_line = line = ""
                # We print the line unless we'd get 2 empty lines
                if len(stripped_line) > 0 or len(last_line) > 0:
                    markdown.append(line + "\n")
                    last_notebook_cell += line + "\n"

                # If the line is a code line we store it,
                # unless there would be two consecutive empty lines
                if status is ST_CODE:
                    if len(stripped_line) > 0 or len(code_store[-1].strip()) > 0:
                        code_store.append(line)

            # We start processing lines AFTER the first #define..
            if squashed_line.startswith("#define"):
                parsing = True
            if squashed_line.startswith("#if"):
                ifdefs_seen += 1
            last_line = stripped_line

        if not markdown:
            # It's probably not C++ or Python, so let's include it all just as raw 'code'
            code_store = line_store

        # Assemble the notebook cells
        for notebook_cell in notebook_cells:
            if notebook_cell[1]:
                # Convert the notebook cell code to a list for easier processing
                cell_code_lines = notebook_cell[0].split("\n")

                output_lines = []
                skip_list = ["unittest", "__main__", "self.assert"]

                for code_line in cell_code_lines:
                    code_line = code_line.rstrip()
                    indentation = len(code_line)
                    code_line = code_line.lstrip()
                    indentation -= len(code_line)

                    # Strip out lines related to unittest and __main__
                    if any(text in code_line for text in skip_list):
                        continue

                    # Strip out class definition lines such as
                    # `class TestPyCellSortingTutorial(AbstractCellBasedTestSuite):`
                    if code_line.startswith("class Test"):
                        continue

                    # Strip out test function definition lines such as
                    # `def test_potts_monolayer_cell_sorting(self):`
                    if code_line.startswith("def test_") and indentation == 4:
                        continue

                    # Replace notebook placeholder strings
                    if "JUPYTER_SETUP" in code_line:
                        code_line = JUPYTER_SETUP

                    elif "JUPYTER_TEARDOWN" in code_line:
                        code_line = JUPYTER_TEARDOWN

                    elif "JUPYTER_SHOW_FIRST" in code_line:
                        code_line = JUPYTER_SHOW_FIRST

                    elif "JUPYTER_SHOW" in code_line:
                        code_line = JUPYTER_SHOW

                    elif "cell_based.VtkSceneModifier_2()" in code_line:
                        code_line = code_line.replace(
                            "cell_based.VtkSceneModifier_2()",
                            "visualization.JupyterSceneModifier_2(nb_manager)",
                        )

                    elif "cell_based.VtkSceneModifier_3()" in code_line:
                        code_line = code_line.replace(
                            "cell_based.VtkSceneModifier_3()",
                            "visualization.JupyterSceneModifier_3(nb_manager)",
                        )

                    # Reduce indentation by 8: 4 for `class Test...` and 4 for `def test_...`
                    indentation = max(indentation - 8, 0)
                    output_lines.append(" " * indentation + code_line)

                # Reassemble the notebook cell into a string, stripping empty lines
                if output_lines:
                    output_string = os.linesep.join(
                        [line for line in output_lines if line.strip()]
                    )
                    if output_string:
                        notebook["cells"].append(new_code_cell(output_string))

            else:
                # Markdown cell
                notebook["cells"].append(new_markdown_cell(notebook_cell[0]))

    return "".join(markdown), notebook, code_store


def write_tutorial(abs_file_path: Path) -> None:
    """
    Convert a tutorial source file to a Hugo markdown file and a Jupyter notebook.

    :param abs_file_path: The Path object of the tutorial source file.
    """
    # Convert the tutorial file to markdown and Jupyter notebook format
    test_markdown, test_notebook, test_code = convert_tutorial_file(abs_file_path)

    # Process the markdown file
    markdown = []

    # Add the Hugo front matter
    markdown.append("---")
    markdown.append(f'title : "{get_title_from_file(abs_file_path)}"')
    markdown.append('summary: ""')
    markdown.append("draft: false")
    markdown.append("images: []")
    markdown.append("toc: true")
    markdown.append('layout: "single"')
    markdown.append("---")

    # Add the revision string
    revision_string = get_revision_string_from_file(abs_file_path)
    markdown.append(revision_string)
    markdown.append("\n\nNote that the code is given in full at the bottom of the page.\n\n")

    # Add the tutorial content
    markdown.append(test_markdown)

    # Append the full source code
    markdown.append(f"\n\n## Full code\n\n{format_code_lines(test_code)}")

    markdown_string = "\n".join(markdown)

    # Postprocess to remove empty lines at the end of code blocks
    markdown_string = markdown_string.replace("\n\n```\n", "\n```\n")

    # Postprocess to remove more than one empty line above a code block
    markdown_string = markdown_string.replace("\n\n\n```python\n", "\n\n```python\n")

    # Postprocess to remove any cases of 3 blank lines
    markdown_string = markdown_string.replace("\n\n\n\n", "\n\n")

    # Postprocess to remove any cases of 2 blank lines
    markdown_string = markdown_string.replace("\n\n\n", "\n\n")

    tutorial_name = get_output_file_name(abs_file_path)
    markdown_file_path = PYCHASTE_TUTORIAL_DIR / f"{tutorial_name}.md"
    with open(markdown_file_path, "w", encoding="utf-8") as file:
        file.write(markdown_string)

    # Process the Jupyter notebook file
    test_notebook["cells"].insert(0, new_markdown_cell(revision_string))

    notebook_file_path = PYCHASTE_TUTORIAL_DIR / f"{tutorial_name}.ipynb"
    with open(notebook_file_path, "w", encoding="utf-8") as file:
        nbformat.write(test_notebook, file)


if __name__ == "__main__":
    PYCHASTE_TUTORIAL_DIR.mkdir(parents=True, exist_ok=True)
    if any(PYCHASTE_TUTORIAL_DIR.iterdir()):
        warnings.warn(f"{PYCHASTE_TUTORIAL_DIR} is not empty")

    for tutorial_file in get_list_of_tutorial_files():
        write_tutorial(tutorial_file)
