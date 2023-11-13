"""Copyright (c) 2005-2023, University of Oxford.
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

from collections import namedtuple
from pathlib import Path


CHASTE_SRC_DIR = Path(__file__).resolve().parent.parent.parent
HOW_TO_TAG = 'HOW_TO_TAG'

HowTo = namedtuple('HowTo', ['primary', 'secondary', 'line_number', 'test_name', 'github_link', 'how_to_text'])


def get_list_of_test_files():
    return [file.resolve() for file in CHASTE_SRC_DIR.rglob('Test*.hpp')]


def stop_condition(line_to_check: str):
    """

    :return:
    """
    line_to_check = line_to_check.strip()
    if line_to_check == '' or line_to_check == '*' or line_to_check == '*/':
        return True
    return False


def process_line(line_to_process: str) -> str:
    """

    :return:
    """
    return ' ' + line_to_process.replace('*', '').replace('/', '').strip()


def parse_how_to_line(line_to_parse: str) -> tuple[str, str]:
    """

    :return:
    """
    _, _, after = line_to_parse.partition(HOW_TO_TAG)
    _component, _, _category = after.strip().partition('/')
    return _component, _category


def form_github_link(file_path: Path, line_num: int) -> str:
    """

    :return:
    """
    return f'https://github.com/Chaste/Chaste/blob/develop/{file_path.relative_to(CHASTE_SRC_DIR)}#{line_num}'


list_of_how_tos = []

for test_file in get_list_of_test_files():

    with open(test_file, 'r') as f:
        for line_number, line in enumerate(f):

            if HOW_TO_TAG in line:

                primary, secondary = parse_how_to_line(line)

                how_to_text = ''
                finished = False
                while not finished:
                    next_line = next(f).strip()
                    how_to_text += process_line(next_line)
                    finished = stop_condition(next_line)

                how_to = HowTo(primary=primary, secondary=secondary, line_number=line_number, test_name=test_file.name,
                               github_link=form_github_link(test_file, line_number), how_to_text=how_to_text.strip())

                list_of_how_tos.append(how_to)

print(list_of_how_tos)
