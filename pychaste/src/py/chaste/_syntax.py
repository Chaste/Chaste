"""Syntax Module"""

__copyright__ = """Copyright (c) 2005-2024, University of Oxford.
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

import inspect
from collections.abc import Iterable
from typing import Dict, Tuple, Type


class TemplateClassDict:
    """
    Allows using class syntax like Foo[2, 2](...) in place of Foo_2_2(...)

    Usage:
    >>> Foo = TemplateClassDict({ ("2", "2"): Foo_2_2, ("3", "3"): Foo_3_3 })
    """

    def __init__(self, template_dict: Dict[Tuple[str, ...], Type]):
        """
        :param template_dict: A dictionary mapping template arg tuples to classes
        """
        self._dict = {}
        for arg_tuple, cls in template_dict.items():
            if not inspect.isclass(cls):
                raise TypeError("Expected class, got {}".format(type(cls)))
            if not isinstance(arg_tuple, Iterable):
                arg_tuple = (arg_tuple,)
            key = tuple(
                arg.__name__ if inspect.isclass(arg) else str(arg) for arg in arg_tuple
            )
            self._dict[key] = cls

    def __getitem__(self, arg_tuple: Tuple[str, ...]) -> Type:
        if not isinstance(arg_tuple, Iterable):
            arg_tuple = (arg_tuple,)
        key = tuple(
            arg.__name__ if inspect.isclass(arg) else str(arg) for arg in arg_tuple
        )
        return self._dict[key]
