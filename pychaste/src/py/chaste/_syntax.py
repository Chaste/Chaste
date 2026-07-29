"""Syntax Module"""

__copyright__ = """Copyright (c) 2005-2026, University of Oxford.
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
import warnings
from collections.abc import Iterable
from typing import Dict, Tuple, Type


class TemplateClassDict:
    """
    Allows using class syntax like Foo[2, 2](...) in place of Foo_2_2(...)

    Usage:
    >>> Foo = TemplateClassDict({ ("2", "2"): Foo_2_2, ("3", "3"): Foo_3_3 })
    """

    def __init__(self, template_dict: Dict[Tuple[str, ...], Type]) -> None:
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


def _normalize_key(key) -> Tuple[str, ...]:
    """Normalize a template-argument subscript key to a tuple of strings.

    A scalar key becomes a 1-tuple; each argument maps to its ``__name__`` (for a
    class) or ``str`` otherwise - matching the suffix cppwg appends when it names
    an instantiated class or method. A string is treated as a single scalar (not
    iterated character by character).
    """
    if isinstance(key, str) or not isinstance(key, Iterable):
        key = (key,)
    return tuple(arg.__name__ if inspect.isclass(arg) else str(arg) for arg in key)


class TemplateMethod:
    """Subscript syntax for a templated method (the method analogue of
    TemplateClassDict).

    Assign it as a class attribute so ``obj.<base>[Arg]()`` dispatches to the
    mangled binding ``obj.<base><Arg>()`` that cppwg generates for the templated
    C++ method - e.g. ``population.AddCellWriter[CellVolumesWriter]()`` calls
    ``population.AddCellWriterCellVolumesWriter()``. When the name is also a plain
    (non-templated) overload, pass it as ``fallback`` so ``obj.<base>(...)`` keeps
    working alongside the subscript form.

    Usage:
    >>> Foo.AddCellWriter = TemplateMethod("AddCellWriter", Foo.AddCellWriter)
    """

    def __init__(self, base_name: str, fallback=None):
        self._base_name = base_name
        self._fallback = fallback

    def __get__(self, obj, owner=None):
        target = obj if obj is not None else owner
        return _BoundTemplateMethod(target, self._base_name, self._fallback)


class _BoundTemplateMethod:
    def __init__(self, target, base_name: str, fallback):
        self._target = target
        self._base_name = base_name
        self._fallback = fallback

    def __getitem__(self, key):
        # Mangled binding is <base>_<arg1>_<arg2>..., matching cppwg's Foo_2 style.
        suffix = "_" + "_".join(_normalize_key(key))
        return getattr(self._target, self._base_name + suffix)

    def __call__(self, *args, **kwargs):
        if self._fallback is None:
            raise TypeError(
                f"{self._base_name} is templated; use {self._base_name}[Arg](...)"
            )
        return self._fallback(self._target, *args, **kwargs)


class DeprecatedClass:
    """
    Warns when a deprecated class is used and switches to the correct class.

    Usage:
    >>> Foo2_2 = DeprecatedClass("Foo2_2", Foo_2_2)
    """
    def __init__(self, old_name: str, new_class: Type):
        self.old_name = old_name
        self.new_class = new_class

        self.new_syntax = self.new_class.__name__
        if "_" in self.new_syntax:
            # Recommend using Foo["2", "2"]() instead of Foo2_2()
            base_name, *params = self.new_syntax.split("_")
            params = [f'"{param}"' for param in params]
            self.new_syntax = f'{base_name}[{", ".join(params)}]'

    def __call__(self, *args, **kwargs):
        warnings.warn(
            f"{self.old_name} is deprecated and will be removed in a future version. "
            f"Please use {self.new_syntax} instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.new_class(*args, **kwargs)
