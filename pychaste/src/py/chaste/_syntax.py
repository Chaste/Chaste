"""Syntax module.

C++ templates only exist once instantiated at concrete arguments, so cppwg wraps
each instantiation as its own class with a mangled name: Node<2> → Node_2,
Node<3> → Node_3. That's correct but we'd rather write Node[2] in Python,
mirroring C++'s Node<2>. This is a helper module to add that subscript syntax,
holding two helpers — TemplateClass (a stub base class using __class_getitem__,
like list[int]) and TemplateMethod (a descriptor) — plus the shared key
normalization that resolves a subscript to the concrete instantiation.
"""

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

from collections.abc import Iterable


def _normalize_key(key) -> tuple[str, ...]:
    """Normalize a template-argument subscript key to a tuple of strings.

    A scalar key becomes a 1-tuple; each argument maps to its ``__name__`` (for a
    class, including a TemplateClass stub) or ``str`` otherwise - matching the
    suffix cppwg appends when it names an instantiated class or method. A string
    is treated as a single scalar (not iterated character by character).
    """
    if isinstance(key, str) or not isinstance(key, Iterable):
        key = (key,)
    return tuple(arg.__name__ if hasattr(arg, "__name__") else str(arg) for arg in key)


class TemplateClass:
    """Base for a stub class that gives a templated class subscript syntax.

    Subclass it with an ``_instantiations`` map from template-argument tuples to
    the concrete wrapped classes; ``Foo[args]`` then resolves the instantiation -
    e.g. ``Node[2]`` -> ``Node_2`` - mirroring how ``list[int]`` works via
    ``__class_getitem__``. Subclassing (rather than an instance) makes ``Foo`` a
    real class object with a ``__name__``, so it can itself be used as a template
    argument, e.g. ``population.AddCellWriter[CellVolumesWriter]()``. Keys are
    normalized once at subclass creation.

    Usage:
    >>> class Foo(TemplateClass):
    ...     _instantiations = {("2", "2"): Foo_2_2, ("3", "3"): Foo_3_3}
    """

    _instantiations = {}

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        cls._instantiations = {
            _normalize_key(args): cls_ for args, cls_ in cls._instantiations.items()
        }

    def __class_getitem__(cls, key):
        return cls._instantiations[_normalize_key(key)]


class TemplateMethod:
    """Subscript syntax for a templated method.

    TemplateMethod is a descriptor: set it as a class attribute, then use the
    ``obj.<base>[T]()`` subscript form to reach the per-instantiation binding
    ``obj.<base>_<T>()`` that cppwg generates. When the name is also a plain
    (non-templated) overload, pass it as ``fallback`` so ``obj.<base>(...)`` keeps
    working alongside the subscript form.

    Usage:
    >>> Foo.Bar = TemplateMethod("Bar")
    >>> foo_obj.Bar[T]()

    If ``Bar`` also has a plain overload, keep it as the fallback:
    >>> Foo.Bar = TemplateMethod("Bar", Foo.Bar)
    >>> foo_obj.Bar(arg)  # the plain overload, via the fallback
    """

    def __init__(self, base_name: str, fallback=None):
        self._base_name = base_name  # e.g. "Bar" for foo_obj.Bar[T]()
        self._fallback = fallback

    def __get__(self, obj, owner=None):
        # Bar is a descriptor on the class, so accessing ``foo_obj.Bar`` triggers
        # __get__, returning a _BoundTemplateMethod bound to the instance (obj).
        # The class (owner) is used when obj is None, e.g. when accessed on the
        # class itself ``Foo.Bar[T]``.
        target = obj if obj is not None else owner
        return _BoundTemplateMethod(target, self._base_name, self._fallback)


class _BoundTemplateMethod:
    def __init__(self, target, base_name: str, fallback):
        self._target = target  # normally an instance, but can be a class
        self._base_name = base_name  # e.g. "Bar" for foo_obj.Bar[T]()
        self._fallback = fallback

    def __getitem__(self, key):
        # The [T] subscript on ``foo_obj.Bar[T]()`` triggers __getitem__,
        # returning the target.Bar_T method, the binding generated by cppwg.
        suffix = "_" + "_".join(_normalize_key(key))  # e.g. _CellVolumesWriter
        return getattr(self._target, self._base_name + suffix)

    def __call__(self, *args, **kwargs):
        # ``foo_obj.Bar(...)`` with no subscript calls the plain overload kept as
        # the fallback; with no fallback the name is purely templated, so point
        # the caller at the subscript form.
        if self._fallback is None:
            raise TypeError(
                f"{self._base_name} is templated; use {self._base_name}[Arg](...)"
            )
        return self._fallback(self._target, *args, **kwargs)
