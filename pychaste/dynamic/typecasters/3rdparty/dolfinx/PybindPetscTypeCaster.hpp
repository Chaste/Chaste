// Copyright (C) 2017-2023 Chris Richardson and Garth N. Wells
//
// This file is part of DOLFINx (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    LGPL-3.0-or-later

// Modifications:
// This file, originally named "caster_petsc.h", has been modified
// to use pybind11 instead of nanobind.

#ifndef PYBINDPETSCTYPECASTER_HPP_
#define PYBINDPETSCTYPECASTER_HPP_

#include <pybind11/cast.h>
#include <pybind11/pybind11.h>

#include <petsc4py/petsc4py.h>
#include <petscmat.h>
#include <petscvec.h>

// pybind11 casters for PETSc/petsc4py objects

// Import petsc4py on demand
#define VERIFY_PETSC4PY_FROMPY(func) \
  if (!func)                         \
  {                                  \
    if (import_petsc4py() != 0)      \
      return false;                  \
  }

#define VERIFY_PETSC4PY_FROMCPP(func) \
  if (!func)                          \
  {                                   \
    if (import_petsc4py() != 0)       \
      return {};                      \
  }

// Macro for casting between PETSc and petsc4py objects
#define PETSC_CASTER_MACRO(TYPE, P4PYTYPE, NAME)                         \
  template <>                                                            \
  class type_caster<_p_##TYPE>                                           \
  {                                                                      \
  public:                                                                \
    PYBIND11_TYPE_CASTER(TYPE, const_name(#NAME));                       \
    bool load(handle src, bool)                                          \
    {                                                                    \
      VERIFY_PETSC4PY_FROMPY(PyPetsc##P4PYTYPE##_Get);                   \
      if (PyObject_TypeCheck(src.ptr(), &PyPetsc##P4PYTYPE##_Type) != 0) \
      {                                                                  \
        value = PyPetsc##P4PYTYPE##_Get(src.ptr());                      \
        return true;                                                     \
      }                                                                  \
      else                                                               \
        return false;                                                    \
    }                                                                    \
                                                                         \
    static handle cast(TYPE src, return_value_policy policy, handle)     \
    {                                                                    \
      VERIFY_PETSC4PY_FROMCPP(PyPetsc##P4PYTYPE##_New);                  \
      if (policy == return_value_policy::take_ownership)                 \
      {                                                                  \
        PyObject *obj = PyPetsc##P4PYTYPE##_New(src);                    \
        PetscObjectDereference((PetscObject)src);                        \
        return pybind11::handle(obj);                                    \
      }                                                                  \
      else if (policy == return_value_policy::automatic_reference or     \
               policy == return_value_policy::reference)                 \
      {                                                                  \
        PyObject *obj = PyPetsc##P4PYTYPE##_New(src);                    \
        return pybind11::handle(obj);                                    \
      }                                                                  \
      else                                                               \
      {                                                                  \
        return {};                                                       \
      }                                                                  \
    }                                                                    \
                                                                         \
    operator TYPE() { return value; }                                    \
  }

namespace pybind11::detail
{
  PETSC_CASTER_MACRO(Mat, Mat, mat);
  PETSC_CASTER_MACRO(Vec, Vec, vec);
} // namespace pybind11::detail

#endif // PYBINDPETSCTYPECASTER_HPP_
