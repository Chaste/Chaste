/*

Copyright (c) 2005-2024, University of Oxford.
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

*/

#ifndef PYTHONVTKOBJECTCONVERTERS_HPP_
#define PYTHONVTKOBJECTCONVERTERS_HPP_

#include <vtkObjectBase.h>
#include <vtkOpenGLRenderer.h>
#include <vtkPythonUtil.h>
#include <vtkSmartPointer.h>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <typeinfo>

/**
 *  VTK Conversion, from SMTK Source with copyright
//=========================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//=========================================================================
*/

#define PYBIND11_VTK_TYPECASTER(VTK_OBJ)                                                                \
  namespace pybind11                                                                                    \
  {                                                                                                     \
    namespace detail                                                                                    \
    {                                                                                                   \
      template <>                                                                                       \
      struct type_caster<vtkSmartPointer<VTK_OBJ>>                                                      \
      {                                                                                                 \
      protected:                                                                                        \
        VTK_OBJ *value;                                                                                 \
                                                                                                        \
      public:                                                                                           \
        static constexpr auto name = _(#VTK_OBJ);                                                       \
        operator VTK_OBJ *() { return value; }                                                          \
        operator VTK_OBJ &() { return *value; }                                                         \
        template <typename _T>                                                                          \
        using cast_op_type = pybind11::detail::cast_op_type<_T>;                                        \
        bool load(handle src, bool)                                                                     \
        {                                                                                               \
          value = dynamic_cast<VTK_OBJ *>(vtkPythonUtil::GetPointerFromObject(src.ptr(), #VTK_OBJ));    \
          if (!value)                                                                                   \
          {                                                                                             \
            PyErr_Clear();                                                                              \
            throw reference_cast_error();                                                               \
          }                                                                                             \
          return value != nullptr;                                                                      \
        }                                                                                               \
        static handle cast(const vtkSmartPointer<VTK_OBJ> &src,                                         \
                           return_value_policy,                                                         \
                           handle)                                                                      \
        {                                                                                               \
          PyObject *obj = vtkPythonUtil::GetObjectFromPointer(const_cast<VTK_OBJ *>(src.GetPointer())); \
          return obj;                                                                                   \
        }                                                                                               \
      };                                                                                                \
    }                                                                                                   \
  }

PYBIND11_DECLARE_HOLDER_TYPE(T, vtkSmartPointer<T>);
// Only needed if the type's `.get()` goes by another name
namespace PYBIND11_NAMESPACE
{
  namespace detail
  {
    template <typename T>
    struct holder_helper<vtkSmartPointer<T>> // <-- specialization
    {
      static const T *get(const vtkSmartPointer<T> &p) { return p.Get(); }
    };
  }
}

PYBIND11_VTK_TYPECASTER(vtkRenderer);
PYBIND11_VTK_TYPECASTER(vtkOpenGLRenderer);
PYBIND11_VTK_TYPECASTER(vtkUnsignedCharArray);

#undef PYBIND11_VTK_TYPECASTER

#endif /*PYTHONVTKOBJECTCONVERTERS_HPP_*/
