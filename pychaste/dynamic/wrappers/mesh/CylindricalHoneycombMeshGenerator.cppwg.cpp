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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"

#include "CylindricalHoneycombMeshGenerator.cppwg.hpp"

namespace py = pybind11;
typedef CylindricalHoneycombMeshGenerator CylindricalHoneycombMeshGenerator;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::shared_ptr<MutableMesh<2, 2>> _boost_shared_ptr_lt_MutableMesh_lt_2_2_gt__gt_;

class CylindricalHoneycombMeshGenerator_Overrides : public CylindricalHoneycombMeshGenerator{
    public:
    using CylindricalHoneycombMeshGenerator::CylindricalHoneycombMeshGenerator;
    ::boost::shared_ptr<MutableMesh<2, 2>> GetMesh() override {
        PYBIND11_OVERRIDE(
            _boost_shared_ptr_lt_MutableMesh_lt_2_2_gt__gt_,
            CylindricalHoneycombMeshGenerator,
            GetMesh,
            );
    }

};
void register_CylindricalHoneycombMeshGenerator_class(py::module &m){
py::class_<CylindricalHoneycombMeshGenerator , CylindricalHoneycombMeshGenerator_Overrides , boost::shared_ptr<CylindricalHoneycombMeshGenerator >  , HoneycombMeshGenerator  >(m, "CylindricalHoneycombMeshGenerator")
        .def(py::init<unsigned int, unsigned int, unsigned int, double >(), py::arg("numNodesAlongWidth"), py::arg("numNodesAlongLength"), py::arg("ghosts") = 3, py::arg("scaleFactor") = 1.)
        .def(
            "GetMesh",
            (::boost::shared_ptr<MutableMesh<2, 2>>(CylindricalHoneycombMeshGenerator::*)()) &CylindricalHoneycombMeshGenerator::GetMesh,
            " "  )
        .def(
            "GetCylindricalMesh",
            (::boost::shared_ptr<Cylindrical2dMesh>(CylindricalHoneycombMeshGenerator::*)()) &CylindricalHoneycombMeshGenerator::GetCylindricalMesh,
            " "  )
    ;
}
