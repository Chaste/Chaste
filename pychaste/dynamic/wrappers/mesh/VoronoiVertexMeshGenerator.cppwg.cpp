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
#include "VoronoiVertexMeshGenerator.hpp"

#include "VoronoiVertexMeshGenerator.cppwg.hpp"

namespace py = pybind11;
typedef VoronoiVertexMeshGenerator VoronoiVertexMeshGenerator;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::shared_ptr<MutableVertexMesh<2, 2>> _boost_shared_ptr_lt_MutableVertexMesh_lt_2_2_gt__gt_;
typedef ::boost::shared_ptr<MutableVertexMesh<2, 2>> _boost_shared_ptr_lt_MutableVertexMesh_lt_2_2_gt__gt_;
typedef ::boost::shared_ptr<Toroidal2dVertexMesh> _boost_shared_ptr_lt_Toroidal2dVertexMesh_gt_;

class VoronoiVertexMeshGenerator_Overrides : public VoronoiVertexMeshGenerator{
    public:
    using VoronoiVertexMeshGenerator::VoronoiVertexMeshGenerator;
    ::boost::shared_ptr<MutableVertexMesh<2, 2>> GetMesh() override {
        PYBIND11_OVERRIDE(
            _boost_shared_ptr_lt_MutableVertexMesh_lt_2_2_gt__gt_,
            VoronoiVertexMeshGenerator,
            GetMesh,
            );
    }
    ::boost::shared_ptr<MutableVertexMesh<2, 2>> GetMeshAfterReMesh() override {
        PYBIND11_OVERRIDE(
            _boost_shared_ptr_lt_MutableVertexMesh_lt_2_2_gt__gt_,
            VoronoiVertexMeshGenerator,
            GetMeshAfterReMesh,
            );
    }
    ::boost::shared_ptr<Toroidal2dVertexMesh> GetToroidalMesh() override {
        PYBIND11_OVERRIDE(
            _boost_shared_ptr_lt_Toroidal2dVertexMesh_gt_,
            VoronoiVertexMeshGenerator,
            GetToroidalMesh,
            );
    }

};
void register_VoronoiVertexMeshGenerator_class(py::module &m){
py::class_<VoronoiVertexMeshGenerator , VoronoiVertexMeshGenerator_Overrides , boost::shared_ptr<VoronoiVertexMeshGenerator >   >(m, "VoronoiVertexMeshGenerator")
        .def(py::init<unsigned int, unsigned int, unsigned int, double >(), py::arg("numElementsX"), py::arg("numElementsY"), py::arg("numRelaxationSteps"), py::arg("elementTargetArea") = 1.)
        .def(py::init< >())
        .def(
            "GenerateVoronoiMesh",
            (void(VoronoiVertexMeshGenerator::*)()) &VoronoiVertexMeshGenerator::GenerateVoronoiMesh,
            " "  )
        .def(
            "GetMesh",
            (::boost::shared_ptr<MutableVertexMesh<2, 2>>(VoronoiVertexMeshGenerator::*)()) &VoronoiVertexMeshGenerator::GetMesh,
            " "  )
        .def(
            "GetMeshAfterReMesh",
            (::boost::shared_ptr<MutableVertexMesh<2, 2>>(VoronoiVertexMeshGenerator::*)()) &VoronoiVertexMeshGenerator::GetMeshAfterReMesh,
            " "  )
        .def(
            "GetToroidalMesh",
            (::boost::shared_ptr<Toroidal2dVertexMesh>(VoronoiVertexMeshGenerator::*)()) &VoronoiVertexMeshGenerator::GetToroidalMesh,
            " "  )
        .def(
            "GetPolygonDistribution",
            (::std::vector<double>(VoronoiVertexMeshGenerator::*)()) &VoronoiVertexMeshGenerator::GetPolygonDistribution,
            " "  )
        .def(
            "GetAreaCoefficientOfVariation",
            (double(VoronoiVertexMeshGenerator::*)()) &VoronoiVertexMeshGenerator::GetAreaCoefficientOfVariation,
            " "  )
        .def(
            "RefreshSeedsAndRegenerateMesh",
            (void(VoronoiVertexMeshGenerator::*)()) &VoronoiVertexMeshGenerator::RefreshSeedsAndRegenerateMesh,
            " "  )
        .def(
            "SetMaxExpectedNumSidesPerPolygon",
            (void(VoronoiVertexMeshGenerator::*)(unsigned int)) &VoronoiVertexMeshGenerator::SetMaxExpectedNumSidesPerPolygon,
            " " , py::arg("maxExpectedNumSidesPerPolygon") )
        .def(
            "GetMaxExpectedNumSidesPerPolygon",
            (unsigned int(VoronoiVertexMeshGenerator::*)()) &VoronoiVertexMeshGenerator::GetMaxExpectedNumSidesPerPolygon,
            " "  )
    ;
}
