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
#include "PythonUblasObjectConverters.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "SphereGeometryBoundaryCondition.hpp"

#include "SphereGeometryBoundaryCondition3.cppwg.hpp"

namespace py = pybind11;
typedef SphereGeometryBoundaryCondition<3 > SphereGeometryBoundaryCondition3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class SphereGeometryBoundaryCondition3_Overrides : public SphereGeometryBoundaryCondition3{
    public:
    using SphereGeometryBoundaryCondition3::SphereGeometryBoundaryCondition;
    void ImposeBoundaryCondition(::std::map<Node<3> *, boost::numeric::ublas::c_vector<double, 3>> const & rOldLocations) override {
        PYBIND11_OVERRIDE(
            void,
            SphereGeometryBoundaryCondition3,
            ImposeBoundaryCondition,
                    rOldLocations);
    }
    bool VerifyBoundaryCondition() override {
        PYBIND11_OVERRIDE(
            bool,
            SphereGeometryBoundaryCondition3,
            VerifyBoundaryCondition,
            );
    }
    void OutputCellPopulationBoundaryConditionParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            SphereGeometryBoundaryCondition3,
            OutputCellPopulationBoundaryConditionParameters,
                    rParamsFile);
    }

};
void register_SphereGeometryBoundaryCondition3_class(py::module &m){
py::class_<SphereGeometryBoundaryCondition3 , SphereGeometryBoundaryCondition3_Overrides , boost::shared_ptr<SphereGeometryBoundaryCondition3 >  , AbstractCellPopulationBoundaryCondition<3>  >(m, "SphereGeometryBoundaryCondition3")
        .def(py::init<::AbstractCellPopulation<3> *, ::boost::numeric::ublas::c_vector<double, 3>, double, double >(), py::arg("pCellPopulation"), py::arg("centre"), py::arg("radius"), py::arg("distance") = 1.0000000000000001E-5)
        .def(
            "rGetCentreOfSphere",
            (::boost::numeric::ublas::c_vector<double, 3> const &(SphereGeometryBoundaryCondition3::*)() const ) &SphereGeometryBoundaryCondition3::rGetCentreOfSphere,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetRadiusOfSphere",
            (double(SphereGeometryBoundaryCondition3::*)() const ) &SphereGeometryBoundaryCondition3::GetRadiusOfSphere,
            " "  )
        .def(
            "ImposeBoundaryCondition",
            (void(SphereGeometryBoundaryCondition3::*)(::std::map<Node<3> *, boost::numeric::ublas::c_vector<double, 3>> const &)) &SphereGeometryBoundaryCondition3::ImposeBoundaryCondition,
            " " , py::arg("rOldLocations") )
        .def(
            "VerifyBoundaryCondition",
            (bool(SphereGeometryBoundaryCondition3::*)()) &SphereGeometryBoundaryCondition3::VerifyBoundaryCondition,
            " "  )
        .def(
            "OutputCellPopulationBoundaryConditionParameters",
            (void(SphereGeometryBoundaryCondition3::*)(::out_stream &)) &SphereGeometryBoundaryCondition3::OutputCellPopulationBoundaryConditionParameters,
            " " , py::arg("rParamsFile") )
    ;
}
