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
#include "ImmersedBoundaryElement.hpp"

#include "ImmersedBoundaryElement2_3.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryElement<2,3 > ImmersedBoundaryElement2_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundaryElement2_3_Overrides : public ImmersedBoundaryElement2_3{
    public:
    using ImmersedBoundaryElement2_3::ImmersedBoundaryElement;
    bool IsElementOnBoundary() const  override {
        PYBIND11_OVERRIDE(
            bool,
            ImmersedBoundaryElement2_3,
            IsElementOnBoundary,
            );
    }

};
void register_ImmersedBoundaryElement2_3_class(py::module &m){
py::class_<ImmersedBoundaryElement2_3 , ImmersedBoundaryElement2_3_Overrides , boost::shared_ptr<ImmersedBoundaryElement2_3 >  , MutableElement<2, 3>  >(m, "ImmersedBoundaryElement2_3")
        .def(py::init<unsigned int, ::std::vector<Node<3> *> const & >(), py::arg("index"), py::arg("rNodes"))
        .def(
            "SetFluidSource",
            (void(ImmersedBoundaryElement2_3::*)(::std::shared_ptr<FluidSource<3>>)) &ImmersedBoundaryElement2_3::SetFluidSource,
            " " , py::arg("fluidSource") )
        .def(
            "GetFluidSource",
            (::std::shared_ptr<FluidSource<3>>(ImmersedBoundaryElement2_3::*)()) &ImmersedBoundaryElement2_3::GetFluidSource,
            " "  )
        .def(
            "rGetCornerNodes",
            (::std::vector<Node<3> *> &(ImmersedBoundaryElement2_3::*)()) &ImmersedBoundaryElement2_3::rGetCornerNodes,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetAverageNodeSpacing",
            (double(ImmersedBoundaryElement2_3::*)()) &ImmersedBoundaryElement2_3::GetAverageNodeSpacing,
            " "  )
        .def(
            "SetAverageNodeSpacing",
            (void(ImmersedBoundaryElement2_3::*)(double)) &ImmersedBoundaryElement2_3::SetAverageNodeSpacing,
            " " , py::arg("averageNodeSpacing") )
        .def(
            "IsElementOnBoundary",
            (bool(ImmersedBoundaryElement2_3::*)() const ) &ImmersedBoundaryElement2_3::IsElementOnBoundary,
            " "  )
        .def(
            "SetIsBoundaryElement",
            (void(ImmersedBoundaryElement2_3::*)(bool)) &ImmersedBoundaryElement2_3::SetIsBoundaryElement,
            " " , py::arg("isBoundaryElement") )
    ;
}
