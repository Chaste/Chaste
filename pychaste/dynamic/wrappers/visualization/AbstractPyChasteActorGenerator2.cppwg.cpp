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
#include "AbstractPyChasteActorGenerator.hpp"

#include "AbstractPyChasteActorGenerator2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractPyChasteActorGenerator<2 > AbstractPyChasteActorGenerator2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractPyChasteActorGenerator2_Overrides : public AbstractPyChasteActorGenerator2{
    public:
    using AbstractPyChasteActorGenerator2::AbstractPyChasteActorGenerator;
    void AddActor(::vtkSmartPointer<vtkRenderer> pRenderer) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractPyChasteActorGenerator2,
            AddActor,
                    pRenderer);
    }

};
void register_AbstractPyChasteActorGenerator2_class(py::module &m){
py::class_<AbstractPyChasteActorGenerator2 , AbstractPyChasteActorGenerator2_Overrides , boost::shared_ptr<AbstractPyChasteActorGenerator2 >   >(m, "AbstractPyChasteActorGenerator2")
        .def(py::init< >())
        .def(
            "GetColorTransferFunction",
            (::vtkSmartPointer<vtkColorTransferFunction>(AbstractPyChasteActorGenerator2::*)()) &AbstractPyChasteActorGenerator2::GetColorTransferFunction,
            " "  )
        .def(
            "GetDiscreteColorTransferFunction",
            (::vtkSmartPointer<vtkColorTransferFunction>(AbstractPyChasteActorGenerator2::*)()) &AbstractPyChasteActorGenerator2::GetDiscreteColorTransferFunction,
            " "  )
        .def(
            "GetScaleBar",
            (::vtkSmartPointer<vtkScalarBarActor>(AbstractPyChasteActorGenerator2::*)()) &AbstractPyChasteActorGenerator2::GetScaleBar,
            " "  )
        .def(
            "AddActor",
            (void(AbstractPyChasteActorGenerator2::*)(::vtkSmartPointer<vtkRenderer>)) &AbstractPyChasteActorGenerator2::AddActor,
            " " , py::arg("pRenderer") )
        .def(
            "SetShowEdges",
            (void(AbstractPyChasteActorGenerator2::*)(bool)) &AbstractPyChasteActorGenerator2::SetShowEdges,
            " " , py::arg("show") )
        .def(
            "SetShowPoints",
            (void(AbstractPyChasteActorGenerator2::*)(bool)) &AbstractPyChasteActorGenerator2::SetShowPoints,
            " " , py::arg("show") )
        .def(
            "SetShowVolume",
            (void(AbstractPyChasteActorGenerator2::*)(bool)) &AbstractPyChasteActorGenerator2::SetShowVolume,
            " " , py::arg("show") )
        .def(
            "SetEdgeColor",
            (void(AbstractPyChasteActorGenerator2::*)(::boost::numeric::ublas::c_vector<double, 3> const &)) &AbstractPyChasteActorGenerator2::SetEdgeColor,
            " " , py::arg("rColor") )
        .def(
            "SetPointColor",
            (void(AbstractPyChasteActorGenerator2::*)(::boost::numeric::ublas::c_vector<double, 3> const &)) &AbstractPyChasteActorGenerator2::SetPointColor,
            " " , py::arg("rColor") )
        .def(
            "SetVolumeColor",
            (void(AbstractPyChasteActorGenerator2::*)(::boost::numeric::ublas::c_vector<double, 3> const &)) &AbstractPyChasteActorGenerator2::SetVolumeColor,
            " " , py::arg("rColor") )
        .def(
            "SetVolumeOpacity",
            (void(AbstractPyChasteActorGenerator2::*)(double)) &AbstractPyChasteActorGenerator2::SetVolumeOpacity,
            " " , py::arg("opacity") )
        .def(
            "SetPointSize",
            (void(AbstractPyChasteActorGenerator2::*)(double)) &AbstractPyChasteActorGenerator2::SetPointSize,
            " " , py::arg("size") )
        .def(
            "SetEdgeSize",
            (void(AbstractPyChasteActorGenerator2::*)(double)) &AbstractPyChasteActorGenerator2::SetEdgeSize,
            " " , py::arg("size") )
        .def(
            "SetDataLabel",
            (void(AbstractPyChasteActorGenerator2::*)(::std::string const &)) &AbstractPyChasteActorGenerator2::SetDataLabel,
            " " , py::arg("rLabel") )
        .def(
            "SetShowScaleBar",
            (void(AbstractPyChasteActorGenerator2::*)(double)) &AbstractPyChasteActorGenerator2::SetShowScaleBar,
            " " , py::arg("show") )
    ;
}
