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

// This file is auto-generated; manual changes will be overwritten.
// To make enduring changes, see pychaste/dynamic/config.yaml.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "PythonVtkObjectConverters.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "PythonUblasObjectConverters.hpp"
#include "VtkScene.hpp"

#include "VtkScene_3.cppwg.hpp"

namespace py = pybind11;
typedef VtkScene<3> VtkScene_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class VtkScene_3_Overrides : public VtkScene_3
{
public:
    using VtkScene_3::VtkScene;
    void ResetRenderer(unsigned int timeStep) override
    {
        PYBIND11_OVERRIDE(
            void,
            VtkScene_3,
            ResetRenderer,
            timeStep);
    }
};

void register_VtkScene_3_class(py::module &m)
{
    py::class_<VtkScene_3, VtkScene_3_Overrides, boost::shared_ptr<VtkScene_3>>(m, "VtkScene_3")
        .def(py::init<>())
        .def("End",
            (void(VtkScene_3::*)()) &VtkScene_3::End,
            " ")
        .def("GetSceneAsCharBuffer",
            (::vtkSmartPointer<vtkUnsignedCharArray>(VtkScene_3::*)()) &VtkScene_3::GetSceneAsCharBuffer,
            " ")
        .def("GetRenderer",
            (::vtkSmartPointer<vtkRenderer>(VtkScene_3::*)()) &VtkScene_3::GetRenderer,
            " ")
        .def("GetCellPopulationActorGenerator",
            (::boost::shared_ptr<CellPopulationPyChasteActorGenerator<3>>(VtkScene_3::*)()) &VtkScene_3::GetCellPopulationActorGenerator,
            " ")
        .def("ResetRenderer",
            (void(VtkScene_3::*)(unsigned int)) &VtkScene_3::ResetRenderer,
            " ", py::arg("timeStep") = 0)
        .def("Start",
            (void(VtkScene_3::*)()) &VtkScene_3::Start,
            " ")
        .def("SetCellPopulation",
            (void(VtkScene_3::*)(::boost::shared_ptr<AbstractCellPopulation<3, 3>>)) &VtkScene_3::SetCellPopulation,
            " ", py::arg("pCellPopulation"))
        .def("SetOutputFilePath",
            (void(VtkScene_3::*)(::std::string const &)) &VtkScene_3::SetOutputFilePath,
            " ", py::arg("rPath"))
        .def("SetIsInteractive",
            (void(VtkScene_3::*)(bool)) &VtkScene_3::SetIsInteractive,
            " ", py::arg("isInteractive"))
        .def("SetSaveAsAnimation",
            (void(VtkScene_3::*)(bool)) &VtkScene_3::SetSaveAsAnimation,
            " ", py::arg("saveAsAnimation"))
        .def("SetSaveAsImages",
            (void(VtkScene_3::*)(bool)) &VtkScene_3::SetSaveAsImages,
            " ", py::arg("saveAsImages"))
        .def("StartInteractiveEventHandler",
            (void(VtkScene_3::*)()) &VtkScene_3::StartInteractiveEventHandler,
            " ")
    ;
}
