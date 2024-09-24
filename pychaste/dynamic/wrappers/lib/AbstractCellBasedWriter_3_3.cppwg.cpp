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
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractCellBasedWriter.hpp"

#include "AbstractCellBasedWriter_3_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCellBasedWriter<3, 3> AbstractCellBasedWriter_3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractCellBasedWriter_3_3_Overrides : public AbstractCellBasedWriter_3_3
{
public:
    using AbstractCellBasedWriter_3_3::AbstractCellBasedWriter;
    void OpenOutputFile(::OutputFileHandler & rOutputFileHandler) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellBasedWriter_3_3,
            OpenOutputFile,
            rOutputFileHandler);
    }
    void WriteTimeStamp() override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellBasedWriter_3_3,
            WriteTimeStamp,
            );
    }
    void WriteNewline() override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellBasedWriter_3_3,
            WriteNewline,
            );
    }
};

void register_AbstractCellBasedWriter_3_3_class(py::module &m)
{
    py::class_<AbstractCellBasedWriter_3_3, AbstractCellBasedWriter_3_3_Overrides, boost::shared_ptr<AbstractCellBasedWriter_3_3>, Identifiable>(m, "AbstractCellBasedWriter_3_3")
        .def(py::init<::std::string const &>(), py::arg("rFileName"))
        .def("CloseFile",
            (void(AbstractCellBasedWriter_3_3::*)()) &AbstractCellBasedWriter_3_3::CloseFile,
            " ")
        .def("OpenOutputFile",
            (void(AbstractCellBasedWriter_3_3::*)(::OutputFileHandler &)) &AbstractCellBasedWriter_3_3::OpenOutputFile,
            " ", py::arg("rOutputFileHandler"))
        .def("OpenOutputFileForAppend",
            (void(AbstractCellBasedWriter_3_3::*)(::OutputFileHandler &)) &AbstractCellBasedWriter_3_3::OpenOutputFileForAppend,
            " ", py::arg("rOutputFileHandler"))
        .def("WriteTimeStamp",
            (void(AbstractCellBasedWriter_3_3::*)()) &AbstractCellBasedWriter_3_3::WriteTimeStamp,
            " ")
        .def("WriteNewline",
            (void(AbstractCellBasedWriter_3_3::*)()) &AbstractCellBasedWriter_3_3::WriteNewline,
            " ")
        .def("SetFileName",
            (void(AbstractCellBasedWriter_3_3::*)(::std::string)) &AbstractCellBasedWriter_3_3::SetFileName,
            " ", py::arg("fileName"))
        .def("GetFileName",
            (::std::string(AbstractCellBasedWriter_3_3::*)()) &AbstractCellBasedWriter_3_3::GetFileName,
            " ")
    ;
}
