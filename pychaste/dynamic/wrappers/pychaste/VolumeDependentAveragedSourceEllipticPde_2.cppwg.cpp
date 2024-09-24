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

// This file is auto-generated. Manual changes will be overwritten. For changes
// to persist, update the configuration in pychaste/dynamic/config.yaml.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "VolumeDependentAveragedSourceEllipticPde.hpp"

#include "VolumeDependentAveragedSourceEllipticPde_2.cppwg.hpp"

namespace py = pybind11;
typedef VolumeDependentAveragedSourceEllipticPde<2> VolumeDependentAveragedSourceEllipticPde_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class VolumeDependentAveragedSourceEllipticPde_2_Overrides : public VolumeDependentAveragedSourceEllipticPde_2
{
public:
    using VolumeDependentAveragedSourceEllipticPde_2::VolumeDependentAveragedSourceEllipticPde;
    void SetupSourceTerms(::TetrahedralMesh<2, 2> & rCoarseMesh, ::std::map<boost::shared_ptr<Cell>, unsigned int> * pCellPdeElementMap) override
    {
        PYBIND11_OVERRIDE(
            void,
            VolumeDependentAveragedSourceEllipticPde_2,
            SetupSourceTerms,
            rCoarseMesh,
            pCellPdeElementMap);
    }
};

void register_VolumeDependentAveragedSourceEllipticPde_2_class(py::module &m)
{
    py::class_<VolumeDependentAveragedSourceEllipticPde_2, VolumeDependentAveragedSourceEllipticPde_2_Overrides, boost::shared_ptr<VolumeDependentAveragedSourceEllipticPde_2>, AveragedSourceEllipticPde<2>>(m, "VolumeDependentAveragedSourceEllipticPde_2")
        .def(py::init<::AbstractCellPopulation<2> &, double>(), py::arg("rCellPopulation"), py::arg("coefficient") = 0.)
        .def("SetupSourceTerms",
            (void(VolumeDependentAveragedSourceEllipticPde_2::*)(::TetrahedralMesh<2, 2> &, ::std::map<boost::shared_ptr<Cell>, unsigned int> *)) &VolumeDependentAveragedSourceEllipticPde_2::SetupSourceTerms,
            " ", py::arg("rCoarseMesh"), py::arg("pCellPdeElementMap") = nullptr)
    ;
}
