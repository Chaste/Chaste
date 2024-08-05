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
#include "AbstractChasteRegion.hpp"

#include "AbstractChasteRegion2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractChasteRegion<2 > AbstractChasteRegion2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractChasteRegion2_Overrides : public AbstractChasteRegion2{
    public:
    using AbstractChasteRegion2::AbstractChasteRegion;
    void Destroy() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractChasteRegion2,
            Destroy,
            );
    }
    bool DoesContain(::ChastePoint<2> const & rPointToCheck) const  override {
        PYBIND11_OVERRIDE_PURE(
            bool,
            AbstractChasteRegion2,
            DoesContain,
                    rPointToCheck);
    }

};
void register_AbstractChasteRegion2_class(py::module &m){
py::class_<AbstractChasteRegion2 , AbstractChasteRegion2_Overrides , boost::shared_ptr<AbstractChasteRegion2 >   >(m, "AbstractChasteRegion2")
        .def(py::init< >())
        .def(
            "Destroy",
            (void(AbstractChasteRegion2::*)()) &AbstractChasteRegion2::Destroy,
            " "  )
        .def(
            "DoesContain",
            (bool(AbstractChasteRegion2::*)(::ChastePoint<2> const &) const ) &AbstractChasteRegion2::DoesContain,
            " " , py::arg("rPointToCheck") )
    ;
}
