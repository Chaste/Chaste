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
#include "RandomNumberGenerator.hpp"

#include "RandomNumberGenerator.cppwg.hpp"

namespace py = pybind11;
typedef RandomNumberGenerator RandomNumberGenerator;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

void register_RandomNumberGenerator_class(py::module &m){
py::class_<RandomNumberGenerator  , boost::shared_ptr<RandomNumberGenerator >   >(m, "RandomNumberGenerator")
        .def(
            "StandardNormalRandomDeviate",
            (double(RandomNumberGenerator::*)()) &RandomNumberGenerator::StandardNormalRandomDeviate,
            " "  )
        .def(
            "NormalRandomDeviate",
            (double(RandomNumberGenerator::*)(double, double)) &RandomNumberGenerator::NormalRandomDeviate,
            " " , py::arg("mean"), py::arg("stdDev") )
        .def(
            "ranf",
            (double(RandomNumberGenerator::*)()) &RandomNumberGenerator::ranf,
            " "  )
        .def(
            "GammaRandomDeviate",
            (double(RandomNumberGenerator::*)(double, double)) &RandomNumberGenerator::GammaRandomDeviate,
            " " , py::arg("shape"), py::arg("scale") )
        .def(
            "ExponentialRandomDeviate",
            (double(RandomNumberGenerator::*)(double)) &RandomNumberGenerator::ExponentialRandomDeviate,
            " " , py::arg("scale") )
        .def(
            "randMod",
            (unsigned int(RandomNumberGenerator::*)(unsigned int)) &RandomNumberGenerator::randMod,
            " " , py::arg("base") )
        .def(
            "Shuffle",
            (void(RandomNumberGenerator::*)(unsigned int, ::std::vector<unsigned int> &)) &RandomNumberGenerator::Shuffle,
            " " , py::arg("num"), py::arg("rValues") )
        .def_static(
            "Instance",
            (::RandomNumberGenerator *(*)()) &RandomNumberGenerator::Instance,
            " "  , py::return_value_policy::reference)
        .def_static(
            "Destroy",
            (void(*)()) &RandomNumberGenerator::Destroy,
            " "  )
        .def(
            "Reseed",
            (void(RandomNumberGenerator::*)(unsigned int)) &RandomNumberGenerator::Reseed,
            " " , py::arg("seed") )
    ;
}
