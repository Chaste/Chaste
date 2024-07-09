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
