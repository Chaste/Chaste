#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "Goldbeter1991OdeSystem.hpp"

#include "Goldbeter1991OdeSystem.cppwg.hpp"

namespace py = pybind11;
typedef Goldbeter1991OdeSystem Goldbeter1991OdeSystem;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class Goldbeter1991OdeSystem_Overrides : public Goldbeter1991OdeSystem{
    public:
    using Goldbeter1991OdeSystem::Goldbeter1991OdeSystem;
    void EvaluateYDerivatives(double time, ::std::vector<double> const & rY, ::std::vector<double> & rDY) override {
        PYBIND11_OVERRIDE(
            void,
            Goldbeter1991OdeSystem,
            EvaluateYDerivatives,
                    time,
        rY,
        rDY);
    }

};
void register_Goldbeter1991OdeSystem_class(py::module &m){
py::class_<Goldbeter1991OdeSystem , Goldbeter1991OdeSystem_Overrides , boost::shared_ptr<Goldbeter1991OdeSystem >  , AbstractOdeSystem  >(m, "Goldbeter1991OdeSystem")
        .def(py::init<::std::vector<double> >(), py::arg("stateVariables") = std::vector<double>())
        .def(
            "EvaluateYDerivatives",
            (void(Goldbeter1991OdeSystem::*)(double, ::std::vector<double> const &, ::std::vector<double> &)) &Goldbeter1991OdeSystem::EvaluateYDerivatives,
            " " , py::arg("time"), py::arg("rY"), py::arg("rDY") )
    ;
}
