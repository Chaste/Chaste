#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "DeltaNotchEdgeOdeSystem.hpp"

#include "DeltaNotchEdgeOdeSystem.cppwg.hpp"

namespace py = pybind11;
typedef DeltaNotchEdgeOdeSystem DeltaNotchEdgeOdeSystem;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class DeltaNotchEdgeOdeSystem_Overrides : public DeltaNotchEdgeOdeSystem{
    public:
    using DeltaNotchEdgeOdeSystem::DeltaNotchEdgeOdeSystem;
    void EvaluateYDerivatives(double time, ::std::vector<double> const & rY, ::std::vector<double> & rDY) override {
        PYBIND11_OVERRIDE(
            void,
            DeltaNotchEdgeOdeSystem,
            EvaluateYDerivatives,
                    time,
        rY,
        rDY);
    }

};
void register_DeltaNotchEdgeOdeSystem_class(py::module &m){
py::class_<DeltaNotchEdgeOdeSystem , DeltaNotchEdgeOdeSystem_Overrides , boost::shared_ptr<DeltaNotchEdgeOdeSystem >  , AbstractOdeSystem  >(m, "DeltaNotchEdgeOdeSystem")
        .def(py::init<::std::vector<double> >(), py::arg("stateVariables") = std::vector<double>())
        .def(
            "EvaluateYDerivatives",
            (void(DeltaNotchEdgeOdeSystem::*)(double, ::std::vector<double> const &, ::std::vector<double> &)) &DeltaNotchEdgeOdeSystem::EvaluateYDerivatives,
            " " , py::arg("time"), py::arg("rY"), py::arg("rDY") )
    ;
}
