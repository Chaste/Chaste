#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "RepulsionForce.hpp"

#include "RepulsionForce3.cppwg.hpp"

namespace py = pybind11;
typedef RepulsionForce<3 > RepulsionForce3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class RepulsionForce3_Overrides : public RepulsionForce3{
    public:
    using RepulsionForce3::RepulsionForce;
    void AddForceContribution(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            RepulsionForce3,
            AddForceContribution,
                    rCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            RepulsionForce3,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_RepulsionForce3_class(py::module &m){
py::class_<RepulsionForce3 , RepulsionForce3_Overrides , boost::shared_ptr<RepulsionForce3 >  , GeneralisedLinearSpringForce<3>  >(m, "RepulsionForce3")
        .def(py::init< >())
        .def(
            "AddForceContribution",
            (void(RepulsionForce3::*)(::AbstractCellPopulation<3> &)) &RepulsionForce3::AddForceContribution,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputForceParameters",
            (void(RepulsionForce3::*)(::out_stream &)) &RepulsionForce3::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
