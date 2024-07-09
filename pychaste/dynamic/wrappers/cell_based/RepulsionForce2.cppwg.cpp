#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "RepulsionForce.hpp"

#include "RepulsionForce2.cppwg.hpp"

namespace py = pybind11;
typedef RepulsionForce<2 > RepulsionForce2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class RepulsionForce2_Overrides : public RepulsionForce2{
    public:
    using RepulsionForce2::RepulsionForce;
    void AddForceContribution(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            RepulsionForce2,
            AddForceContribution,
                    rCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            RepulsionForce2,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_RepulsionForce2_class(py::module &m){
py::class_<RepulsionForce2 , RepulsionForce2_Overrides , boost::shared_ptr<RepulsionForce2 >  , GeneralisedLinearSpringForce<2>  >(m, "RepulsionForce2")
        .def(py::init< >())
        .def(
            "AddForceContribution",
            (void(RepulsionForce2::*)(::AbstractCellPopulation<2> &)) &RepulsionForce2::AddForceContribution,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputForceParameters",
            (void(RepulsionForce2::*)(::out_stream &)) &RepulsionForce2::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
