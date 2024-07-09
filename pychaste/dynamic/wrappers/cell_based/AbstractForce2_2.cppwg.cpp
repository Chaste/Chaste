#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractForce.hpp"

#include "AbstractForce2_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractForce<2,2 > AbstractForce2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractForce2_2_Overrides : public AbstractForce2_2{
    public:
    using AbstractForce2_2::AbstractForce;
    void AddForceContribution(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractForce2_2,
            AddForceContribution,
                    rCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractForce2_2,
            OutputForceParameters,
                    rParamsFile);
    }
    void WriteDataToVisualizerSetupFile(::out_stream & pVizSetupFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractForce2_2,
            WriteDataToVisualizerSetupFile,
                    pVizSetupFile);
    }

};
void register_AbstractForce2_2_class(py::module &m){
py::class_<AbstractForce2_2 , AbstractForce2_2_Overrides , boost::shared_ptr<AbstractForce2_2 >   >(m, "AbstractForce2_2")
        .def(py::init< >())
        .def(
            "AddForceContribution",
            (void(AbstractForce2_2::*)(::AbstractCellPopulation<2> &)) &AbstractForce2_2::AddForceContribution,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputForceInfo",
            (void(AbstractForce2_2::*)(::out_stream &)) &AbstractForce2_2::OutputForceInfo,
            " " , py::arg("rParamsFile") )
        .def(
            "OutputForceParameters",
            (void(AbstractForce2_2::*)(::out_stream &)) &AbstractForce2_2::OutputForceParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "WriteDataToVisualizerSetupFile",
            (void(AbstractForce2_2::*)(::out_stream &)) &AbstractForce2_2::WriteDataToVisualizerSetupFile,
            " " , py::arg("pVizSetupFile") )
    ;
}
