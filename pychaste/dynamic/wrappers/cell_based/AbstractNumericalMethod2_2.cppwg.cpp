#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractNumericalMethod.hpp"

#include "AbstractNumericalMethod2_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractNumericalMethod<2,2 > AbstractNumericalMethod2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractNumericalMethod2_2_Overrides : public AbstractNumericalMethod2_2{
    public:
    using AbstractNumericalMethod2_2::AbstractNumericalMethod;
    void UpdateAllNodePositions(double dt) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractNumericalMethod2_2,
            UpdateAllNodePositions,
                    dt);
    }
    void OutputNumericalMethodParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractNumericalMethod2_2,
            OutputNumericalMethodParameters,
                    rParamsFile);
    }

};
void register_AbstractNumericalMethod2_2_class(py::module &m){
py::class_<AbstractNumericalMethod2_2 , AbstractNumericalMethod2_2_Overrides , boost::shared_ptr<AbstractNumericalMethod2_2 >   >(m, "AbstractNumericalMethod2_2")
        .def(py::init< >())
        .def(
            "SetCellPopulation",
            (void(AbstractNumericalMethod2_2::*)(::AbstractOffLatticeCellPopulation<2> *)) &AbstractNumericalMethod2_2::SetCellPopulation,
            " " , py::arg("pPopulation") )
        .def(
            "SetForceCollection",
            (void(AbstractNumericalMethod2_2::*)(::std::vector<boost::shared_ptr<AbstractForce<2, 2>>> *)) &AbstractNumericalMethod2_2::SetForceCollection,
            " " , py::arg("pForces") )
        .def(
            "SetBoundaryConditions",
            (void(AbstractNumericalMethod2_2::*)(::std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<2, 2>>> *)) &AbstractNumericalMethod2_2::SetBoundaryConditions,
            " " , py::arg("pBoundaryConditions") )
        .def(
            "SetUseAdaptiveTimestep",
            (void(AbstractNumericalMethod2_2::*)(bool)) &AbstractNumericalMethod2_2::SetUseAdaptiveTimestep,
            " " , py::arg("useAdaptiveTimestep") )
        .def(
            "SetUseUpdateNodeLocation",
            (void(AbstractNumericalMethod2_2::*)(bool)) &AbstractNumericalMethod2_2::SetUseUpdateNodeLocation,
            " " , py::arg("useUpdateNodeLocation") )
        .def(
            "GetUseUpdateNodeLocation",
            (bool(AbstractNumericalMethod2_2::*)()) &AbstractNumericalMethod2_2::GetUseUpdateNodeLocation,
            " "  )
        .def(
            "HasAdaptiveTimestep",
            (bool(AbstractNumericalMethod2_2::*)()) &AbstractNumericalMethod2_2::HasAdaptiveTimestep,
            " "  )
        .def(
            "UpdateAllNodePositions",
            (void(AbstractNumericalMethod2_2::*)(double)) &AbstractNumericalMethod2_2::UpdateAllNodePositions,
            " " , py::arg("dt") )
        .def(
            "OutputNumericalMethodInfo",
            (void(AbstractNumericalMethod2_2::*)(::out_stream &)) &AbstractNumericalMethod2_2::OutputNumericalMethodInfo,
            " " , py::arg("rParamsFile") )
        .def(
            "OutputNumericalMethodParameters",
            (void(AbstractNumericalMethod2_2::*)(::out_stream &)) &AbstractNumericalMethod2_2::OutputNumericalMethodParameters,
            " " , py::arg("rParamsFile") )
    ;
}
