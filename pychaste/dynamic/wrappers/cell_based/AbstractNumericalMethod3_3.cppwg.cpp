#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractNumericalMethod.hpp"

#include "AbstractNumericalMethod3_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractNumericalMethod<3,3 > AbstractNumericalMethod3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractNumericalMethod3_3_Overrides : public AbstractNumericalMethod3_3{
    public:
    using AbstractNumericalMethod3_3::AbstractNumericalMethod;
    void UpdateAllNodePositions(double dt) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractNumericalMethod3_3,
            UpdateAllNodePositions,
                    dt);
    }
    void OutputNumericalMethodParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractNumericalMethod3_3,
            OutputNumericalMethodParameters,
                    rParamsFile);
    }

};
void register_AbstractNumericalMethod3_3_class(py::module &m){
py::class_<AbstractNumericalMethod3_3 , AbstractNumericalMethod3_3_Overrides , boost::shared_ptr<AbstractNumericalMethod3_3 >   >(m, "AbstractNumericalMethod3_3")
        .def(py::init< >())
        .def(
            "SetCellPopulation",
            (void(AbstractNumericalMethod3_3::*)(::AbstractOffLatticeCellPopulation<3> *)) &AbstractNumericalMethod3_3::SetCellPopulation,
            " " , py::arg("pPopulation") )
        .def(
            "SetForceCollection",
            (void(AbstractNumericalMethod3_3::*)(::std::vector<boost::shared_ptr<AbstractForce<3, 3>>> *)) &AbstractNumericalMethod3_3::SetForceCollection,
            " " , py::arg("pForces") )
        .def(
            "SetBoundaryConditions",
            (void(AbstractNumericalMethod3_3::*)(::std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<3, 3>>> *)) &AbstractNumericalMethod3_3::SetBoundaryConditions,
            " " , py::arg("pBoundaryConditions") )
        .def(
            "SetUseAdaptiveTimestep",
            (void(AbstractNumericalMethod3_3::*)(bool)) &AbstractNumericalMethod3_3::SetUseAdaptiveTimestep,
            " " , py::arg("useAdaptiveTimestep") )
        .def(
            "SetUseUpdateNodeLocation",
            (void(AbstractNumericalMethod3_3::*)(bool)) &AbstractNumericalMethod3_3::SetUseUpdateNodeLocation,
            " " , py::arg("useUpdateNodeLocation") )
        .def(
            "GetUseUpdateNodeLocation",
            (bool(AbstractNumericalMethod3_3::*)()) &AbstractNumericalMethod3_3::GetUseUpdateNodeLocation,
            " "  )
        .def(
            "HasAdaptiveTimestep",
            (bool(AbstractNumericalMethod3_3::*)()) &AbstractNumericalMethod3_3::HasAdaptiveTimestep,
            " "  )
        .def(
            "UpdateAllNodePositions",
            (void(AbstractNumericalMethod3_3::*)(double)) &AbstractNumericalMethod3_3::UpdateAllNodePositions,
            " " , py::arg("dt") )
        .def(
            "OutputNumericalMethodInfo",
            (void(AbstractNumericalMethod3_3::*)(::out_stream &)) &AbstractNumericalMethod3_3::OutputNumericalMethodInfo,
            " " , py::arg("rParamsFile") )
        .def(
            "OutputNumericalMethodParameters",
            (void(AbstractNumericalMethod3_3::*)(::out_stream &)) &AbstractNumericalMethod3_3::OutputNumericalMethodParameters,
            " " , py::arg("rParamsFile") )
    ;
}
