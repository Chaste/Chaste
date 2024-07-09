#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <petsc/private/vecimpl.h>
#include <petsc/private/matimpl.h>
#include "PythonPetscObjectConverters.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractBoxDomainPdeModifier.hpp"

#include "AbstractBoxDomainPdeModifier3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractBoxDomainPdeModifier<3 > AbstractBoxDomainPdeModifier3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractBoxDomainPdeModifier3_Overrides : public AbstractBoxDomainPdeModifier3{
    public:
    using AbstractBoxDomainPdeModifier3::AbstractBoxDomainPdeModifier;
    void SetupSolve(::AbstractCellPopulation<3> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractBoxDomainPdeModifier3,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractBoxDomainPdeModifier3,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_AbstractBoxDomainPdeModifier3_class(py::module &m){
py::class_<AbstractBoxDomainPdeModifier3 , AbstractBoxDomainPdeModifier3_Overrides , boost::shared_ptr<AbstractBoxDomainPdeModifier3 >  , AbstractPdeModifier<3>  >(m, "AbstractBoxDomainPdeModifier3")
        .def(
            "GetStepSize",
            (double(AbstractBoxDomainPdeModifier3::*)()) &AbstractBoxDomainPdeModifier3::GetStepSize,
            " "  )
        .def(
            "SetBcsOnBoxBoundary",
            (void(AbstractBoxDomainPdeModifier3::*)(bool)) &AbstractBoxDomainPdeModifier3::SetBcsOnBoxBoundary,
            " " , py::arg("setBcsOnBoxBoundary") )
        .def(
            "AreBcsSetOnBoxBoundary",
            (bool(AbstractBoxDomainPdeModifier3::*)()) &AbstractBoxDomainPdeModifier3::AreBcsSetOnBoxBoundary,
            " "  )
        .def(
            "SetupSolve",
            (void(AbstractBoxDomainPdeModifier3::*)(::AbstractCellPopulation<3> &, ::std::string)) &AbstractBoxDomainPdeModifier3::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "GenerateFeMesh",
            (void(AbstractBoxDomainPdeModifier3::*)(::boost::shared_ptr<ChasteCuboid<3>>, double)) &AbstractBoxDomainPdeModifier3::GenerateFeMesh,
            " " , py::arg("pMeshCuboid"), py::arg("stepSize") )
        .def(
            "UpdateCellData",
            (void(AbstractBoxDomainPdeModifier3::*)(::AbstractCellPopulation<3> &)) &AbstractBoxDomainPdeModifier3::UpdateCellData,
            " " , py::arg("rCellPopulation") )
        .def(
            "InitialiseCellPdeElementMap",
            (void(AbstractBoxDomainPdeModifier3::*)(::AbstractCellPopulation<3> &)) &AbstractBoxDomainPdeModifier3::InitialiseCellPdeElementMap,
            " " , py::arg("rCellPopulation") )
        .def(
            "UpdateCellPdeElementMap",
            (void(AbstractBoxDomainPdeModifier3::*)(::AbstractCellPopulation<3> &)) &AbstractBoxDomainPdeModifier3::UpdateCellPdeElementMap,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputSimulationModifierParameters",
            (void(AbstractBoxDomainPdeModifier3::*)(::out_stream &)) &AbstractBoxDomainPdeModifier3::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
