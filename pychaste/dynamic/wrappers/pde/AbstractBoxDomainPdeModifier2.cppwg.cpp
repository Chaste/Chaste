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

#include "AbstractBoxDomainPdeModifier2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractBoxDomainPdeModifier<2 > AbstractBoxDomainPdeModifier2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractBoxDomainPdeModifier2_Overrides : public AbstractBoxDomainPdeModifier2{
    public:
    using AbstractBoxDomainPdeModifier2::AbstractBoxDomainPdeModifier;
    void SetupSolve(::AbstractCellPopulation<2> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractBoxDomainPdeModifier2,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractBoxDomainPdeModifier2,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_AbstractBoxDomainPdeModifier2_class(py::module &m){
py::class_<AbstractBoxDomainPdeModifier2 , AbstractBoxDomainPdeModifier2_Overrides , boost::shared_ptr<AbstractBoxDomainPdeModifier2 >  , AbstractPdeModifier<2>  >(m, "AbstractBoxDomainPdeModifier2")
        .def(
            "GetStepSize",
            (double(AbstractBoxDomainPdeModifier2::*)()) &AbstractBoxDomainPdeModifier2::GetStepSize,
            " "  )
        .def(
            "SetBcsOnBoxBoundary",
            (void(AbstractBoxDomainPdeModifier2::*)(bool)) &AbstractBoxDomainPdeModifier2::SetBcsOnBoxBoundary,
            " " , py::arg("setBcsOnBoxBoundary") )
        .def(
            "AreBcsSetOnBoxBoundary",
            (bool(AbstractBoxDomainPdeModifier2::*)()) &AbstractBoxDomainPdeModifier2::AreBcsSetOnBoxBoundary,
            " "  )
        .def(
            "SetupSolve",
            (void(AbstractBoxDomainPdeModifier2::*)(::AbstractCellPopulation<2> &, ::std::string)) &AbstractBoxDomainPdeModifier2::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "GenerateFeMesh",
            (void(AbstractBoxDomainPdeModifier2::*)(::boost::shared_ptr<ChasteCuboid<2>>, double)) &AbstractBoxDomainPdeModifier2::GenerateFeMesh,
            " " , py::arg("pMeshCuboid"), py::arg("stepSize") )
        .def(
            "UpdateCellData",
            (void(AbstractBoxDomainPdeModifier2::*)(::AbstractCellPopulation<2> &)) &AbstractBoxDomainPdeModifier2::UpdateCellData,
            " " , py::arg("rCellPopulation") )
        .def(
            "InitialiseCellPdeElementMap",
            (void(AbstractBoxDomainPdeModifier2::*)(::AbstractCellPopulation<2> &)) &AbstractBoxDomainPdeModifier2::InitialiseCellPdeElementMap,
            " " , py::arg("rCellPopulation") )
        .def(
            "UpdateCellPdeElementMap",
            (void(AbstractBoxDomainPdeModifier2::*)(::AbstractCellPopulation<2> &)) &AbstractBoxDomainPdeModifier2::UpdateCellPdeElementMap,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputSimulationModifierParameters",
            (void(AbstractBoxDomainPdeModifier2::*)(::out_stream &)) &AbstractBoxDomainPdeModifier2::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
