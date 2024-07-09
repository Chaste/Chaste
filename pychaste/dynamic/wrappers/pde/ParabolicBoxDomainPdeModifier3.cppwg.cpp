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
#include "ParabolicBoxDomainPdeModifier.hpp"

#include "ParabolicBoxDomainPdeModifier3.cppwg.hpp"

namespace py = pybind11;
typedef ParabolicBoxDomainPdeModifier<3 > ParabolicBoxDomainPdeModifier3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::std::shared_ptr<BoundaryConditionsContainer<3, 3, 1>> _std_shared_ptr_lt_BoundaryConditionsContainer_lt_3_3_1_gt__gt_;

class ParabolicBoxDomainPdeModifier3_Overrides : public ParabolicBoxDomainPdeModifier3{
    public:
    using ParabolicBoxDomainPdeModifier3::ParabolicBoxDomainPdeModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ParabolicBoxDomainPdeModifier3,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<3> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            ParabolicBoxDomainPdeModifier3,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ParabolicBoxDomainPdeModifier3,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_ParabolicBoxDomainPdeModifier3_class(py::module &m){
py::class_<ParabolicBoxDomainPdeModifier3 , ParabolicBoxDomainPdeModifier3_Overrides , boost::shared_ptr<ParabolicBoxDomainPdeModifier3 >  , AbstractBoxDomainPdeModifier<3>  >(m, "ParabolicBoxDomainPdeModifier3")
        .def(py::init<::boost::shared_ptr<AbstractLinearPde<3, 3>>, ::boost::shared_ptr<AbstractBoundaryCondition<3>>, bool, ::boost::shared_ptr<ChasteCuboid<3>>, double, ::Vec >(), py::arg("pPde") = boost::shared_ptr<AbstractLinearPde<3, 3>>(), py::arg("pBoundaryCondition") = boost::shared_ptr<AbstractBoundaryCondition<3>>(), py::arg("isNeumannBoundaryCondition") = true, py::arg("pMeshCuboid") = boost::shared_ptr<ChasteCuboid<3>>(), py::arg("stepSize") = 1., py::arg("solution") = nullptr)
        .def(
            "UpdateAtEndOfTimeStep",
            (void(ParabolicBoxDomainPdeModifier3::*)(::AbstractCellPopulation<3> &)) &ParabolicBoxDomainPdeModifier3::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(ParabolicBoxDomainPdeModifier3::*)(::AbstractCellPopulation<3> &, ::std::string)) &ParabolicBoxDomainPdeModifier3::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "SetupInitialSolutionVector",
            (void(ParabolicBoxDomainPdeModifier3::*)(::AbstractCellPopulation<3> &)) &ParabolicBoxDomainPdeModifier3::SetupInitialSolutionVector,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputSimulationModifierParameters",
            (void(ParabolicBoxDomainPdeModifier3::*)(::out_stream &)) &ParabolicBoxDomainPdeModifier3::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
