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
#include "ParabolicGrowingDomainPdeModifier.hpp"

#include "ParabolicGrowingDomainPdeModifier2.cppwg.hpp"

namespace py = pybind11;
typedef ParabolicGrowingDomainPdeModifier<2 > ParabolicGrowingDomainPdeModifier2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::std::shared_ptr<BoundaryConditionsContainer<2, 2, 1>> _std_shared_ptr_lt_BoundaryConditionsContainer_lt_2_2_1_gt__gt_;

class ParabolicGrowingDomainPdeModifier2_Overrides : public ParabolicGrowingDomainPdeModifier2{
    public:
    using ParabolicGrowingDomainPdeModifier2::ParabolicGrowingDomainPdeModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ParabolicGrowingDomainPdeModifier2,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<2> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            ParabolicGrowingDomainPdeModifier2,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ParabolicGrowingDomainPdeModifier2,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_ParabolicGrowingDomainPdeModifier2_class(py::module &m){
py::class_<ParabolicGrowingDomainPdeModifier2 , ParabolicGrowingDomainPdeModifier2_Overrides , boost::shared_ptr<ParabolicGrowingDomainPdeModifier2 >  , AbstractGrowingDomainPdeModifier<2>  >(m, "ParabolicGrowingDomainPdeModifier2")
        .def(py::init<::boost::shared_ptr<AbstractLinearPde<2, 2>>, ::boost::shared_ptr<AbstractBoundaryCondition<2>>, bool, ::Vec >(), py::arg("pPde") = boost::shared_ptr<AbstractLinearPde<2, 2>>(), py::arg("pBoundaryCondition") = boost::shared_ptr<AbstractBoundaryCondition<2>>(), py::arg("isNeumannBoundaryCondition") = true, py::arg("solution") = nullptr)
        .def(
            "UpdateAtEndOfTimeStep",
            (void(ParabolicGrowingDomainPdeModifier2::*)(::AbstractCellPopulation<2> &)) &ParabolicGrowingDomainPdeModifier2::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(ParabolicGrowingDomainPdeModifier2::*)(::AbstractCellPopulation<2> &, ::std::string)) &ParabolicGrowingDomainPdeModifier2::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "UpdateSolutionVector",
            (void(ParabolicGrowingDomainPdeModifier2::*)(::AbstractCellPopulation<2> &)) &ParabolicGrowingDomainPdeModifier2::UpdateSolutionVector,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputSimulationModifierParameters",
            (void(ParabolicGrowingDomainPdeModifier2::*)(::out_stream &)) &ParabolicGrowingDomainPdeModifier2::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
