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

#include "ParabolicBoxDomainPdeModifier2.cppwg.hpp"

namespace py = pybind11;
typedef ParabolicBoxDomainPdeModifier<2 > ParabolicBoxDomainPdeModifier2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::std::shared_ptr<BoundaryConditionsContainer<2, 2, 1>> _std_shared_ptr_lt_BoundaryConditionsContainer_lt_2_2_1_gt__gt_;

class ParabolicBoxDomainPdeModifier2_Overrides : public ParabolicBoxDomainPdeModifier2{
    public:
    using ParabolicBoxDomainPdeModifier2::ParabolicBoxDomainPdeModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ParabolicBoxDomainPdeModifier2,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<2> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            ParabolicBoxDomainPdeModifier2,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ParabolicBoxDomainPdeModifier2,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_ParabolicBoxDomainPdeModifier2_class(py::module &m){
py::class_<ParabolicBoxDomainPdeModifier2 , ParabolicBoxDomainPdeModifier2_Overrides , boost::shared_ptr<ParabolicBoxDomainPdeModifier2 >  , AbstractBoxDomainPdeModifier<2>  >(m, "ParabolicBoxDomainPdeModifier2")
        .def(py::init<::boost::shared_ptr<AbstractLinearPde<2, 2>>, ::boost::shared_ptr<AbstractBoundaryCondition<2>>, bool, ::boost::shared_ptr<ChasteCuboid<2>>, double, ::Vec >(), py::arg("pPde") = boost::shared_ptr<AbstractLinearPde<2, 2>>(), py::arg("pBoundaryCondition") = boost::shared_ptr<AbstractBoundaryCondition<2>>(), py::arg("isNeumannBoundaryCondition") = true, py::arg("pMeshCuboid") = boost::shared_ptr<ChasteCuboid<2>>(), py::arg("stepSize") = 1., py::arg("solution") = nullptr)
        .def(
            "UpdateAtEndOfTimeStep",
            (void(ParabolicBoxDomainPdeModifier2::*)(::AbstractCellPopulation<2> &)) &ParabolicBoxDomainPdeModifier2::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(ParabolicBoxDomainPdeModifier2::*)(::AbstractCellPopulation<2> &, ::std::string)) &ParabolicBoxDomainPdeModifier2::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "SetupInitialSolutionVector",
            (void(ParabolicBoxDomainPdeModifier2::*)(::AbstractCellPopulation<2> &)) &ParabolicBoxDomainPdeModifier2::SetupInitialSolutionVector,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputSimulationModifierParameters",
            (void(ParabolicBoxDomainPdeModifier2::*)(::out_stream &)) &ParabolicBoxDomainPdeModifier2::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
