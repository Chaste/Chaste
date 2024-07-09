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
#include "EllipticGrowingDomainPdeModifier.hpp"

#include "EllipticGrowingDomainPdeModifier3.cppwg.hpp"

namespace py = pybind11;
typedef EllipticGrowingDomainPdeModifier<3 > EllipticGrowingDomainPdeModifier3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::std::shared_ptr<BoundaryConditionsContainer<3, 3, 1>> _std_shared_ptr_lt_BoundaryConditionsContainer_lt_3_3_1_gt__gt_;

class EllipticGrowingDomainPdeModifier3_Overrides : public EllipticGrowingDomainPdeModifier3{
    public:
    using EllipticGrowingDomainPdeModifier3::EllipticGrowingDomainPdeModifier;
    void UpdateAtEndOfTimeStep(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            EllipticGrowingDomainPdeModifier3,
            UpdateAtEndOfTimeStep,
                    rCellPopulation);
    }
    void SetupSolve(::AbstractCellPopulation<3> & rCellPopulation, ::std::string outputDirectory) override {
        PYBIND11_OVERRIDE(
            void,
            EllipticGrowingDomainPdeModifier3,
            SetupSolve,
                    rCellPopulation,
        outputDirectory);
    }
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            EllipticGrowingDomainPdeModifier3,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_EllipticGrowingDomainPdeModifier3_class(py::module &m){
py::class_<EllipticGrowingDomainPdeModifier3 , EllipticGrowingDomainPdeModifier3_Overrides , boost::shared_ptr<EllipticGrowingDomainPdeModifier3 >  , AbstractGrowingDomainPdeModifier<3>  >(m, "EllipticGrowingDomainPdeModifier3")
        .def(py::init<::boost::shared_ptr<AbstractLinearPde<3, 3>>, ::boost::shared_ptr<AbstractBoundaryCondition<3>>, bool, ::Vec >(), py::arg("pPde") = boost::shared_ptr<AbstractLinearPde<3, 3>>(), py::arg("pBoundaryCondition") = boost::shared_ptr<AbstractBoundaryCondition<3>>(), py::arg("isNeumannBoundaryCondition") = true, py::arg("solution") = nullptr)
        .def(
            "UpdateAtEndOfTimeStep",
            (void(EllipticGrowingDomainPdeModifier3::*)(::AbstractCellPopulation<3> &)) &EllipticGrowingDomainPdeModifier3::UpdateAtEndOfTimeStep,
            " " , py::arg("rCellPopulation") )
        .def(
            "SetupSolve",
            (void(EllipticGrowingDomainPdeModifier3::*)(::AbstractCellPopulation<3> &, ::std::string)) &EllipticGrowingDomainPdeModifier3::SetupSolve,
            " " , py::arg("rCellPopulation"), py::arg("outputDirectory") )
        .def(
            "OutputSimulationModifierParameters",
            (void(EllipticGrowingDomainPdeModifier3::*)(::out_stream &)) &EllipticGrowingDomainPdeModifier3::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
