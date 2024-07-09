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
#include "AbstractGrowingDomainPdeModifier.hpp"

#include "AbstractGrowingDomainPdeModifier2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractGrowingDomainPdeModifier<2 > AbstractGrowingDomainPdeModifier2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractGrowingDomainPdeModifier2_Overrides : public AbstractGrowingDomainPdeModifier2{
    public:
    using AbstractGrowingDomainPdeModifier2::AbstractGrowingDomainPdeModifier;
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractGrowingDomainPdeModifier2,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_AbstractGrowingDomainPdeModifier2_class(py::module &m){
py::class_<AbstractGrowingDomainPdeModifier2 , AbstractGrowingDomainPdeModifier2_Overrides , boost::shared_ptr<AbstractGrowingDomainPdeModifier2 >  , AbstractPdeModifier<2>  >(m, "AbstractGrowingDomainPdeModifier2")
        .def(
            "GenerateFeMesh",
            (void(AbstractGrowingDomainPdeModifier2::*)(::AbstractCellPopulation<2> &)) &AbstractGrowingDomainPdeModifier2::GenerateFeMesh,
            " " , py::arg("rCellPopulation") )
        .def(
            "UpdateCellData",
            (void(AbstractGrowingDomainPdeModifier2::*)(::AbstractCellPopulation<2> &)) &AbstractGrowingDomainPdeModifier2::UpdateCellData,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputSimulationModifierParameters",
            (void(AbstractGrowingDomainPdeModifier2::*)(::out_stream &)) &AbstractGrowingDomainPdeModifier2::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
