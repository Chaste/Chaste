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

#include "AbstractGrowingDomainPdeModifier3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractGrowingDomainPdeModifier<3 > AbstractGrowingDomainPdeModifier3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractGrowingDomainPdeModifier3_Overrides : public AbstractGrowingDomainPdeModifier3{
    public:
    using AbstractGrowingDomainPdeModifier3::AbstractGrowingDomainPdeModifier;
    void OutputSimulationModifierParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractGrowingDomainPdeModifier3,
            OutputSimulationModifierParameters,
                    rParamsFile);
    }

};
void register_AbstractGrowingDomainPdeModifier3_class(py::module &m){
py::class_<AbstractGrowingDomainPdeModifier3 , AbstractGrowingDomainPdeModifier3_Overrides , boost::shared_ptr<AbstractGrowingDomainPdeModifier3 >  , AbstractPdeModifier<3>  >(m, "AbstractGrowingDomainPdeModifier3")
        .def(
            "GenerateFeMesh",
            (void(AbstractGrowingDomainPdeModifier3::*)(::AbstractCellPopulation<3> &)) &AbstractGrowingDomainPdeModifier3::GenerateFeMesh,
            " " , py::arg("rCellPopulation") )
        .def(
            "UpdateCellData",
            (void(AbstractGrowingDomainPdeModifier3::*)(::AbstractCellPopulation<3> &)) &AbstractGrowingDomainPdeModifier3::UpdateCellData,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputSimulationModifierParameters",
            (void(AbstractGrowingDomainPdeModifier3::*)(::out_stream &)) &AbstractGrowingDomainPdeModifier3::OutputSimulationModifierParameters,
            " " , py::arg("rParamsFile") )
    ;
}
