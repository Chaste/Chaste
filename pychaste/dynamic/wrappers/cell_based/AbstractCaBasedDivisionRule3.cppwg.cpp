#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractCaBasedDivisionRule.hpp"

#include "AbstractCaBasedDivisionRule3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCaBasedDivisionRule<3 > AbstractCaBasedDivisionRule3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;

class AbstractCaBasedDivisionRule3_Overrides : public AbstractCaBasedDivisionRule3{
    public:
    using AbstractCaBasedDivisionRule3::AbstractCaBasedDivisionRule;
    bool IsRoomToDivide(::CellPtr pParentCell, ::CaBasedCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            bool,
            AbstractCaBasedDivisionRule3,
            IsRoomToDivide,
                    pParentCell,
        rCellPopulation);
    }
    unsigned int CalculateDaughterNodeIndex(::CellPtr pNewCell, ::CellPtr pParentCell, ::CaBasedCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            unsignedint,
            AbstractCaBasedDivisionRule3,
            CalculateDaughterNodeIndex,
                    pNewCell,
        pParentCell,
        rCellPopulation);
    }
    void OutputCellCaBasedDivisionRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCaBasedDivisionRule3,
            OutputCellCaBasedDivisionRuleParameters,
                    rParamsFile);
    }

};
void register_AbstractCaBasedDivisionRule3_class(py::module &m){
py::class_<AbstractCaBasedDivisionRule3 , AbstractCaBasedDivisionRule3_Overrides , boost::shared_ptr<AbstractCaBasedDivisionRule3 >   >(m, "AbstractCaBasedDivisionRule3")
        .def(py::init< >())
        .def(
            "IsRoomToDivide",
            (bool(AbstractCaBasedDivisionRule3::*)(::CellPtr, ::CaBasedCellPopulation<3> &)) &AbstractCaBasedDivisionRule3::IsRoomToDivide,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
        .def(
            "CalculateDaughterNodeIndex",
            (unsigned int(AbstractCaBasedDivisionRule3::*)(::CellPtr, ::CellPtr, ::CaBasedCellPopulation<3> &)) &AbstractCaBasedDivisionRule3::CalculateDaughterNodeIndex,
            " " , py::arg("pNewCell"), py::arg("pParentCell"), py::arg("rCellPopulation") )
        .def(
            "OutputCellCaBasedDivisionRuleInfo",
            (void(AbstractCaBasedDivisionRule3::*)(::out_stream &)) &AbstractCaBasedDivisionRule3::OutputCellCaBasedDivisionRuleInfo,
            " " , py::arg("rParamsFile") )
    ;
}
