#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractCaBasedDivisionRule.hpp"

#include "AbstractCaBasedDivisionRule2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCaBasedDivisionRule<2 > AbstractCaBasedDivisionRule2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;

class AbstractCaBasedDivisionRule2_Overrides : public AbstractCaBasedDivisionRule2{
    public:
    using AbstractCaBasedDivisionRule2::AbstractCaBasedDivisionRule;
    bool IsRoomToDivide(::CellPtr pParentCell, ::CaBasedCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            bool,
            AbstractCaBasedDivisionRule2,
            IsRoomToDivide,
                    pParentCell,
        rCellPopulation);
    }
    unsigned int CalculateDaughterNodeIndex(::CellPtr pNewCell, ::CellPtr pParentCell, ::CaBasedCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            unsignedint,
            AbstractCaBasedDivisionRule2,
            CalculateDaughterNodeIndex,
                    pNewCell,
        pParentCell,
        rCellPopulation);
    }
    void OutputCellCaBasedDivisionRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCaBasedDivisionRule2,
            OutputCellCaBasedDivisionRuleParameters,
                    rParamsFile);
    }

};
void register_AbstractCaBasedDivisionRule2_class(py::module &m){
py::class_<AbstractCaBasedDivisionRule2 , AbstractCaBasedDivisionRule2_Overrides , boost::shared_ptr<AbstractCaBasedDivisionRule2 >   >(m, "AbstractCaBasedDivisionRule2")
        .def(py::init< >())
        .def(
            "IsRoomToDivide",
            (bool(AbstractCaBasedDivisionRule2::*)(::CellPtr, ::CaBasedCellPopulation<2> &)) &AbstractCaBasedDivisionRule2::IsRoomToDivide,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
        .def(
            "CalculateDaughterNodeIndex",
            (unsigned int(AbstractCaBasedDivisionRule2::*)(::CellPtr, ::CellPtr, ::CaBasedCellPopulation<2> &)) &AbstractCaBasedDivisionRule2::CalculateDaughterNodeIndex,
            " " , py::arg("pNewCell"), py::arg("pParentCell"), py::arg("rCellPopulation") )
        .def(
            "OutputCellCaBasedDivisionRuleInfo",
            (void(AbstractCaBasedDivisionRule2::*)(::out_stream &)) &AbstractCaBasedDivisionRule2::OutputCellCaBasedDivisionRuleInfo,
            " " , py::arg("rParamsFile") )
    ;
}
