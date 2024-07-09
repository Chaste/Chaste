#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ExclusionCaBasedDivisionRule.hpp"

#include "ExclusionCaBasedDivisionRule2.cppwg.hpp"

namespace py = pybind11;
typedef ExclusionCaBasedDivisionRule<2 > ExclusionCaBasedDivisionRule2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;

class ExclusionCaBasedDivisionRule2_Overrides : public ExclusionCaBasedDivisionRule2{
    public:
    using ExclusionCaBasedDivisionRule2::ExclusionCaBasedDivisionRule;
    bool IsRoomToDivide(::CellPtr pParentCell, ::CaBasedCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            bool,
            ExclusionCaBasedDivisionRule2,
            IsRoomToDivide,
                    pParentCell,
        rCellPopulation);
    }
    unsigned int CalculateDaughterNodeIndex(::CellPtr pNewCell, ::CellPtr pParentCell, ::CaBasedCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            unsignedint,
            ExclusionCaBasedDivisionRule2,
            CalculateDaughterNodeIndex,
                    pNewCell,
        pParentCell,
        rCellPopulation);
    }

};
void register_ExclusionCaBasedDivisionRule2_class(py::module &m){
py::class_<ExclusionCaBasedDivisionRule2 , ExclusionCaBasedDivisionRule2_Overrides , boost::shared_ptr<ExclusionCaBasedDivisionRule2 >  , AbstractCaBasedDivisionRule<2>  >(m, "ExclusionCaBasedDivisionRule2")
        .def(py::init< >())
        .def(
            "IsRoomToDivide",
            (bool(ExclusionCaBasedDivisionRule2::*)(::CellPtr, ::CaBasedCellPopulation<2> &)) &ExclusionCaBasedDivisionRule2::IsRoomToDivide,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
        .def(
            "CalculateDaughterNodeIndex",
            (unsigned int(ExclusionCaBasedDivisionRule2::*)(::CellPtr, ::CellPtr, ::CaBasedCellPopulation<2> &)) &ExclusionCaBasedDivisionRule2::CalculateDaughterNodeIndex,
            " " , py::arg("pNewCell"), py::arg("pParentCell"), py::arg("rCellPopulation") )
    ;
}
