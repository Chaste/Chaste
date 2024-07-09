#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ExclusionCaBasedDivisionRule.hpp"

#include "ExclusionCaBasedDivisionRule3.cppwg.hpp"

namespace py = pybind11;
typedef ExclusionCaBasedDivisionRule<3 > ExclusionCaBasedDivisionRule3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;

class ExclusionCaBasedDivisionRule3_Overrides : public ExclusionCaBasedDivisionRule3{
    public:
    using ExclusionCaBasedDivisionRule3::ExclusionCaBasedDivisionRule;
    bool IsRoomToDivide(::CellPtr pParentCell, ::CaBasedCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            bool,
            ExclusionCaBasedDivisionRule3,
            IsRoomToDivide,
                    pParentCell,
        rCellPopulation);
    }
    unsigned int CalculateDaughterNodeIndex(::CellPtr pNewCell, ::CellPtr pParentCell, ::CaBasedCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            unsignedint,
            ExclusionCaBasedDivisionRule3,
            CalculateDaughterNodeIndex,
                    pNewCell,
        pParentCell,
        rCellPopulation);
    }

};
void register_ExclusionCaBasedDivisionRule3_class(py::module &m){
py::class_<ExclusionCaBasedDivisionRule3 , ExclusionCaBasedDivisionRule3_Overrides , boost::shared_ptr<ExclusionCaBasedDivisionRule3 >  , AbstractCaBasedDivisionRule<3>  >(m, "ExclusionCaBasedDivisionRule3")
        .def(py::init< >())
        .def(
            "IsRoomToDivide",
            (bool(ExclusionCaBasedDivisionRule3::*)(::CellPtr, ::CaBasedCellPopulation<3> &)) &ExclusionCaBasedDivisionRule3::IsRoomToDivide,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
        .def(
            "CalculateDaughterNodeIndex",
            (unsigned int(ExclusionCaBasedDivisionRule3::*)(::CellPtr, ::CellPtr, ::CaBasedCellPopulation<3> &)) &ExclusionCaBasedDivisionRule3::CalculateDaughterNodeIndex,
            " " , py::arg("pNewCell"), py::arg("pParentCell"), py::arg("rCellPopulation") )
    ;
}
