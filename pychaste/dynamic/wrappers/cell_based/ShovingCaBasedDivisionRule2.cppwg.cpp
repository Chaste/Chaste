#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ShovingCaBasedDivisionRule.hpp"

#include "ShovingCaBasedDivisionRule2.cppwg.hpp"

namespace py = pybind11;
typedef ShovingCaBasedDivisionRule<2 > ShovingCaBasedDivisionRule2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;

class ShovingCaBasedDivisionRule2_Overrides : public ShovingCaBasedDivisionRule2{
    public:
    using ShovingCaBasedDivisionRule2::ShovingCaBasedDivisionRule;
    bool IsRoomToDivide(::CellPtr pParentCell, ::CaBasedCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            bool,
            ShovingCaBasedDivisionRule2,
            IsRoomToDivide,
                    pParentCell,
        rCellPopulation);
    }
    unsigned int CalculateDaughterNodeIndex(::CellPtr pNewCell, ::CellPtr pParentCell, ::CaBasedCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            unsignedint,
            ShovingCaBasedDivisionRule2,
            CalculateDaughterNodeIndex,
                    pNewCell,
        pParentCell,
        rCellPopulation);
    }

};
void register_ShovingCaBasedDivisionRule2_class(py::module &m){
py::class_<ShovingCaBasedDivisionRule2 , ShovingCaBasedDivisionRule2_Overrides , boost::shared_ptr<ShovingCaBasedDivisionRule2 >  , AbstractCaBasedDivisionRule<2>  >(m, "ShovingCaBasedDivisionRule2")
        .def(py::init< >())
        .def(
            "IsNodeOnBoundary",
            (void(ShovingCaBasedDivisionRule2::*)(unsigned int)) &ShovingCaBasedDivisionRule2::IsNodeOnBoundary,
            " " , py::arg("numNeighbours") )
        .def(
            "IsRoomToDivide",
            (bool(ShovingCaBasedDivisionRule2::*)(::CellPtr, ::CaBasedCellPopulation<2> &)) &ShovingCaBasedDivisionRule2::IsRoomToDivide,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
        .def(
            "CalculateDaughterNodeIndex",
            (unsigned int(ShovingCaBasedDivisionRule2::*)(::CellPtr, ::CellPtr, ::CaBasedCellPopulation<2> &)) &ShovingCaBasedDivisionRule2::CalculateDaughterNodeIndex,
            " " , py::arg("pNewCell"), py::arg("pParentCell"), py::arg("rCellPopulation") )
    ;
}
