#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ShovingCaBasedDivisionRule.hpp"

#include "ShovingCaBasedDivisionRule3.cppwg.hpp"

namespace py = pybind11;
typedef ShovingCaBasedDivisionRule<3 > ShovingCaBasedDivisionRule3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;

class ShovingCaBasedDivisionRule3_Overrides : public ShovingCaBasedDivisionRule3{
    public:
    using ShovingCaBasedDivisionRule3::ShovingCaBasedDivisionRule;
    bool IsRoomToDivide(::CellPtr pParentCell, ::CaBasedCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            bool,
            ShovingCaBasedDivisionRule3,
            IsRoomToDivide,
                    pParentCell,
        rCellPopulation);
    }
    unsigned int CalculateDaughterNodeIndex(::CellPtr pNewCell, ::CellPtr pParentCell, ::CaBasedCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            unsignedint,
            ShovingCaBasedDivisionRule3,
            CalculateDaughterNodeIndex,
                    pNewCell,
        pParentCell,
        rCellPopulation);
    }

};
void register_ShovingCaBasedDivisionRule3_class(py::module &m){
py::class_<ShovingCaBasedDivisionRule3 , ShovingCaBasedDivisionRule3_Overrides , boost::shared_ptr<ShovingCaBasedDivisionRule3 >  , AbstractCaBasedDivisionRule<3>  >(m, "ShovingCaBasedDivisionRule3")
        .def(py::init< >())
        .def(
            "IsNodeOnBoundary",
            (void(ShovingCaBasedDivisionRule3::*)(unsigned int)) &ShovingCaBasedDivisionRule3::IsNodeOnBoundary,
            " " , py::arg("numNeighbours") )
        .def(
            "IsRoomToDivide",
            (bool(ShovingCaBasedDivisionRule3::*)(::CellPtr, ::CaBasedCellPopulation<3> &)) &ShovingCaBasedDivisionRule3::IsRoomToDivide,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
        .def(
            "CalculateDaughterNodeIndex",
            (unsigned int(ShovingCaBasedDivisionRule3::*)(::CellPtr, ::CellPtr, ::CaBasedCellPopulation<3> &)) &ShovingCaBasedDivisionRule3::CalculateDaughterNodeIndex,
            " " , py::arg("pNewCell"), py::arg("pParentCell"), py::arg("rCellPopulation") )
    ;
}
