#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "AbstractCellPopulation.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ImmersedBoundaryBoundaryCellWriter.hpp"

#include "ImmersedBoundaryBoundaryCellWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryBoundaryCellWriter<2,2 > ImmersedBoundaryBoundaryCellWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundaryBoundaryCellWriter2_2_Overrides : public ImmersedBoundaryBoundaryCellWriter2_2{
    public:
    using ImmersedBoundaryBoundaryCellWriter2_2::ImmersedBoundaryBoundaryCellWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            ImmersedBoundaryBoundaryCellWriter2_2,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryBoundaryCellWriter2_2,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_ImmersedBoundaryBoundaryCellWriter2_2_class(py::module &m){
py::class_<ImmersedBoundaryBoundaryCellWriter2_2 , ImmersedBoundaryBoundaryCellWriter2_2_Overrides , boost::shared_ptr<ImmersedBoundaryBoundaryCellWriter2_2 >  , AbstractCellWriter<2, 2>  >(m, "ImmersedBoundaryBoundaryCellWriter2_2")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(ImmersedBoundaryBoundaryCellWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &ImmersedBoundaryBoundaryCellWriter2_2::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(ImmersedBoundaryBoundaryCellWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &ImmersedBoundaryBoundaryCellWriter2_2::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
