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
#include "ImmersedBoundaryNeighbourNumberWriter.hpp"

#include "ImmersedBoundaryNeighbourNumberWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryNeighbourNumberWriter<2,2 > ImmersedBoundaryNeighbourNumberWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundaryNeighbourNumberWriter2_2_Overrides : public ImmersedBoundaryNeighbourNumberWriter2_2{
    public:
    using ImmersedBoundaryNeighbourNumberWriter2_2::ImmersedBoundaryNeighbourNumberWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            ImmersedBoundaryNeighbourNumberWriter2_2,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryNeighbourNumberWriter2_2,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_ImmersedBoundaryNeighbourNumberWriter2_2_class(py::module &m){
py::class_<ImmersedBoundaryNeighbourNumberWriter2_2 , ImmersedBoundaryNeighbourNumberWriter2_2_Overrides , boost::shared_ptr<ImmersedBoundaryNeighbourNumberWriter2_2 >  , AbstractCellWriter<2, 2>  >(m, "ImmersedBoundaryNeighbourNumberWriter2_2")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(ImmersedBoundaryNeighbourNumberWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &ImmersedBoundaryNeighbourNumberWriter2_2::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(ImmersedBoundaryNeighbourNumberWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &ImmersedBoundaryNeighbourNumberWriter2_2::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
