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

#include "ImmersedBoundaryNeighbourNumberWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryNeighbourNumberWriter<3,3 > ImmersedBoundaryNeighbourNumberWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundaryNeighbourNumberWriter3_3_Overrides : public ImmersedBoundaryNeighbourNumberWriter3_3{
    public:
    using ImmersedBoundaryNeighbourNumberWriter3_3::ImmersedBoundaryNeighbourNumberWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            ImmersedBoundaryNeighbourNumberWriter3_3,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryNeighbourNumberWriter3_3,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_ImmersedBoundaryNeighbourNumberWriter3_3_class(py::module &m){
py::class_<ImmersedBoundaryNeighbourNumberWriter3_3 , ImmersedBoundaryNeighbourNumberWriter3_3_Overrides , boost::shared_ptr<ImmersedBoundaryNeighbourNumberWriter3_3 >  , AbstractCellWriter<3, 3>  >(m, "ImmersedBoundaryNeighbourNumberWriter3_3")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(ImmersedBoundaryNeighbourNumberWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &ImmersedBoundaryNeighbourNumberWriter3_3::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(ImmersedBoundaryNeighbourNumberWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &ImmersedBoundaryNeighbourNumberWriter3_3::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
