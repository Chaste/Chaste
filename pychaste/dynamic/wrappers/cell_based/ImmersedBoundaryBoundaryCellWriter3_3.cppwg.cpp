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

#include "ImmersedBoundaryBoundaryCellWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryBoundaryCellWriter<3,3 > ImmersedBoundaryBoundaryCellWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundaryBoundaryCellWriter3_3_Overrides : public ImmersedBoundaryBoundaryCellWriter3_3{
    public:
    using ImmersedBoundaryBoundaryCellWriter3_3::ImmersedBoundaryBoundaryCellWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            ImmersedBoundaryBoundaryCellWriter3_3,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ImmersedBoundaryBoundaryCellWriter3_3,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_ImmersedBoundaryBoundaryCellWriter3_3_class(py::module &m){
py::class_<ImmersedBoundaryBoundaryCellWriter3_3 , ImmersedBoundaryBoundaryCellWriter3_3_Overrides , boost::shared_ptr<ImmersedBoundaryBoundaryCellWriter3_3 >  , AbstractCellWriter<3, 3>  >(m, "ImmersedBoundaryBoundaryCellWriter3_3")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(ImmersedBoundaryBoundaryCellWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &ImmersedBoundaryBoundaryCellWriter3_3::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(ImmersedBoundaryBoundaryCellWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &ImmersedBoundaryBoundaryCellWriter3_3::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
