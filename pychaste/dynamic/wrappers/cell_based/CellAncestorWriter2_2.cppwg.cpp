#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "AbstractCellPopulation.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellAncestorWriter.hpp"

#include "CellAncestorWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef CellAncestorWriter<2,2 > CellAncestorWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellAncestorWriter2_2_Overrides : public CellAncestorWriter2_2{
    public:
    using CellAncestorWriter2_2::CellAncestorWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellAncestorWriter2_2,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellAncestorWriter2_2,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellAncestorWriter2_2_class(py::module &m){
py::class_<CellAncestorWriter2_2 , CellAncestorWriter2_2_Overrides , boost::shared_ptr<CellAncestorWriter2_2 >  , AbstractCellWriter<2, 2>  >(m, "CellAncestorWriter2_2")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(CellAncestorWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellAncestorWriter2_2::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellAncestorWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellAncestorWriter2_2::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
