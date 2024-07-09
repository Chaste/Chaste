#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "AbstractCellPopulation.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellIdWriter.hpp"

#include "CellIdWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef CellIdWriter<2,2 > CellIdWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellIdWriter2_2_Overrides : public CellIdWriter2_2{
    public:
    using CellIdWriter2_2::CellIdWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellIdWriter2_2,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellIdWriter2_2,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellIdWriter2_2_class(py::module &m){
py::class_<CellIdWriter2_2 , CellIdWriter2_2_Overrides , boost::shared_ptr<CellIdWriter2_2 >  , AbstractCellWriter<2, 2>  >(m, "CellIdWriter2_2")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(CellIdWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellIdWriter2_2::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellIdWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellIdWriter2_2::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
