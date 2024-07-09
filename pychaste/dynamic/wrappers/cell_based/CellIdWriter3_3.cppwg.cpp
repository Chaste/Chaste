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

#include "CellIdWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef CellIdWriter<3,3 > CellIdWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellIdWriter3_3_Overrides : public CellIdWriter3_3{
    public:
    using CellIdWriter3_3::CellIdWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellIdWriter3_3,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellIdWriter3_3,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellIdWriter3_3_class(py::module &m){
py::class_<CellIdWriter3_3 , CellIdWriter3_3_Overrides , boost::shared_ptr<CellIdWriter3_3 >  , AbstractCellWriter<3, 3>  >(m, "CellIdWriter3_3")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(CellIdWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellIdWriter3_3::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellIdWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellIdWriter3_3::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
