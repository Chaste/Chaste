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

#include "CellAncestorWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef CellAncestorWriter<3,3 > CellAncestorWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellAncestorWriter3_3_Overrides : public CellAncestorWriter3_3{
    public:
    using CellAncestorWriter3_3::CellAncestorWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellAncestorWriter3_3,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellAncestorWriter3_3,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellAncestorWriter3_3_class(py::module &m){
py::class_<CellAncestorWriter3_3 , CellAncestorWriter3_3_Overrides , boost::shared_ptr<CellAncestorWriter3_3 >  , AbstractCellWriter<3, 3>  >(m, "CellAncestorWriter3_3")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(CellAncestorWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellAncestorWriter3_3::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellAncestorWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellAncestorWriter3_3::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
