#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "AbstractCellPopulation.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellAgesWriter.hpp"

#include "CellAgesWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef CellAgesWriter<2,2 > CellAgesWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellAgesWriter2_2_Overrides : public CellAgesWriter2_2{
    public:
    using CellAgesWriter2_2::CellAgesWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellAgesWriter2_2,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellAgesWriter2_2,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellAgesWriter2_2_class(py::module &m){
py::class_<CellAgesWriter2_2 , CellAgesWriter2_2_Overrides , boost::shared_ptr<CellAgesWriter2_2 >  , AbstractCellWriter<2, 2>  >(m, "CellAgesWriter2_2")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(CellAgesWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellAgesWriter2_2::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellAgesWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellAgesWriter2_2::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
