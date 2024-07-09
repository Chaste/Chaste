#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "AbstractCellPopulation.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellLocationIndexWriter.hpp"

#include "CellLocationIndexWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef CellLocationIndexWriter<2,2 > CellLocationIndexWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellLocationIndexWriter2_2_Overrides : public CellLocationIndexWriter2_2{
    public:
    using CellLocationIndexWriter2_2::CellLocationIndexWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellLocationIndexWriter2_2,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellLocationIndexWriter2_2,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellLocationIndexWriter2_2_class(py::module &m){
py::class_<CellLocationIndexWriter2_2 , CellLocationIndexWriter2_2_Overrides , boost::shared_ptr<CellLocationIndexWriter2_2 >  , AbstractCellWriter<2, 2>  >(m, "CellLocationIndexWriter2_2")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(CellLocationIndexWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellLocationIndexWriter2_2::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellLocationIndexWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellLocationIndexWriter2_2::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
