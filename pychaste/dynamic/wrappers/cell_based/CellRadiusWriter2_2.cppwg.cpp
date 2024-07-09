#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "AbstractCellPopulation.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellRadiusWriter.hpp"

#include "CellRadiusWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef CellRadiusWriter<2,2 > CellRadiusWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellRadiusWriter2_2_Overrides : public CellRadiusWriter2_2{
    public:
    using CellRadiusWriter2_2::CellRadiusWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellRadiusWriter2_2,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellRadiusWriter2_2,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellRadiusWriter2_2_class(py::module &m){
py::class_<CellRadiusWriter2_2 , CellRadiusWriter2_2_Overrides , boost::shared_ptr<CellRadiusWriter2_2 >  , AbstractCellWriter<2, 2>  >(m, "CellRadiusWriter2_2")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(CellRadiusWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellRadiusWriter2_2::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellRadiusWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellRadiusWriter2_2::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
