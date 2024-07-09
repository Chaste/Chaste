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

#include "CellLocationIndexWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef CellLocationIndexWriter<3,3 > CellLocationIndexWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellLocationIndexWriter3_3_Overrides : public CellLocationIndexWriter3_3{
    public:
    using CellLocationIndexWriter3_3::CellLocationIndexWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellLocationIndexWriter3_3,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellLocationIndexWriter3_3,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellLocationIndexWriter3_3_class(py::module &m){
py::class_<CellLocationIndexWriter3_3 , CellLocationIndexWriter3_3_Overrides , boost::shared_ptr<CellLocationIndexWriter3_3 >  , AbstractCellWriter<3, 3>  >(m, "CellLocationIndexWriter3_3")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(CellLocationIndexWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellLocationIndexWriter3_3::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellLocationIndexWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellLocationIndexWriter3_3::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
