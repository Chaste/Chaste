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

#include "CellAgesWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef CellAgesWriter<3,3 > CellAgesWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellAgesWriter3_3_Overrides : public CellAgesWriter3_3{
    public:
    using CellAgesWriter3_3::CellAgesWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellAgesWriter3_3,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellAgesWriter3_3,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellAgesWriter3_3_class(py::module &m){
py::class_<CellAgesWriter3_3 , CellAgesWriter3_3_Overrides , boost::shared_ptr<CellAgesWriter3_3 >  , AbstractCellWriter<3, 3>  >(m, "CellAgesWriter3_3")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(CellAgesWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellAgesWriter3_3::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellAgesWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellAgesWriter3_3::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
