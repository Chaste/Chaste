#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "AbstractCellPopulation.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellCycleModelProteinConcentrationsWriter.hpp"

#include "CellCycleModelProteinConcentrationsWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef CellCycleModelProteinConcentrationsWriter<3,3 > CellCycleModelProteinConcentrationsWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellCycleModelProteinConcentrationsWriter3_3_Overrides : public CellCycleModelProteinConcentrationsWriter3_3{
    public:
    using CellCycleModelProteinConcentrationsWriter3_3::CellCycleModelProteinConcentrationsWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellCycleModelProteinConcentrationsWriter3_3,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellCycleModelProteinConcentrationsWriter3_3,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellCycleModelProteinConcentrationsWriter3_3_class(py::module &m){
py::class_<CellCycleModelProteinConcentrationsWriter3_3 , CellCycleModelProteinConcentrationsWriter3_3_Overrides , boost::shared_ptr<CellCycleModelProteinConcentrationsWriter3_3 >  , AbstractCellWriter<3, 3>  >(m, "CellCycleModelProteinConcentrationsWriter3_3")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(CellCycleModelProteinConcentrationsWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellCycleModelProteinConcentrationsWriter3_3::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellCycleModelProteinConcentrationsWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellCycleModelProteinConcentrationsWriter3_3::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
