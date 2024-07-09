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

#include "CellCycleModelProteinConcentrationsWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef CellCycleModelProteinConcentrationsWriter<2,2 > CellCycleModelProteinConcentrationsWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellCycleModelProteinConcentrationsWriter2_2_Overrides : public CellCycleModelProteinConcentrationsWriter2_2{
    public:
    using CellCycleModelProteinConcentrationsWriter2_2::CellCycleModelProteinConcentrationsWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellCycleModelProteinConcentrationsWriter2_2,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellCycleModelProteinConcentrationsWriter2_2,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellCycleModelProteinConcentrationsWriter2_2_class(py::module &m){
py::class_<CellCycleModelProteinConcentrationsWriter2_2 , CellCycleModelProteinConcentrationsWriter2_2_Overrides , boost::shared_ptr<CellCycleModelProteinConcentrationsWriter2_2 >  , AbstractCellWriter<2, 2>  >(m, "CellCycleModelProteinConcentrationsWriter2_2")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(CellCycleModelProteinConcentrationsWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellCycleModelProteinConcentrationsWriter2_2::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellCycleModelProteinConcentrationsWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellCycleModelProteinConcentrationsWriter2_2::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
