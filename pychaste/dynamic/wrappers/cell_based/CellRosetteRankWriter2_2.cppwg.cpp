#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "AbstractCellPopulation.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellRosetteRankWriter.hpp"

#include "CellRosetteRankWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef CellRosetteRankWriter<2,2 > CellRosetteRankWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellRosetteRankWriter2_2_Overrides : public CellRosetteRankWriter2_2{
    public:
    using CellRosetteRankWriter2_2::CellRosetteRankWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellRosetteRankWriter2_2,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellRosetteRankWriter2_2,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellRosetteRankWriter2_2_class(py::module &m){
py::class_<CellRosetteRankWriter2_2 , CellRosetteRankWriter2_2_Overrides , boost::shared_ptr<CellRosetteRankWriter2_2 >  , AbstractCellWriter<2, 2>  >(m, "CellRosetteRankWriter2_2")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(CellRosetteRankWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellRosetteRankWriter2_2::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellRosetteRankWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellRosetteRankWriter2_2::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
