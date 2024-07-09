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

#include "CellRosetteRankWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef CellRosetteRankWriter<3,3 > CellRosetteRankWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellRosetteRankWriter3_3_Overrides : public CellRosetteRankWriter3_3{
    public:
    using CellRosetteRankWriter3_3::CellRosetteRankWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellRosetteRankWriter3_3,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellRosetteRankWriter3_3,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellRosetteRankWriter3_3_class(py::module &m){
py::class_<CellRosetteRankWriter3_3 , CellRosetteRankWriter3_3_Overrides , boost::shared_ptr<CellRosetteRankWriter3_3 >  , AbstractCellWriter<3, 3>  >(m, "CellRosetteRankWriter3_3")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(CellRosetteRankWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellRosetteRankWriter3_3::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellRosetteRankWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellRosetteRankWriter3_3::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
