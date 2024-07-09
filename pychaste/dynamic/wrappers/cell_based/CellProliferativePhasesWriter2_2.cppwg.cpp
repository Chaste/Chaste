#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "AbstractCellPopulation.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellProliferativePhasesWriter.hpp"

#include "CellProliferativePhasesWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef CellProliferativePhasesWriter<2,2 > CellProliferativePhasesWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellProliferativePhasesWriter2_2_Overrides : public CellProliferativePhasesWriter2_2{
    public:
    using CellProliferativePhasesWriter2_2::CellProliferativePhasesWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellProliferativePhasesWriter2_2,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellProliferativePhasesWriter2_2,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellProliferativePhasesWriter2_2_class(py::module &m){
py::class_<CellProliferativePhasesWriter2_2 , CellProliferativePhasesWriter2_2_Overrides , boost::shared_ptr<CellProliferativePhasesWriter2_2 >  , AbstractCellWriter<2, 2>  >(m, "CellProliferativePhasesWriter2_2")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(CellProliferativePhasesWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellProliferativePhasesWriter2_2::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellProliferativePhasesWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellProliferativePhasesWriter2_2::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
