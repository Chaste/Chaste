#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "AbstractCellPopulation.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellProliferativeTypesWriter.hpp"

#include "CellProliferativeTypesWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef CellProliferativeTypesWriter<2,2 > CellProliferativeTypesWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellProliferativeTypesWriter2_2_Overrides : public CellProliferativeTypesWriter2_2{
    public:
    using CellProliferativeTypesWriter2_2::CellProliferativeTypesWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellProliferativeTypesWriter2_2,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellProliferativeTypesWriter2_2,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellProliferativeTypesWriter2_2_class(py::module &m){
py::class_<CellProliferativeTypesWriter2_2 , CellProliferativeTypesWriter2_2_Overrides , boost::shared_ptr<CellProliferativeTypesWriter2_2 >  , AbstractCellWriter<2, 2>  >(m, "CellProliferativeTypesWriter2_2")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(CellProliferativeTypesWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellProliferativeTypesWriter2_2::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellProliferativeTypesWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellProliferativeTypesWriter2_2::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
