#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "AbstractCellPopulation.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "LegacyCellProliferativeTypesWriter.hpp"

#include "LegacyCellProliferativeTypesWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef LegacyCellProliferativeTypesWriter<2,2 > LegacyCellProliferativeTypesWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class LegacyCellProliferativeTypesWriter2_2_Overrides : public LegacyCellProliferativeTypesWriter2_2{
    public:
    using LegacyCellProliferativeTypesWriter2_2::LegacyCellProliferativeTypesWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            LegacyCellProliferativeTypesWriter2_2,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            LegacyCellProliferativeTypesWriter2_2,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_LegacyCellProliferativeTypesWriter2_2_class(py::module &m){
py::class_<LegacyCellProliferativeTypesWriter2_2 , LegacyCellProliferativeTypesWriter2_2_Overrides , boost::shared_ptr<LegacyCellProliferativeTypesWriter2_2 >  , AbstractCellWriter<2, 2>  >(m, "LegacyCellProliferativeTypesWriter2_2")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(LegacyCellProliferativeTypesWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &LegacyCellProliferativeTypesWriter2_2::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(LegacyCellProliferativeTypesWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &LegacyCellProliferativeTypesWriter2_2::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
