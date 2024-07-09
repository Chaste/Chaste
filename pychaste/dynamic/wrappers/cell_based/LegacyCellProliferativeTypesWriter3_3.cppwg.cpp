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

#include "LegacyCellProliferativeTypesWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef LegacyCellProliferativeTypesWriter<3,3 > LegacyCellProliferativeTypesWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class LegacyCellProliferativeTypesWriter3_3_Overrides : public LegacyCellProliferativeTypesWriter3_3{
    public:
    using LegacyCellProliferativeTypesWriter3_3::LegacyCellProliferativeTypesWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            LegacyCellProliferativeTypesWriter3_3,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            LegacyCellProliferativeTypesWriter3_3,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_LegacyCellProliferativeTypesWriter3_3_class(py::module &m){
py::class_<LegacyCellProliferativeTypesWriter3_3 , LegacyCellProliferativeTypesWriter3_3_Overrides , boost::shared_ptr<LegacyCellProliferativeTypesWriter3_3 >  , AbstractCellWriter<3, 3>  >(m, "LegacyCellProliferativeTypesWriter3_3")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(LegacyCellProliferativeTypesWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &LegacyCellProliferativeTypesWriter3_3::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(LegacyCellProliferativeTypesWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &LegacyCellProliferativeTypesWriter3_3::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
