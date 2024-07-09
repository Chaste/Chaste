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

#include "CellProliferativeTypesWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef CellProliferativeTypesWriter<3,3 > CellProliferativeTypesWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellProliferativeTypesWriter3_3_Overrides : public CellProliferativeTypesWriter3_3{
    public:
    using CellProliferativeTypesWriter3_3::CellProliferativeTypesWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellProliferativeTypesWriter3_3,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellProliferativeTypesWriter3_3,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellProliferativeTypesWriter3_3_class(py::module &m){
py::class_<CellProliferativeTypesWriter3_3 , CellProliferativeTypesWriter3_3_Overrides , boost::shared_ptr<CellProliferativeTypesWriter3_3 >  , AbstractCellWriter<3, 3>  >(m, "CellProliferativeTypesWriter3_3")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(CellProliferativeTypesWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellProliferativeTypesWriter3_3::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellProliferativeTypesWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellProliferativeTypesWriter3_3::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
