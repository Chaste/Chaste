#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "AbstractCellPopulation.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellRadiusWriter.hpp"

#include "CellRadiusWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef CellRadiusWriter<3,3 > CellRadiusWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellRadiusWriter3_3_Overrides : public CellRadiusWriter3_3{
    public:
    using CellRadiusWriter3_3::CellRadiusWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellRadiusWriter3_3,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellRadiusWriter3_3,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellRadiusWriter3_3_class(py::module &m){
py::class_<CellRadiusWriter3_3 , CellRadiusWriter3_3_Overrides , boost::shared_ptr<CellRadiusWriter3_3 >  , AbstractCellWriter<3, 3>  >(m, "CellRadiusWriter3_3")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(CellRadiusWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellRadiusWriter3_3::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellRadiusWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellRadiusWriter3_3::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
