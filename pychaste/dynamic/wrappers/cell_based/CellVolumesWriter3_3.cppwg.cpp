#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "AbstractCellPopulation.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellVolumesWriter.hpp"

#include "CellVolumesWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef CellVolumesWriter<3,3 > CellVolumesWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellVolumesWriter3_3_Overrides : public CellVolumesWriter3_3{
    public:
    using CellVolumesWriter3_3::CellVolumesWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellVolumesWriter3_3,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellVolumesWriter3_3,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellVolumesWriter3_3_class(py::module &m){
py::class_<CellVolumesWriter3_3 , CellVolumesWriter3_3_Overrides , boost::shared_ptr<CellVolumesWriter3_3 >  , AbstractCellWriter<3, 3>  >(m, "CellVolumesWriter3_3")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(CellVolumesWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellVolumesWriter3_3::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellVolumesWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellVolumesWriter3_3::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
