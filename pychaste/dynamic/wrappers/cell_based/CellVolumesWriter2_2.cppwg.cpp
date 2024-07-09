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

#include "CellVolumesWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef CellVolumesWriter<2,2 > CellVolumesWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellVolumesWriter2_2_Overrides : public CellVolumesWriter2_2{
    public:
    using CellVolumesWriter2_2::CellVolumesWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellVolumesWriter2_2,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellVolumesWriter2_2,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellVolumesWriter2_2_class(py::module &m){
py::class_<CellVolumesWriter2_2 , CellVolumesWriter2_2_Overrides , boost::shared_ptr<CellVolumesWriter2_2 >  , AbstractCellWriter<2, 2>  >(m, "CellVolumesWriter2_2")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(CellVolumesWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellVolumesWriter2_2::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellVolumesWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellVolumesWriter2_2::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
