#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "AbstractCellPopulation.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellDeltaNotchWriter.hpp"

#include "CellDeltaNotchWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef CellDeltaNotchWriter<2,2 > CellDeltaNotchWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellDeltaNotchWriter2_2_Overrides : public CellDeltaNotchWriter2_2{
    public:
    using CellDeltaNotchWriter2_2::CellDeltaNotchWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellDeltaNotchWriter2_2,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellDeltaNotchWriter2_2,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellDeltaNotchWriter2_2_class(py::module &m){
py::class_<CellDeltaNotchWriter2_2 , CellDeltaNotchWriter2_2_Overrides , boost::shared_ptr<CellDeltaNotchWriter2_2 >  , AbstractCellWriter<2, 2>  >(m, "CellDeltaNotchWriter2_2")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(CellDeltaNotchWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellDeltaNotchWriter2_2::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellDeltaNotchWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellDeltaNotchWriter2_2::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
