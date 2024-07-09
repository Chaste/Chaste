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

#include "CellDeltaNotchWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef CellDeltaNotchWriter<3,3 > CellDeltaNotchWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellDeltaNotchWriter3_3_Overrides : public CellDeltaNotchWriter3_3{
    public:
    using CellDeltaNotchWriter3_3::CellDeltaNotchWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellDeltaNotchWriter3_3,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellDeltaNotchWriter3_3,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellDeltaNotchWriter3_3_class(py::module &m){
py::class_<CellDeltaNotchWriter3_3 , CellDeltaNotchWriter3_3_Overrides , boost::shared_ptr<CellDeltaNotchWriter3_3 >  , AbstractCellWriter<3, 3>  >(m, "CellDeltaNotchWriter3_3")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(CellDeltaNotchWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellDeltaNotchWriter3_3::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellDeltaNotchWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellDeltaNotchWriter3_3::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
