#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "AbstractCellPopulation.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellDataItemWriter.hpp"

#include "CellDataItemWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef CellDataItemWriter<3,3 > CellDataItemWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellDataItemWriter3_3_Overrides : public CellDataItemWriter3_3{
    public:
    using CellDataItemWriter3_3::CellDataItemWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellDataItemWriter3_3,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellDataItemWriter3_3,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellDataItemWriter3_3_class(py::module &m){
py::class_<CellDataItemWriter3_3 , CellDataItemWriter3_3_Overrides , boost::shared_ptr<CellDataItemWriter3_3 >  , AbstractCellWriter<3, 3>  >(m, "CellDataItemWriter3_3")
        .def(py::init<::std::string >(), py::arg("cellDataVariableName") = "")
        .def(
            "GetCellDataForVtkOutput",
            (double(CellDataItemWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellDataItemWriter3_3::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellDataItemWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellDataItemWriter3_3::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "GetCellDataVariableName",
            (::std::string(CellDataItemWriter3_3::*)() const ) &CellDataItemWriter3_3::GetCellDataVariableName,
            " "  )
    ;
}
