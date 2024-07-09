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

#include "CellDataItemWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef CellDataItemWriter<2,2 > CellDataItemWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellDataItemWriter2_2_Overrides : public CellDataItemWriter2_2{
    public:
    using CellDataItemWriter2_2::CellDataItemWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellDataItemWriter2_2,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellDataItemWriter2_2,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellDataItemWriter2_2_class(py::module &m){
py::class_<CellDataItemWriter2_2 , CellDataItemWriter2_2_Overrides , boost::shared_ptr<CellDataItemWriter2_2 >  , AbstractCellWriter<2, 2>  >(m, "CellDataItemWriter2_2")
        .def(py::init<::std::string >(), py::arg("cellDataVariableName") = "")
        .def(
            "GetCellDataForVtkOutput",
            (double(CellDataItemWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellDataItemWriter2_2::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellDataItemWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellDataItemWriter2_2::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "GetCellDataVariableName",
            (::std::string(CellDataItemWriter2_2::*)() const ) &CellDataItemWriter2_2::GetCellDataVariableName,
            " "  )
    ;
}
