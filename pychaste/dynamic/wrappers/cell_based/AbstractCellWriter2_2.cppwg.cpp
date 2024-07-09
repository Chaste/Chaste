#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "AbstractCellPopulation.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractCellWriter.hpp"

#include "AbstractCellWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCellWriter<2,2 > AbstractCellWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;

class AbstractCellWriter2_2_Overrides : public AbstractCellWriter2_2{
    public:
    using AbstractCellWriter2_2::AbstractCellWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            AbstractCellWriter2_2,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    ::boost::numeric::ublas::c_vector<double, 2> GetVectorCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            AbstractCellWriter2_2,
            GetVectorCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellWriter2_2,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_AbstractCellWriter2_2_class(py::module &m){
py::class_<AbstractCellWriter2_2 , AbstractCellWriter2_2_Overrides , boost::shared_ptr<AbstractCellWriter2_2 >  , AbstractCellBasedWriter<2, 2>  >(m, "AbstractCellWriter2_2")
        .def(py::init<::std::string const & >(), py::arg("rFileName"))
        .def(
            "GetCellDataForVtkOutput",
            (double(AbstractCellWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &AbstractCellWriter2_2::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "GetVectorCellDataForVtkOutput",
            (::boost::numeric::ublas::c_vector<double, 2>(AbstractCellWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &AbstractCellWriter2_2::GetVectorCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(AbstractCellWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &AbstractCellWriter2_2::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "GetOutputScalarData",
            (bool(AbstractCellWriter2_2::*)()) &AbstractCellWriter2_2::GetOutputScalarData,
            " "  )
        .def(
            "GetOutputVectorData",
            (bool(AbstractCellWriter2_2::*)()) &AbstractCellWriter2_2::GetOutputVectorData,
            " "  )
        .def(
            "SetVtkCellDataName",
            (void(AbstractCellWriter2_2::*)(::std::string)) &AbstractCellWriter2_2::SetVtkCellDataName,
            " " , py::arg("vtkCellDataName") )
        .def(
            "SetVtkVectorCellDataName",
            (void(AbstractCellWriter2_2::*)(::std::string)) &AbstractCellWriter2_2::SetVtkVectorCellDataName,
            " " , py::arg("vtkCellDataName") )
        .def(
            "GetVtkCellDataName",
            (::std::string(AbstractCellWriter2_2::*)()) &AbstractCellWriter2_2::GetVtkCellDataName,
            " "  )
        .def(
            "GetVtkVectorCellDataName",
            (::std::string(AbstractCellWriter2_2::*)()) &AbstractCellWriter2_2::GetVtkVectorCellDataName,
            " "  )
    ;
}
