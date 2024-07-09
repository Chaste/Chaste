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

#include "AbstractCellWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCellWriter<3,3 > AbstractCellWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;

class AbstractCellWriter3_3_Overrides : public AbstractCellWriter3_3{
    public:
    using AbstractCellWriter3_3::AbstractCellWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            AbstractCellWriter3_3,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    ::boost::numeric::ublas::c_vector<double, 3> GetVectorCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            AbstractCellWriter3_3,
            GetVectorCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellWriter3_3,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_AbstractCellWriter3_3_class(py::module &m){
py::class_<AbstractCellWriter3_3 , AbstractCellWriter3_3_Overrides , boost::shared_ptr<AbstractCellWriter3_3 >  , AbstractCellBasedWriter<3, 3>  >(m, "AbstractCellWriter3_3")
        .def(py::init<::std::string const & >(), py::arg("rFileName"))
        .def(
            "GetCellDataForVtkOutput",
            (double(AbstractCellWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &AbstractCellWriter3_3::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "GetVectorCellDataForVtkOutput",
            (::boost::numeric::ublas::c_vector<double, 3>(AbstractCellWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &AbstractCellWriter3_3::GetVectorCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(AbstractCellWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &AbstractCellWriter3_3::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "GetOutputScalarData",
            (bool(AbstractCellWriter3_3::*)()) &AbstractCellWriter3_3::GetOutputScalarData,
            " "  )
        .def(
            "GetOutputVectorData",
            (bool(AbstractCellWriter3_3::*)()) &AbstractCellWriter3_3::GetOutputVectorData,
            " "  )
        .def(
            "SetVtkCellDataName",
            (void(AbstractCellWriter3_3::*)(::std::string)) &AbstractCellWriter3_3::SetVtkCellDataName,
            " " , py::arg("vtkCellDataName") )
        .def(
            "SetVtkVectorCellDataName",
            (void(AbstractCellWriter3_3::*)(::std::string)) &AbstractCellWriter3_3::SetVtkVectorCellDataName,
            " " , py::arg("vtkCellDataName") )
        .def(
            "GetVtkCellDataName",
            (::std::string(AbstractCellWriter3_3::*)()) &AbstractCellWriter3_3::GetVtkCellDataName,
            " "  )
        .def(
            "GetVtkVectorCellDataName",
            (::std::string(AbstractCellWriter3_3::*)()) &AbstractCellWriter3_3::GetVtkVectorCellDataName,
            " "  )
    ;
}
