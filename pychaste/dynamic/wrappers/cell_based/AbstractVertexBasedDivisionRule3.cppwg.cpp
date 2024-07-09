#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractVertexBasedDivisionRule.hpp"

#include "AbstractVertexBasedDivisionRule3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractVertexBasedDivisionRule<3 > AbstractVertexBasedDivisionRule3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;

class AbstractVertexBasedDivisionRule3_Overrides : public AbstractVertexBasedDivisionRule3{
    public:
    using AbstractVertexBasedDivisionRule3::AbstractVertexBasedDivisionRule;
    ::boost::numeric::ublas::c_vector<double, 3> CalculateCellDivisionVector(::CellPtr pParentCell, ::VertexBasedCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            AbstractVertexBasedDivisionRule3,
            CalculateCellDivisionVector,
                    pParentCell,
        rCellPopulation);
    }
    void OutputCellVertexBasedDivisionRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractVertexBasedDivisionRule3,
            OutputCellVertexBasedDivisionRuleParameters,
                    rParamsFile);
    }

};
void register_AbstractVertexBasedDivisionRule3_class(py::module &m){
py::class_<AbstractVertexBasedDivisionRule3 , AbstractVertexBasedDivisionRule3_Overrides , boost::shared_ptr<AbstractVertexBasedDivisionRule3 >   >(m, "AbstractVertexBasedDivisionRule3")
        .def(py::init< >())
        .def(
            "CalculateCellDivisionVector",
            (::boost::numeric::ublas::c_vector<double, 3>(AbstractVertexBasedDivisionRule3::*)(::CellPtr, ::VertexBasedCellPopulation<3> &)) &AbstractVertexBasedDivisionRule3::CalculateCellDivisionVector,
            " " , py::arg("pParentCell"), py::arg("rCellPopulation") )
        .def(
            "OutputCellVertexBasedDivisionRuleInfo",
            (void(AbstractVertexBasedDivisionRule3::*)(::out_stream &)) &AbstractVertexBasedDivisionRule3::OutputCellVertexBasedDivisionRuleInfo,
            " " , py::arg("rParamsFile") )
    ;
}
