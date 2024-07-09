#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractUpdateRule.hpp"

#include "AbstractUpdateRule2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractUpdateRule<2 > AbstractUpdateRule2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractUpdateRule2_Overrides : public AbstractUpdateRule2{
    public:
    using AbstractUpdateRule2::AbstractUpdateRule;
    void OutputUpdateRuleParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractUpdateRule2,
            OutputUpdateRuleParameters,
                    rParamsFile);
    }

};
void register_AbstractUpdateRule2_class(py::module &m){
py::class_<AbstractUpdateRule2 , AbstractUpdateRule2_Overrides , boost::shared_ptr<AbstractUpdateRule2 >   >(m, "AbstractUpdateRule2")
        .def(py::init< >())
        .def(
            "OutputUpdateRuleInfo",
            (void(AbstractUpdateRule2::*)(::out_stream &)) &AbstractUpdateRule2::OutputUpdateRuleInfo,
            " " , py::arg("rParamsFile") )
        .def(
            "OutputUpdateRuleParameters",
            (void(AbstractUpdateRule2::*)(::out_stream &)) &AbstractUpdateRule2::OutputUpdateRuleParameters,
            " " , py::arg("rParamsFile") )
    ;
}
