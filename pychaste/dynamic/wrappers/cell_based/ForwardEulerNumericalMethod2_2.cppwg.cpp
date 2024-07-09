#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ForwardEulerNumericalMethod.hpp"

#include "ForwardEulerNumericalMethod2_2.cppwg.hpp"

namespace py = pybind11;
typedef ForwardEulerNumericalMethod<2,2 > ForwardEulerNumericalMethod2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ForwardEulerNumericalMethod2_2_Overrides : public ForwardEulerNumericalMethod2_2{
    public:
    using ForwardEulerNumericalMethod2_2::ForwardEulerNumericalMethod;
    void UpdateAllNodePositions(double dt) override {
        PYBIND11_OVERRIDE(
            void,
            ForwardEulerNumericalMethod2_2,
            UpdateAllNodePositions,
                    dt);
    }
    void OutputNumericalMethodParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ForwardEulerNumericalMethod2_2,
            OutputNumericalMethodParameters,
                    rParamsFile);
    }

};
void register_ForwardEulerNumericalMethod2_2_class(py::module &m){
py::class_<ForwardEulerNumericalMethod2_2 , ForwardEulerNumericalMethod2_2_Overrides , boost::shared_ptr<ForwardEulerNumericalMethod2_2 >  , AbstractNumericalMethod<2>  >(m, "ForwardEulerNumericalMethod2_2")
        .def(py::init< >())
        .def(
            "UpdateAllNodePositions",
            (void(ForwardEulerNumericalMethod2_2::*)(double)) &ForwardEulerNumericalMethod2_2::UpdateAllNodePositions,
            " " , py::arg("dt") )
        .def(
            "OutputNumericalMethodParameters",
            (void(ForwardEulerNumericalMethod2_2::*)(::out_stream &)) &ForwardEulerNumericalMethod2_2::OutputNumericalMethodParameters,
            " " , py::arg("rParamsFile") )
    ;
}
