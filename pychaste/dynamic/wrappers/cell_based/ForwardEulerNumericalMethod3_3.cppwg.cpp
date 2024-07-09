#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ForwardEulerNumericalMethod.hpp"

#include "ForwardEulerNumericalMethod3_3.cppwg.hpp"

namespace py = pybind11;
typedef ForwardEulerNumericalMethod<3,3 > ForwardEulerNumericalMethod3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ForwardEulerNumericalMethod3_3_Overrides : public ForwardEulerNumericalMethod3_3{
    public:
    using ForwardEulerNumericalMethod3_3::ForwardEulerNumericalMethod;
    void UpdateAllNodePositions(double dt) override {
        PYBIND11_OVERRIDE(
            void,
            ForwardEulerNumericalMethod3_3,
            UpdateAllNodePositions,
                    dt);
    }
    void OutputNumericalMethodParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ForwardEulerNumericalMethod3_3,
            OutputNumericalMethodParameters,
                    rParamsFile);
    }

};
void register_ForwardEulerNumericalMethod3_3_class(py::module &m){
py::class_<ForwardEulerNumericalMethod3_3 , ForwardEulerNumericalMethod3_3_Overrides , boost::shared_ptr<ForwardEulerNumericalMethod3_3 >  , AbstractNumericalMethod<3>  >(m, "ForwardEulerNumericalMethod3_3")
        .def(py::init< >())
        .def(
            "UpdateAllNodePositions",
            (void(ForwardEulerNumericalMethod3_3::*)(double)) &ForwardEulerNumericalMethod3_3::UpdateAllNodePositions,
            " " , py::arg("dt") )
        .def(
            "OutputNumericalMethodParameters",
            (void(ForwardEulerNumericalMethod3_3::*)(::out_stream &)) &ForwardEulerNumericalMethod3_3::OutputNumericalMethodParameters,
            " " , py::arg("rParamsFile") )
    ;
}
