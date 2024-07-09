#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "PlanarPolarisedFarhadifarForce.hpp"

#include "PlanarPolarisedFarhadifarForce3.cppwg.hpp"

namespace py = pybind11;
typedef PlanarPolarisedFarhadifarForce<3 > PlanarPolarisedFarhadifarForce3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class PlanarPolarisedFarhadifarForce3_Overrides : public PlanarPolarisedFarhadifarForce3{
    public:
    using PlanarPolarisedFarhadifarForce3::PlanarPolarisedFarhadifarForce;
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            PlanarPolarisedFarhadifarForce3,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_PlanarPolarisedFarhadifarForce3_class(py::module &m){
py::class_<PlanarPolarisedFarhadifarForce3 , PlanarPolarisedFarhadifarForce3_Overrides , boost::shared_ptr<PlanarPolarisedFarhadifarForce3 >  , FarhadifarForce<3>  >(m, "PlanarPolarisedFarhadifarForce3")
        .def(py::init< >())
        .def(
            "GetPlanarPolarisedLineTensionMultiplier",
            (double(PlanarPolarisedFarhadifarForce3::*)()) &PlanarPolarisedFarhadifarForce3::GetPlanarPolarisedLineTensionMultiplier,
            " "  )
        .def(
            "SetPlanarPolarisedLineTensionMultiplier",
            (void(PlanarPolarisedFarhadifarForce3::*)(double)) &PlanarPolarisedFarhadifarForce3::SetPlanarPolarisedLineTensionMultiplier,
            " " , py::arg("planarPolarisedLineTensionMultiplier") )
        .def(
            "GetBoundaryLineTensionParameter",
            (double(PlanarPolarisedFarhadifarForce3::*)()) &PlanarPolarisedFarhadifarForce3::GetBoundaryLineTensionParameter,
            " "  )
        .def(
            "OutputForceParameters",
            (void(PlanarPolarisedFarhadifarForce3::*)(::out_stream &)) &PlanarPolarisedFarhadifarForce3::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
