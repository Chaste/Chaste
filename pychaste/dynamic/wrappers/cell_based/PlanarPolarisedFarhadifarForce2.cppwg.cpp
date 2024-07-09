#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "PlanarPolarisedFarhadifarForce.hpp"

#include "PlanarPolarisedFarhadifarForce2.cppwg.hpp"

namespace py = pybind11;
typedef PlanarPolarisedFarhadifarForce<2 > PlanarPolarisedFarhadifarForce2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class PlanarPolarisedFarhadifarForce2_Overrides : public PlanarPolarisedFarhadifarForce2{
    public:
    using PlanarPolarisedFarhadifarForce2::PlanarPolarisedFarhadifarForce;
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            PlanarPolarisedFarhadifarForce2,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_PlanarPolarisedFarhadifarForce2_class(py::module &m){
py::class_<PlanarPolarisedFarhadifarForce2 , PlanarPolarisedFarhadifarForce2_Overrides , boost::shared_ptr<PlanarPolarisedFarhadifarForce2 >  , FarhadifarForce<2>  >(m, "PlanarPolarisedFarhadifarForce2")
        .def(py::init< >())
        .def(
            "GetPlanarPolarisedLineTensionMultiplier",
            (double(PlanarPolarisedFarhadifarForce2::*)()) &PlanarPolarisedFarhadifarForce2::GetPlanarPolarisedLineTensionMultiplier,
            " "  )
        .def(
            "SetPlanarPolarisedLineTensionMultiplier",
            (void(PlanarPolarisedFarhadifarForce2::*)(double)) &PlanarPolarisedFarhadifarForce2::SetPlanarPolarisedLineTensionMultiplier,
            " " , py::arg("planarPolarisedLineTensionMultiplier") )
        .def(
            "GetBoundaryLineTensionParameter",
            (double(PlanarPolarisedFarhadifarForce2::*)()) &PlanarPolarisedFarhadifarForce2::GetBoundaryLineTensionParameter,
            " "  )
        .def(
            "OutputForceParameters",
            (void(PlanarPolarisedFarhadifarForce2::*)(::out_stream &)) &PlanarPolarisedFarhadifarForce2::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
