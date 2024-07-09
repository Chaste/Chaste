#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ImmersedBoundaryElement.hpp"

#include "ImmersedBoundaryElement3_3.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryElement<3,3 > ImmersedBoundaryElement3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundaryElement3_3_Overrides : public ImmersedBoundaryElement3_3{
    public:
    using ImmersedBoundaryElement3_3::ImmersedBoundaryElement;
    bool IsElementOnBoundary() const  override {
        PYBIND11_OVERRIDE(
            bool,
            ImmersedBoundaryElement3_3,
            IsElementOnBoundary,
            );
    }

};
void register_ImmersedBoundaryElement3_3_class(py::module &m){
py::class_<ImmersedBoundaryElement3_3 , ImmersedBoundaryElement3_3_Overrides , boost::shared_ptr<ImmersedBoundaryElement3_3 >  , MutableElement<3, 3>  >(m, "ImmersedBoundaryElement3_3")
        .def(py::init<unsigned int, ::std::vector<Node<3> *> const & >(), py::arg("index"), py::arg("rNodes"))
        .def(
            "SetFluidSource",
            (void(ImmersedBoundaryElement3_3::*)(::std::shared_ptr<FluidSource<3>>)) &ImmersedBoundaryElement3_3::SetFluidSource,
            " " , py::arg("fluidSource") )
        .def(
            "GetFluidSource",
            (::std::shared_ptr<FluidSource<3>>(ImmersedBoundaryElement3_3::*)()) &ImmersedBoundaryElement3_3::GetFluidSource,
            " "  )
        .def(
            "rGetCornerNodes",
            (::std::vector<Node<3> *> &(ImmersedBoundaryElement3_3::*)()) &ImmersedBoundaryElement3_3::rGetCornerNodes,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetAverageNodeSpacing",
            (double(ImmersedBoundaryElement3_3::*)()) &ImmersedBoundaryElement3_3::GetAverageNodeSpacing,
            " "  )
        .def(
            "SetAverageNodeSpacing",
            (void(ImmersedBoundaryElement3_3::*)(double)) &ImmersedBoundaryElement3_3::SetAverageNodeSpacing,
            " " , py::arg("averageNodeSpacing") )
        .def(
            "IsElementOnBoundary",
            (bool(ImmersedBoundaryElement3_3::*)() const ) &ImmersedBoundaryElement3_3::IsElementOnBoundary,
            " "  )
        .def(
            "SetIsBoundaryElement",
            (void(ImmersedBoundaryElement3_3::*)(bool)) &ImmersedBoundaryElement3_3::SetIsBoundaryElement,
            " " , py::arg("isBoundaryElement") )
    ;
}
