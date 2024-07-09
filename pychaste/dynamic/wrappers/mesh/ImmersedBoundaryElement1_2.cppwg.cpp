#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ImmersedBoundaryElement.hpp"

#include "ImmersedBoundaryElement1_2.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryElement<1,2 > ImmersedBoundaryElement1_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundaryElement1_2_Overrides : public ImmersedBoundaryElement1_2{
    public:
    using ImmersedBoundaryElement1_2::ImmersedBoundaryElement;
    bool IsElementOnBoundary() const  override {
        PYBIND11_OVERRIDE(
            bool,
            ImmersedBoundaryElement1_2,
            IsElementOnBoundary,
            );
    }

};
void register_ImmersedBoundaryElement1_2_class(py::module &m){
py::class_<ImmersedBoundaryElement1_2 , ImmersedBoundaryElement1_2_Overrides , boost::shared_ptr<ImmersedBoundaryElement1_2 >  , MutableElement<1, 2>  >(m, "ImmersedBoundaryElement1_2")
        .def(py::init<unsigned int, ::std::vector<Node<2> *> const & >(), py::arg("index"), py::arg("rNodes"))
        .def(
            "SetFluidSource",
            (void(ImmersedBoundaryElement1_2::*)(::std::shared_ptr<FluidSource<2>>)) &ImmersedBoundaryElement1_2::SetFluidSource,
            " " , py::arg("fluidSource") )
        .def(
            "GetFluidSource",
            (::std::shared_ptr<FluidSource<2>>(ImmersedBoundaryElement1_2::*)()) &ImmersedBoundaryElement1_2::GetFluidSource,
            " "  )
        .def(
            "rGetCornerNodes",
            (::std::vector<Node<2> *> &(ImmersedBoundaryElement1_2::*)()) &ImmersedBoundaryElement1_2::rGetCornerNodes,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetAverageNodeSpacing",
            (double(ImmersedBoundaryElement1_2::*)()) &ImmersedBoundaryElement1_2::GetAverageNodeSpacing,
            " "  )
        .def(
            "SetAverageNodeSpacing",
            (void(ImmersedBoundaryElement1_2::*)(double)) &ImmersedBoundaryElement1_2::SetAverageNodeSpacing,
            " " , py::arg("averageNodeSpacing") )
        .def(
            "IsElementOnBoundary",
            (bool(ImmersedBoundaryElement1_2::*)() const ) &ImmersedBoundaryElement1_2::IsElementOnBoundary,
            " "  )
        .def(
            "SetIsBoundaryElement",
            (void(ImmersedBoundaryElement1_2::*)(bool)) &ImmersedBoundaryElement1_2::SetIsBoundaryElement,
            " " , py::arg("isBoundaryElement") )
    ;
}
