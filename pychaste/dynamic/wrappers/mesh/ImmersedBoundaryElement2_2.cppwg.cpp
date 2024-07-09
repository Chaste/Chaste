#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ImmersedBoundaryElement.hpp"

#include "ImmersedBoundaryElement2_2.cppwg.hpp"

namespace py = pybind11;
typedef ImmersedBoundaryElement<2,2 > ImmersedBoundaryElement2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ImmersedBoundaryElement2_2_Overrides : public ImmersedBoundaryElement2_2{
    public:
    using ImmersedBoundaryElement2_2::ImmersedBoundaryElement;
    bool IsElementOnBoundary() const  override {
        PYBIND11_OVERRIDE(
            bool,
            ImmersedBoundaryElement2_2,
            IsElementOnBoundary,
            );
    }

};
void register_ImmersedBoundaryElement2_2_class(py::module &m){
py::class_<ImmersedBoundaryElement2_2 , ImmersedBoundaryElement2_2_Overrides , boost::shared_ptr<ImmersedBoundaryElement2_2 >  , MutableElement<2, 2>  >(m, "ImmersedBoundaryElement2_2")
        .def(py::init<unsigned int, ::std::vector<Node<2> *> const & >(), py::arg("index"), py::arg("rNodes"))
        .def(
            "SetFluidSource",
            (void(ImmersedBoundaryElement2_2::*)(::std::shared_ptr<FluidSource<2>>)) &ImmersedBoundaryElement2_2::SetFluidSource,
            " " , py::arg("fluidSource") )
        .def(
            "GetFluidSource",
            (::std::shared_ptr<FluidSource<2>>(ImmersedBoundaryElement2_2::*)()) &ImmersedBoundaryElement2_2::GetFluidSource,
            " "  )
        .def(
            "rGetCornerNodes",
            (::std::vector<Node<2> *> &(ImmersedBoundaryElement2_2::*)()) &ImmersedBoundaryElement2_2::rGetCornerNodes,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetAverageNodeSpacing",
            (double(ImmersedBoundaryElement2_2::*)()) &ImmersedBoundaryElement2_2::GetAverageNodeSpacing,
            " "  )
        .def(
            "SetAverageNodeSpacing",
            (void(ImmersedBoundaryElement2_2::*)(double)) &ImmersedBoundaryElement2_2::SetAverageNodeSpacing,
            " " , py::arg("averageNodeSpacing") )
        .def(
            "IsElementOnBoundary",
            (bool(ImmersedBoundaryElement2_2::*)() const ) &ImmersedBoundaryElement2_2::IsElementOnBoundary,
            " "  )
        .def(
            "SetIsBoundaryElement",
            (void(ImmersedBoundaryElement2_2::*)(bool)) &ImmersedBoundaryElement2_2::SetIsBoundaryElement,
            " " , py::arg("isBoundaryElement") )
    ;
}
