#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractPyChasteActorGenerator.hpp"

#include "AbstractPyChasteActorGenerator2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractPyChasteActorGenerator<2 > AbstractPyChasteActorGenerator2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractPyChasteActorGenerator2_Overrides : public AbstractPyChasteActorGenerator2{
    public:
    using AbstractPyChasteActorGenerator2::AbstractPyChasteActorGenerator;
    void AddActor(::vtkSmartPointer<vtkRenderer> pRenderer) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractPyChasteActorGenerator2,
            AddActor,
                    pRenderer);
    }

};
void register_AbstractPyChasteActorGenerator2_class(py::module &m){
py::class_<AbstractPyChasteActorGenerator2 , AbstractPyChasteActorGenerator2_Overrides , boost::shared_ptr<AbstractPyChasteActorGenerator2 >   >(m, "AbstractPyChasteActorGenerator2")
        .def(py::init< >())
        .def(
            "GetColorTransferFunction",
            (::vtkSmartPointer<vtkColorTransferFunction>(AbstractPyChasteActorGenerator2::*)()) &AbstractPyChasteActorGenerator2::GetColorTransferFunction,
            " "  )
        .def(
            "GetDiscreteColorTransferFunction",
            (::vtkSmartPointer<vtkColorTransferFunction>(AbstractPyChasteActorGenerator2::*)()) &AbstractPyChasteActorGenerator2::GetDiscreteColorTransferFunction,
            " "  )
        .def(
            "GetScaleBar",
            (::vtkSmartPointer<vtkScalarBarActor>(AbstractPyChasteActorGenerator2::*)()) &AbstractPyChasteActorGenerator2::GetScaleBar,
            " "  )
        .def(
            "AddActor",
            (void(AbstractPyChasteActorGenerator2::*)(::vtkSmartPointer<vtkRenderer>)) &AbstractPyChasteActorGenerator2::AddActor,
            " " , py::arg("pRenderer") )
        .def(
            "SetShowEdges",
            (void(AbstractPyChasteActorGenerator2::*)(bool)) &AbstractPyChasteActorGenerator2::SetShowEdges,
            " " , py::arg("show") )
        .def(
            "SetShowPoints",
            (void(AbstractPyChasteActorGenerator2::*)(bool)) &AbstractPyChasteActorGenerator2::SetShowPoints,
            " " , py::arg("show") )
        .def(
            "SetShowVolume",
            (void(AbstractPyChasteActorGenerator2::*)(bool)) &AbstractPyChasteActorGenerator2::SetShowVolume,
            " " , py::arg("show") )
        .def(
            "SetEdgeColor",
            (void(AbstractPyChasteActorGenerator2::*)(::boost::numeric::ublas::c_vector<double, 3> const &)) &AbstractPyChasteActorGenerator2::SetEdgeColor,
            " " , py::arg("rColor") )
        .def(
            "SetPointColor",
            (void(AbstractPyChasteActorGenerator2::*)(::boost::numeric::ublas::c_vector<double, 3> const &)) &AbstractPyChasteActorGenerator2::SetPointColor,
            " " , py::arg("rColor") )
        .def(
            "SetVolumeColor",
            (void(AbstractPyChasteActorGenerator2::*)(::boost::numeric::ublas::c_vector<double, 3> const &)) &AbstractPyChasteActorGenerator2::SetVolumeColor,
            " " , py::arg("rColor") )
        .def(
            "SetVolumeOpacity",
            (void(AbstractPyChasteActorGenerator2::*)(double)) &AbstractPyChasteActorGenerator2::SetVolumeOpacity,
            " " , py::arg("opacity") )
        .def(
            "SetPointSize",
            (void(AbstractPyChasteActorGenerator2::*)(double)) &AbstractPyChasteActorGenerator2::SetPointSize,
            " " , py::arg("size") )
        .def(
            "SetEdgeSize",
            (void(AbstractPyChasteActorGenerator2::*)(double)) &AbstractPyChasteActorGenerator2::SetEdgeSize,
            " " , py::arg("size") )
        .def(
            "SetDataLabel",
            (void(AbstractPyChasteActorGenerator2::*)(::std::string const &)) &AbstractPyChasteActorGenerator2::SetDataLabel,
            " " , py::arg("rLabel") )
        .def(
            "SetShowScaleBar",
            (void(AbstractPyChasteActorGenerator2::*)(double)) &AbstractPyChasteActorGenerator2::SetShowScaleBar,
            " " , py::arg("show") )
    ;
}
