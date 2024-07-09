#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractPyChasteActorGenerator.hpp"

#include "AbstractPyChasteActorGenerator3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractPyChasteActorGenerator<3 > AbstractPyChasteActorGenerator3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractPyChasteActorGenerator3_Overrides : public AbstractPyChasteActorGenerator3{
    public:
    using AbstractPyChasteActorGenerator3::AbstractPyChasteActorGenerator;
    void AddActor(::vtkSmartPointer<vtkRenderer> pRenderer) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractPyChasteActorGenerator3,
            AddActor,
                    pRenderer);
    }

};
void register_AbstractPyChasteActorGenerator3_class(py::module &m){
py::class_<AbstractPyChasteActorGenerator3 , AbstractPyChasteActorGenerator3_Overrides , boost::shared_ptr<AbstractPyChasteActorGenerator3 >   >(m, "AbstractPyChasteActorGenerator3")
        .def(py::init< >())
        .def(
            "GetColorTransferFunction",
            (::vtkSmartPointer<vtkColorTransferFunction>(AbstractPyChasteActorGenerator3::*)()) &AbstractPyChasteActorGenerator3::GetColorTransferFunction,
            " "  )
        .def(
            "GetDiscreteColorTransferFunction",
            (::vtkSmartPointer<vtkColorTransferFunction>(AbstractPyChasteActorGenerator3::*)()) &AbstractPyChasteActorGenerator3::GetDiscreteColorTransferFunction,
            " "  )
        .def(
            "GetScaleBar",
            (::vtkSmartPointer<vtkScalarBarActor>(AbstractPyChasteActorGenerator3::*)()) &AbstractPyChasteActorGenerator3::GetScaleBar,
            " "  )
        .def(
            "AddActor",
            (void(AbstractPyChasteActorGenerator3::*)(::vtkSmartPointer<vtkRenderer>)) &AbstractPyChasteActorGenerator3::AddActor,
            " " , py::arg("pRenderer") )
        .def(
            "SetShowEdges",
            (void(AbstractPyChasteActorGenerator3::*)(bool)) &AbstractPyChasteActorGenerator3::SetShowEdges,
            " " , py::arg("show") )
        .def(
            "SetShowPoints",
            (void(AbstractPyChasteActorGenerator3::*)(bool)) &AbstractPyChasteActorGenerator3::SetShowPoints,
            " " , py::arg("show") )
        .def(
            "SetShowVolume",
            (void(AbstractPyChasteActorGenerator3::*)(bool)) &AbstractPyChasteActorGenerator3::SetShowVolume,
            " " , py::arg("show") )
        .def(
            "SetEdgeColor",
            (void(AbstractPyChasteActorGenerator3::*)(::boost::numeric::ublas::c_vector<double, 3> const &)) &AbstractPyChasteActorGenerator3::SetEdgeColor,
            " " , py::arg("rColor") )
        .def(
            "SetPointColor",
            (void(AbstractPyChasteActorGenerator3::*)(::boost::numeric::ublas::c_vector<double, 3> const &)) &AbstractPyChasteActorGenerator3::SetPointColor,
            " " , py::arg("rColor") )
        .def(
            "SetVolumeColor",
            (void(AbstractPyChasteActorGenerator3::*)(::boost::numeric::ublas::c_vector<double, 3> const &)) &AbstractPyChasteActorGenerator3::SetVolumeColor,
            " " , py::arg("rColor") )
        .def(
            "SetVolumeOpacity",
            (void(AbstractPyChasteActorGenerator3::*)(double)) &AbstractPyChasteActorGenerator3::SetVolumeOpacity,
            " " , py::arg("opacity") )
        .def(
            "SetPointSize",
            (void(AbstractPyChasteActorGenerator3::*)(double)) &AbstractPyChasteActorGenerator3::SetPointSize,
            " " , py::arg("size") )
        .def(
            "SetEdgeSize",
            (void(AbstractPyChasteActorGenerator3::*)(double)) &AbstractPyChasteActorGenerator3::SetEdgeSize,
            " " , py::arg("size") )
        .def(
            "SetDataLabel",
            (void(AbstractPyChasteActorGenerator3::*)(::std::string const &)) &AbstractPyChasteActorGenerator3::SetDataLabel,
            " " , py::arg("rLabel") )
        .def(
            "SetShowScaleBar",
            (void(AbstractPyChasteActorGenerator3::*)(double)) &AbstractPyChasteActorGenerator3::SetShowScaleBar,
            " " , py::arg("show") )
    ;
}
