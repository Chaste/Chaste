#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellPopulationPyChasteActorGenerator.hpp"

#include "CellPopulationPyChasteActorGenerator3.cppwg.hpp"

namespace py = pybind11;
typedef CellPopulationPyChasteActorGenerator<3 > CellPopulationPyChasteActorGenerator3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellPopulationPyChasteActorGenerator3_Overrides : public CellPopulationPyChasteActorGenerator3{
    public:
    using CellPopulationPyChasteActorGenerator3::CellPopulationPyChasteActorGenerator;
    void AddActor(::vtkSmartPointer<vtkRenderer> pRenderer) override {
        PYBIND11_OVERRIDE(
            void,
            CellPopulationPyChasteActorGenerator3,
            AddActor,
                    pRenderer);
    }

};
void register_CellPopulationPyChasteActorGenerator3_class(py::module &m){
py::class_<CellPopulationPyChasteActorGenerator3 , CellPopulationPyChasteActorGenerator3_Overrides , boost::shared_ptr<CellPopulationPyChasteActorGenerator3 >  , AbstractPyChasteActorGenerator<3>  >(m, "CellPopulationPyChasteActorGenerator3")
        .def(py::init< >())
        .def(
            "AddActor",
            (void(CellPopulationPyChasteActorGenerator3::*)(::vtkSmartPointer<vtkRenderer>)) &CellPopulationPyChasteActorGenerator3::AddActor,
            " " , py::arg("pRenderer") )
        .def(
            "AddMeshBasedCellPopulationActor",
            (void(CellPopulationPyChasteActorGenerator3::*)(::vtkSmartPointer<vtkRenderer>)) &CellPopulationPyChasteActorGenerator3::AddMeshBasedCellPopulationActor,
            " " , py::arg("pRenderer") )
        .def(
            "AddVertexBasedCellPopulationActor",
            (void(CellPopulationPyChasteActorGenerator3::*)(::vtkSmartPointer<vtkRenderer>)) &CellPopulationPyChasteActorGenerator3::AddVertexBasedCellPopulationActor,
            " " , py::arg("pRenderer") )
        .def(
            "AddImmersedBoundaryCellPopulationActor",
            (void(CellPopulationPyChasteActorGenerator3::*)(::vtkSmartPointer<vtkRenderer>)) &CellPopulationPyChasteActorGenerator3::AddImmersedBoundaryCellPopulationActor,
            " " , py::arg("pRenderer") )
        .def(
            "AddCaBasedCellPopulationActor",
            (void(CellPopulationPyChasteActorGenerator3::*)(::vtkSmartPointer<vtkRenderer>)) &CellPopulationPyChasteActorGenerator3::AddCaBasedCellPopulationActor,
            " " , py::arg("pRenderer") )
        .def(
            "AddPottsBasedCellPopulationActor",
            (void(CellPopulationPyChasteActorGenerator3::*)(::vtkSmartPointer<vtkRenderer>)) &CellPopulationPyChasteActorGenerator3::AddPottsBasedCellPopulationActor,
            " " , py::arg("pRenderer") )
        .def(
            "SetCellPopulation",
            (void(CellPopulationPyChasteActorGenerator3::*)(::boost::shared_ptr<AbstractCellPopulation<3, 3>>)) &CellPopulationPyChasteActorGenerator3::SetCellPopulation,
            " " , py::arg("pCellPopulation") )
        .def(
            "SetShowVoronoiMeshEdges",
            (void(CellPopulationPyChasteActorGenerator3::*)(bool)) &CellPopulationPyChasteActorGenerator3::SetShowVoronoiMeshEdges,
            " " , py::arg("showEdges") )
        .def(
            "SetShowMutableMeshEdges",
            (void(CellPopulationPyChasteActorGenerator3::*)(bool)) &CellPopulationPyChasteActorGenerator3::SetShowMutableMeshEdges,
            " " , py::arg("showEdges") )
        .def(
            "SetShowPottsMeshEdges",
            (void(CellPopulationPyChasteActorGenerator3::*)(bool)) &CellPopulationPyChasteActorGenerator3::SetShowPottsMeshEdges,
            " " , py::arg("showEdges") )
        .def(
            "SetShowPottsMeshOutlines",
            (void(CellPopulationPyChasteActorGenerator3::*)(bool)) &CellPopulationPyChasteActorGenerator3::SetShowPottsMeshOutlines,
            " " , py::arg("showOutlines") )
        .def(
            "SetColorByCellType",
            (void(CellPopulationPyChasteActorGenerator3::*)(bool)) &CellPopulationPyChasteActorGenerator3::SetColorByCellType,
            " " , py::arg("colorByCellType") )
        .def(
            "SetColorByCellMutationState",
            (void(CellPopulationPyChasteActorGenerator3::*)(bool)) &CellPopulationPyChasteActorGenerator3::SetColorByCellMutationState,
            " " , py::arg("colorByCellMutationState") )
        .def(
            "SetColorByCellLabel",
            (void(CellPopulationPyChasteActorGenerator3::*)(bool)) &CellPopulationPyChasteActorGenerator3::SetColorByCellLabel,
            " " , py::arg("colorByCellLabel") )
        .def(
            "SetColorByUserDefined",
            (void(CellPopulationPyChasteActorGenerator3::*)(bool)) &CellPopulationPyChasteActorGenerator3::SetColorByUserDefined,
            " " , py::arg("colorByCellUserDefined") )
        .def(
            "SetColorByCellData",
            (void(CellPopulationPyChasteActorGenerator3::*)(bool)) &CellPopulationPyChasteActorGenerator3::SetColorByCellData,
            " " , py::arg("colorByCellData") )
        .def(
            "SetShowCellCentres",
            (void(CellPopulationPyChasteActorGenerator3::*)(bool)) &CellPopulationPyChasteActorGenerator3::SetShowCellCentres,
            " " , py::arg("showCentres") )
    ;
}
