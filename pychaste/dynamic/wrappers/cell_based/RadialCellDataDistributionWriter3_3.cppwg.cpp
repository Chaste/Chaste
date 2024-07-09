#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "RadialCellDataDistributionWriter.hpp"

#include "RadialCellDataDistributionWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef RadialCellDataDistributionWriter<3,3 > RadialCellDataDistributionWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class RadialCellDataDistributionWriter3_3_Overrides : public RadialCellDataDistributionWriter3_3{
    public:
    using RadialCellDataDistributionWriter3_3::RadialCellDataDistributionWriter;
    void Visit(::MeshBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            RadialCellDataDistributionWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::CaBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            RadialCellDataDistributionWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::NodeBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            RadialCellDataDistributionWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::PottsBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            RadialCellDataDistributionWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::VertexBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            RadialCellDataDistributionWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::ImmersedBoundaryCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            RadialCellDataDistributionWriter3_3,
            Visit,
                    pCellPopulation);
    }

};
void register_RadialCellDataDistributionWriter3_3_class(py::module &m){
py::class_<RadialCellDataDistributionWriter3_3 , RadialCellDataDistributionWriter3_3_Overrides , boost::shared_ptr<RadialCellDataDistributionWriter3_3 >  , AbstractCellPopulationWriter<3, 3>  >(m, "RadialCellDataDistributionWriter3_3")
        .def(py::init< >())
        .def(
            "VisitAnyPopulation",
            (void(RadialCellDataDistributionWriter3_3::*)(::AbstractCellPopulation<3> *)) &RadialCellDataDistributionWriter3_3::VisitAnyPopulation,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(RadialCellDataDistributionWriter3_3::*)(::MeshBasedCellPopulation<3> *)) &RadialCellDataDistributionWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(RadialCellDataDistributionWriter3_3::*)(::CaBasedCellPopulation<3> *)) &RadialCellDataDistributionWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(RadialCellDataDistributionWriter3_3::*)(::NodeBasedCellPopulation<3> *)) &RadialCellDataDistributionWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(RadialCellDataDistributionWriter3_3::*)(::PottsBasedCellPopulation<3> *)) &RadialCellDataDistributionWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(RadialCellDataDistributionWriter3_3::*)(::VertexBasedCellPopulation<3> *)) &RadialCellDataDistributionWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(RadialCellDataDistributionWriter3_3::*)(::ImmersedBoundaryCellPopulation<3> *)) &RadialCellDataDistributionWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "SetVariableName",
            (void(RadialCellDataDistributionWriter3_3::*)(::std::string)) &RadialCellDataDistributionWriter3_3::SetVariableName,
            " " , py::arg("variableName") )
        .def(
            "GetVariableName",
            (::std::string(RadialCellDataDistributionWriter3_3::*)() const ) &RadialCellDataDistributionWriter3_3::GetVariableName,
            " "  )
        .def(
            "SetNumRadialBins",
            (void(RadialCellDataDistributionWriter3_3::*)(unsigned int)) &RadialCellDataDistributionWriter3_3::SetNumRadialBins,
            " " , py::arg("numRadialBins") )
        .def(
            "GetNumRadialBins",
            (unsigned int(RadialCellDataDistributionWriter3_3::*)() const ) &RadialCellDataDistributionWriter3_3::GetNumRadialBins,
            " "  )
    ;
}
