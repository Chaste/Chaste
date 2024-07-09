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

#include "RadialCellDataDistributionWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef RadialCellDataDistributionWriter<2,2 > RadialCellDataDistributionWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class RadialCellDataDistributionWriter2_2_Overrides : public RadialCellDataDistributionWriter2_2{
    public:
    using RadialCellDataDistributionWriter2_2::RadialCellDataDistributionWriter;
    void Visit(::MeshBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            RadialCellDataDistributionWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::CaBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            RadialCellDataDistributionWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::NodeBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            RadialCellDataDistributionWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::PottsBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            RadialCellDataDistributionWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::VertexBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            RadialCellDataDistributionWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::ImmersedBoundaryCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            RadialCellDataDistributionWriter2_2,
            Visit,
                    pCellPopulation);
    }

};
void register_RadialCellDataDistributionWriter2_2_class(py::module &m){
py::class_<RadialCellDataDistributionWriter2_2 , RadialCellDataDistributionWriter2_2_Overrides , boost::shared_ptr<RadialCellDataDistributionWriter2_2 >  , AbstractCellPopulationWriter<2, 2>  >(m, "RadialCellDataDistributionWriter2_2")
        .def(py::init< >())
        .def(
            "VisitAnyPopulation",
            (void(RadialCellDataDistributionWriter2_2::*)(::AbstractCellPopulation<2> *)) &RadialCellDataDistributionWriter2_2::VisitAnyPopulation,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(RadialCellDataDistributionWriter2_2::*)(::MeshBasedCellPopulation<2> *)) &RadialCellDataDistributionWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(RadialCellDataDistributionWriter2_2::*)(::CaBasedCellPopulation<2> *)) &RadialCellDataDistributionWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(RadialCellDataDistributionWriter2_2::*)(::NodeBasedCellPopulation<2> *)) &RadialCellDataDistributionWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(RadialCellDataDistributionWriter2_2::*)(::PottsBasedCellPopulation<2> *)) &RadialCellDataDistributionWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(RadialCellDataDistributionWriter2_2::*)(::VertexBasedCellPopulation<2> *)) &RadialCellDataDistributionWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(RadialCellDataDistributionWriter2_2::*)(::ImmersedBoundaryCellPopulation<2> *)) &RadialCellDataDistributionWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "SetVariableName",
            (void(RadialCellDataDistributionWriter2_2::*)(::std::string)) &RadialCellDataDistributionWriter2_2::SetVariableName,
            " " , py::arg("variableName") )
        .def(
            "GetVariableName",
            (::std::string(RadialCellDataDistributionWriter2_2::*)() const ) &RadialCellDataDistributionWriter2_2::GetVariableName,
            " "  )
        .def(
            "SetNumRadialBins",
            (void(RadialCellDataDistributionWriter2_2::*)(unsigned int)) &RadialCellDataDistributionWriter2_2::SetNumRadialBins,
            " " , py::arg("numRadialBins") )
        .def(
            "GetNumRadialBins",
            (unsigned int(RadialCellDataDistributionWriter2_2::*)() const ) &RadialCellDataDistributionWriter2_2::GetNumRadialBins,
            " "  )
    ;
}
