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
#include "AbstractCellPopulationEventWriter.hpp"

#include "AbstractCellPopulationEventWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCellPopulationEventWriter<3,3 > AbstractCellPopulationEventWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractCellPopulationEventWriter3_3_Overrides : public AbstractCellPopulationEventWriter3_3{
    public:
    using AbstractCellPopulationEventWriter3_3::AbstractCellPopulationEventWriter;
    void WriteHeader(::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellPopulationEventWriter3_3,
            WriteHeader,
                    pCellPopulation);
    }
    void Visit(::MeshBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellPopulationEventWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::CaBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellPopulationEventWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::NodeBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellPopulationEventWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::PottsBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellPopulationEventWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::VertexBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellPopulationEventWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::ImmersedBoundaryCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellPopulationEventWriter3_3,
            Visit,
                    pCellPopulation);
    }

};
void register_AbstractCellPopulationEventWriter3_3_class(py::module &m){
py::class_<AbstractCellPopulationEventWriter3_3 , AbstractCellPopulationEventWriter3_3_Overrides , boost::shared_ptr<AbstractCellPopulationEventWriter3_3 >  , AbstractCellBasedWriter<3, 3>  >(m, "AbstractCellPopulationEventWriter3_3")
        .def(py::init<::std::string const & >(), py::arg("rFileName"))
        .def(
            "WriteHeader",
            (void(AbstractCellPopulationEventWriter3_3::*)(::AbstractCellPopulation<3> *)) &AbstractCellPopulationEventWriter3_3::WriteHeader,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(AbstractCellPopulationEventWriter3_3::*)(::MeshBasedCellPopulation<3> *)) &AbstractCellPopulationEventWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(AbstractCellPopulationEventWriter3_3::*)(::CaBasedCellPopulation<3> *)) &AbstractCellPopulationEventWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(AbstractCellPopulationEventWriter3_3::*)(::NodeBasedCellPopulation<3> *)) &AbstractCellPopulationEventWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(AbstractCellPopulationEventWriter3_3::*)(::PottsBasedCellPopulation<3> *)) &AbstractCellPopulationEventWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(AbstractCellPopulationEventWriter3_3::*)(::VertexBasedCellPopulation<3> *)) &AbstractCellPopulationEventWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(AbstractCellPopulationEventWriter3_3::*)(::ImmersedBoundaryCellPopulation<3> *)) &AbstractCellPopulationEventWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
    ;
}
