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
#include "AbstractCellPopulationWriter.hpp"

#include "AbstractCellPopulationWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCellPopulationWriter<2,2 > AbstractCellPopulationWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractCellPopulationWriter2_2_Overrides : public AbstractCellPopulationWriter2_2{
    public:
    using AbstractCellPopulationWriter2_2::AbstractCellPopulationWriter;
    void WriteHeader(::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCellPopulationWriter2_2,
            WriteHeader,
                    pCellPopulation);
    }
    void Visit(::MeshBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellPopulationWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::CaBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellPopulationWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::NodeBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellPopulationWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::PottsBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellPopulationWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::VertexBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellPopulationWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::ImmersedBoundaryCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCellPopulationWriter2_2,
            Visit,
                    pCellPopulation);
    }

};
void register_AbstractCellPopulationWriter2_2_class(py::module &m){
py::class_<AbstractCellPopulationWriter2_2 , AbstractCellPopulationWriter2_2_Overrides , boost::shared_ptr<AbstractCellPopulationWriter2_2 >  , AbstractCellBasedWriter<2, 2>  >(m, "AbstractCellPopulationWriter2_2")
        .def(py::init<::std::string const & >(), py::arg("rFileName"))
        .def(
            "WriteHeader",
            (void(AbstractCellPopulationWriter2_2::*)(::AbstractCellPopulation<2> *)) &AbstractCellPopulationWriter2_2::WriteHeader,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(AbstractCellPopulationWriter2_2::*)(::MeshBasedCellPopulation<2> *)) &AbstractCellPopulationWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(AbstractCellPopulationWriter2_2::*)(::CaBasedCellPopulation<2> *)) &AbstractCellPopulationWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(AbstractCellPopulationWriter2_2::*)(::NodeBasedCellPopulation<2> *)) &AbstractCellPopulationWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(AbstractCellPopulationWriter2_2::*)(::PottsBasedCellPopulation<2> *)) &AbstractCellPopulationWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(AbstractCellPopulationWriter2_2::*)(::VertexBasedCellPopulation<2> *)) &AbstractCellPopulationWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(AbstractCellPopulationWriter2_2::*)(::ImmersedBoundaryCellPopulation<2> *)) &AbstractCellPopulationWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
    ;
}
