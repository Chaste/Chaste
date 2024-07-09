#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractBoundaryCondition.hpp"

#include "AbstractBoundaryCondition2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractBoundaryCondition<2 > AbstractBoundaryCondition2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractBoundaryCondition2_Overrides : public AbstractBoundaryCondition2{
    public:
    using AbstractBoundaryCondition2::AbstractBoundaryCondition;
    double GetValue(::ChastePoint<2> const & rX) const  override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractBoundaryCondition2,
            GetValue,
                    rX);
    }

};
void register_AbstractBoundaryCondition2_class(py::module &m){
py::class_<AbstractBoundaryCondition2 , AbstractBoundaryCondition2_Overrides , boost::shared_ptr<AbstractBoundaryCondition2 >   >(m, "AbstractBoundaryCondition2")
        .def(py::init< >())
        .def(
            "GetValue",
            (double(AbstractBoundaryCondition2::*)(::ChastePoint<2> const &) const ) &AbstractBoundaryCondition2::GetValue,
            " " , py::arg("rX") )
    ;
}
