#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractChasteRegion.hpp"

#include "AbstractChasteRegion2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractChasteRegion<2 > AbstractChasteRegion2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractChasteRegion2_Overrides : public AbstractChasteRegion2{
    public:
    using AbstractChasteRegion2::AbstractChasteRegion;
    void Destroy() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractChasteRegion2,
            Destroy,
            );
    }
    bool DoesContain(::ChastePoint<2> const & rPointToCheck) const  override {
        PYBIND11_OVERRIDE_PURE(
            bool,
            AbstractChasteRegion2,
            DoesContain,
                    rPointToCheck);
    }

};
void register_AbstractChasteRegion2_class(py::module &m){
py::class_<AbstractChasteRegion2 , AbstractChasteRegion2_Overrides , boost::shared_ptr<AbstractChasteRegion2 >   >(m, "AbstractChasteRegion2")
        .def(py::init< >())
        .def(
            "Destroy",
            (void(AbstractChasteRegion2::*)()) &AbstractChasteRegion2::Destroy,
            " "  )
        .def(
            "DoesContain",
            (bool(AbstractChasteRegion2::*)(::ChastePoint<2> const &) const ) &AbstractChasteRegion2::DoesContain,
            " " , py::arg("rPointToCheck") )
    ;
}
