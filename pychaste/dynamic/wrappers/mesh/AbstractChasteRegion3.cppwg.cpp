#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractChasteRegion.hpp"

#include "AbstractChasteRegion3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractChasteRegion<3 > AbstractChasteRegion3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractChasteRegion3_Overrides : public AbstractChasteRegion3{
    public:
    using AbstractChasteRegion3::AbstractChasteRegion;
    void Destroy() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractChasteRegion3,
            Destroy,
            );
    }
    bool DoesContain(::ChastePoint<3> const & rPointToCheck) const  override {
        PYBIND11_OVERRIDE_PURE(
            bool,
            AbstractChasteRegion3,
            DoesContain,
                    rPointToCheck);
    }

};
void register_AbstractChasteRegion3_class(py::module &m){
py::class_<AbstractChasteRegion3 , AbstractChasteRegion3_Overrides , boost::shared_ptr<AbstractChasteRegion3 >   >(m, "AbstractChasteRegion3")
        .def(py::init< >())
        .def(
            "Destroy",
            (void(AbstractChasteRegion3::*)()) &AbstractChasteRegion3::Destroy,
            " "  )
        .def(
            "DoesContain",
            (bool(AbstractChasteRegion3::*)(::ChastePoint<3> const &) const ) &AbstractChasteRegion3::DoesContain,
            " " , py::arg("rPointToCheck") )
    ;
}
