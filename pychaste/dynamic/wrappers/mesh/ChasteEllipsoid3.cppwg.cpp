#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ChasteEllipsoid.hpp"

#include "ChasteEllipsoid3.cppwg.hpp"

namespace py = pybind11;
typedef ChasteEllipsoid<3 > ChasteEllipsoid3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ChasteEllipsoid3_Overrides : public ChasteEllipsoid3{
    public:
    using ChasteEllipsoid3::ChasteEllipsoid;
    bool DoesContain(::ChastePoint<3> const & rPointToCheck) const  override {
        PYBIND11_OVERRIDE(
            bool,
            ChasteEllipsoid3,
            DoesContain,
                    rPointToCheck);
    }

};
void register_ChasteEllipsoid3_class(py::module &m){
py::class_<ChasteEllipsoid3 , ChasteEllipsoid3_Overrides , boost::shared_ptr<ChasteEllipsoid3 >  , AbstractChasteRegion<3>  >(m, "ChasteEllipsoid3")
        .def(py::init<::ChastePoint<3> &, ::ChastePoint<3> & >(), py::arg("rCentre"), py::arg("rRadii"))
        .def(
            "DoesContain",
            (bool(ChasteEllipsoid3::*)(::ChastePoint<3> const &) const ) &ChasteEllipsoid3::DoesContain,
            " " , py::arg("rPointToCheck") )
        .def(
            "rGetCentre",
            (::ChastePoint<3> const &(ChasteEllipsoid3::*)() const ) &ChasteEllipsoid3::rGetCentre,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "rGetRadii",
            (::ChastePoint<3> const &(ChasteEllipsoid3::*)() const ) &ChasteEllipsoid3::rGetRadii,
            " "  , py::return_value_policy::reference_internal)
    ;
}
