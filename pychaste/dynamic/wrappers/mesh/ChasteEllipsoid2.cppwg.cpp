#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ChasteEllipsoid.hpp"

#include "ChasteEllipsoid2.cppwg.hpp"

namespace py = pybind11;
typedef ChasteEllipsoid<2 > ChasteEllipsoid2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ChasteEllipsoid2_Overrides : public ChasteEllipsoid2{
    public:
    using ChasteEllipsoid2::ChasteEllipsoid;
    bool DoesContain(::ChastePoint<2> const & rPointToCheck) const  override {
        PYBIND11_OVERRIDE(
            bool,
            ChasteEllipsoid2,
            DoesContain,
                    rPointToCheck);
    }

};
void register_ChasteEllipsoid2_class(py::module &m){
py::class_<ChasteEllipsoid2 , ChasteEllipsoid2_Overrides , boost::shared_ptr<ChasteEllipsoid2 >  , AbstractChasteRegion<2>  >(m, "ChasteEllipsoid2")
        .def(py::init<::ChastePoint<2> &, ::ChastePoint<2> & >(), py::arg("rCentre"), py::arg("rRadii"))
        .def(
            "DoesContain",
            (bool(ChasteEllipsoid2::*)(::ChastePoint<2> const &) const ) &ChasteEllipsoid2::DoesContain,
            " " , py::arg("rPointToCheck") )
        .def(
            "rGetCentre",
            (::ChastePoint<2> const &(ChasteEllipsoid2::*)() const ) &ChasteEllipsoid2::rGetCentre,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "rGetRadii",
            (::ChastePoint<2> const &(ChasteEllipsoid2::*)() const ) &ChasteEllipsoid2::rGetRadii,
            " "  , py::return_value_policy::reference_internal)
    ;
}
