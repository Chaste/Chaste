#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellBasedParabolicPdeSolver.hpp"

#include "CellBasedParabolicPdeSolver3.cppwg.hpp"

namespace py = pybind11;
typedef CellBasedParabolicPdeSolver<3 > CellBasedParabolicPdeSolver3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 4> _boost_numeric_ublas_c_vector_lt_double_4_gt_;
typedef ::boost::numeric::ublas::c_matrix<double, 4, 4> _boost_numeric_ublas_c_matrix_lt_double_4_4_gt_;

class CellBasedParabolicPdeSolver3_Overrides : public CellBasedParabolicPdeSolver3{
    public:
    using CellBasedParabolicPdeSolver3::CellBasedParabolicPdeSolver;
    ::boost::numeric::ublas::c_vector<double, 4> ComputeVectorTerm(::boost::numeric::ublas::c_vector<double, 4> & rPhi, ::boost::numeric::ublas::c_matrix<double, 3, 4> & rGradPhi, ::ChastePoint<3> & rX, ::boost::numeric::ublas::c_vector<double, 1> & rU, ::boost::numeric::ublas::c_matrix<double, 1, 3> & rGradU, ::Element<3, 3> * pElement) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_4_gt_,
            CellBasedParabolicPdeSolver3,
            ComputeVectorTerm,
                    rPhi,
        rGradPhi,
        rX,
        rU,
        rGradU,
        pElement);
    }
    ::boost::numeric::ublas::c_matrix<double, 4, 4> ComputeMatrixTerm(::boost::numeric::ublas::c_vector<double, 4> & rPhi, ::boost::numeric::ublas::c_matrix<double, 3, 4> & rGradPhi, ::ChastePoint<3> & rX, ::boost::numeric::ublas::c_vector<double, 1> & rU, ::boost::numeric::ublas::c_matrix<double, 1, 3> & rGradU, ::Element<3, 3> * pElement) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_matrix_lt_double_4_4_gt_,
            CellBasedParabolicPdeSolver3,
            ComputeMatrixTerm,
                    rPhi,
        rGradPhi,
        rX,
        rU,
        rGradU,
        pElement);
    }
    void ResetInterpolatedQuantities() override {
        PYBIND11_OVERRIDE(
            void,
            CellBasedParabolicPdeSolver3,
            ResetInterpolatedQuantities,
            );
    }
    void IncrementInterpolatedQuantities(double phiI, ::Node<3> const * arg1) override {
        PYBIND11_OVERRIDE(
            void,
            CellBasedParabolicPdeSolver3,
            IncrementInterpolatedQuantities,
                    phiI,
        arg1);
    }

};
void register_CellBasedParabolicPdeSolver3_class(py::module &m){
py::class_<CellBasedParabolicPdeSolver3 , CellBasedParabolicPdeSolver3_Overrides , boost::shared_ptr<CellBasedParabolicPdeSolver3 >   >(m, "CellBasedParabolicPdeSolver3")
        .def(py::init<::TetrahedralMesh<3, 3> *, ::AbstractLinearParabolicPde<3> *, ::BoundaryConditionsContainer<3, 3, 1> * >(), py::arg("pMesh"), py::arg("pPde"), py::arg("pBoundaryConditions"))
    ;
}
