#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <petsc/private/vecimpl.h>
#include <petsc/private/matimpl.h>
#include "PythonPetscObjectConverters.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellBasedEllipticPdeSolver.hpp"

#include "CellBasedEllipticPdeSolver3.cppwg.hpp"

namespace py = pybind11;
typedef CellBasedEllipticPdeSolver<3 > CellBasedEllipticPdeSolver3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 4> _boost_numeric_ublas_c_vector_lt_double_4_gt_;
typedef ::boost::numeric::ublas::c_matrix<double, 4, 4> _boost_numeric_ublas_c_matrix_lt_double_4_4_gt_;

class CellBasedEllipticPdeSolver3_Overrides : public CellBasedEllipticPdeSolver3{
    public:
    using CellBasedEllipticPdeSolver3::CellBasedEllipticPdeSolver;
    ::boost::numeric::ublas::c_vector<double, 4> ComputeVectorTerm(::boost::numeric::ublas::c_vector<double, 4> & rPhi, ::boost::numeric::ublas::c_matrix<double, 3, 4> & rGradPhi, ::ChastePoint<3> & rX, ::boost::numeric::ublas::c_vector<double, 1> & rU, ::boost::numeric::ublas::c_matrix<double, 1, 3> & rGradU, ::Element<3, 3> * pElement) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_4_gt_,
            CellBasedEllipticPdeSolver3,
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
            CellBasedEllipticPdeSolver3,
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
            CellBasedEllipticPdeSolver3,
            ResetInterpolatedQuantities,
            );
    }
    void IncrementInterpolatedQuantities(double phiI, ::Node<3> const * pNode) override {
        PYBIND11_OVERRIDE(
            void,
            CellBasedEllipticPdeSolver3,
            IncrementInterpolatedQuantities,
                    phiI,
        pNode);
    }
    void InitialiseForSolve(::Vec initialSolution) override {
        PYBIND11_OVERRIDE(
            void,
            CellBasedEllipticPdeSolver3,
            InitialiseForSolve,
                    initialSolution);
    }

};
void register_CellBasedEllipticPdeSolver3_class(py::module &m){
py::class_<CellBasedEllipticPdeSolver3 , CellBasedEllipticPdeSolver3_Overrides , boost::shared_ptr<CellBasedEllipticPdeSolver3 >   >(m, "CellBasedEllipticPdeSolver3")
        .def(py::init<::TetrahedralMesh<3, 3> *, ::AbstractLinearEllipticPde<3, 3> *, ::BoundaryConditionsContainer<3, 3, 1> * >(), py::arg("pMesh"), py::arg("pPde"), py::arg("pBoundaryConditions"))
    ;
}
