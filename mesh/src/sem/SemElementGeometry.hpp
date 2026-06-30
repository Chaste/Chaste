/*

Copyright (c) 2005-2026, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef SEMELEMENTGEOMETRY_HPP_
#define SEMELEMENTGEOMETRY_HPP_

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <vector>

#include "ChasteCuboid.hpp"
#include "Exception.hpp"
#include "Node.hpp"
#include "UblasVectorInclude.hpp"

template<unsigned DIM>
class SemMesh;

#ifdef CHASTE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1
#include <vtkCellArray.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDelaunay2D.h>
#include <vtkDelaunay3D.h>
#include <vtkIdList.h>
#include <vtkMassProperties.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkTriangleFilter.h>
#endif // CHASTE_VTK

/**
 * Generated surface data for one SEM element.
 *
 * Measure is length in 1D, area in 2D, and volume in 3D.
 */
template<unsigned DIM>
struct SemElementSurface
{
    std::vector<c_vector<double, DIM> > Points;        /**< Surface point coordinates. */
    std::vector<std::array<unsigned, 2> > Lines;       /**< Line segments (pairs of point indices). */
    std::vector<std::array<unsigned, 3> > Triangles;   /**< Triangular faces (triples of point indices). */
    double Measure = 0.0;                              /**< Length (1D), area (2D), or volume (3D). */
    double LocalSpacing = 0.0;                         /**< Characteristic inter-node spacing. */
    double Alpha = 0.0;                                /**< Alpha parameter used for surface generation. */
    double ExpansionRadius = 0.0;                      /**< Outward expansion radius applied to surface points. */
};

/**
 * Utility methods for deriving a SEM element boundary from its point cloud.
 */
template<unsigned DIM>
class SemElementGeometry
{
public:
    /**
     * Generate a surface/boundary for a SEM element.
     *
     * @param rMesh the SEM mesh
     * @param elementIndex the SEM element index
     * @param alphaMultiplier multiplier applied to local node spacing
     * @param expansionMultiplier multiplier applied to local node spacing for outward expansion
     * @param useExpandedSurface whether to expand the generated shape
     * @return generated surface data
     */
    static SemElementSurface<DIM> GenerateSurface(SemMesh<DIM>& rMesh,
                                                  unsigned elementIndex,
                                                  double alphaMultiplier,
                                                  double expansionMultiplier,
                                                  bool useExpandedSurface)
    {
        if (alphaMultiplier <= 0.0)
        {
            EXCEPTION("SEM surface alpha multiplier must be positive");
        }
        if (expansionMultiplier < 0.0)
        {
            EXCEPTION("SEM surface expansion multiplier must be non-negative");
        }

        if constexpr (DIM == 1)
        {
            std::vector<c_vector<double, DIM> > points = GetNodeLocations(rMesh, elementIndex);
            const double local_spacing = CalculateLocalSpacing(rMesh, elementIndex);
            const double alpha = alphaMultiplier * local_spacing;
            const double expansion_radius = useExpandedSurface ? expansionMultiplier * local_spacing : 0.0;
            return GenerateSurface1d(points, local_spacing, alpha, expansion_radius);
        }
        else
        {
            std::vector<c_vector<double, DIM> > points = GetNodeLocations(rMesh, elementIndex);
            EnsureSufficientPoints(points);

            const double local_spacing = CalculateLocalSpacing(rMesh, elementIndex);
            const double alpha = alphaMultiplier * local_spacing;
            const double expansion_radius = useExpandedSurface ? expansionMultiplier * local_spacing : 0.0;
#ifdef CHASTE_VTK
            ExpandPoints(points, rMesh.GetCentroidOfElement(elementIndex), expansion_radius);

            if constexpr (DIM == 2)
            {
                return GenerateSurface2d(points, local_spacing, alpha, expansion_radius);
            }
            else
            {
                return GenerateSurface3d(points, local_spacing, alpha, expansion_radius);
            }
#else
            (void)points;
            (void)local_spacing;
            (void)alpha;
            (void)expansion_radius;
            EXCEPTION("SEM alpha-shape surfaces in 2D and 3D require Chaste to be compiled with VTK");
#endif // CHASTE_VTK
        }
    }

    /**
     * Estimate local node spacing from the element bounding box and node count.
     *
     * @param rMesh the SEM mesh
     * @param elementIndex the SEM element index
     * @return characteristic local node spacing
     */
    static double CalculateLocalSpacing(SemMesh<DIM>& rMesh, unsigned elementIndex)
    {
        auto p_element = rMesh.GetElement(elementIndex);
        const unsigned num_nodes = p_element->GetNumNodes();
        if (num_nodes < 2u)
        {
            EXCEPTION("A SEM element must contain at least two distinct nodes to calculate local spacing");
        }

        ChasteCuboid<DIM> bounding_box = rMesh.CalculateBoundingBoxOfElement(elementIndex);

        if constexpr (DIM == 1)
        {
            const double width = bounding_box.GetWidth(0u);
            if (width <= 0.0)
            {
                EXCEPTION("A SEM element must contain at least two distinct nodes to calculate local spacing");
            }
            return width / static_cast<double>(num_nodes - 1u);
        }
        else if constexpr (DIM == 2)
        {
            const double area = bounding_box.GetWidth(0u) * bounding_box.GetWidth(1u);
            const double denominator = std::sqrt(static_cast<double>(num_nodes)) - 1.0;
            if (area <= 0.0 || denominator <= 0.0)
            {
                EXCEPTION("A SEM element must contain at least two distinct nodes to calculate local spacing");
            }
            return std::sqrt(area) / denominator;
        }
        else
        {
            const double volume = bounding_box.GetWidth(0u) * bounding_box.GetWidth(1u) * bounding_box.GetWidth(2u);
            const double denominator = std::cbrt(static_cast<double>(num_nodes)) - 1.0;
            if (volume <= 0.0 || denominator <= 0.0)
            {
                EXCEPTION("A SEM element must contain at least two distinct nodes to calculate local spacing");
            }
            return std::cbrt(volume) / denominator;
        }
    }

private:
    /**
     * Collect node positions of a SEM element into a vector.
     *
     * @param rMesh the SEM mesh
     * @param elementIndex the SEM element index
     * @return vector of node locations for the element
     */
    static std::vector<c_vector<double, DIM> > GetNodeLocations(SemMesh<DIM>& rMesh, unsigned elementIndex)
    {
        auto p_element = rMesh.GetElement(elementIndex);
        std::vector<c_vector<double, DIM> > points;
        points.reserve(p_element->GetNumNodes());
        for (unsigned node_index = 0u; node_index < p_element->GetNumNodes(); ++node_index)
        {
            points.push_back(p_element->GetNode(node_index)->rGetLocation());
        }
        return points;
    }

    /**
     * Expand each point outward from the centroid by a fixed radius.
     *
     * @param rPoints the point cloud to modify in place
     * @param rCentroid centroid of the point cloud
     * @param expansionRadius outward displacement to apply (no-op if <= 0)
     */
    static void ExpandPoints(std::vector<c_vector<double, DIM> >& rPoints,
                             const c_vector<double, DIM>& rCentroid,
                             double expansionRadius)
    {
        if (expansionRadius <= 0.0)
        {
            return;
        }

        for (auto& r_point : rPoints)
        {
            c_vector<double, DIM> direction = r_point - rCentroid;
            const double length = std::sqrt(inner_prod(direction, direction));
            if (length > 1e-12)
            {
                r_point += expansionRadius * direction / length;
            }
        }
    }

    /**
     * Generate a surface for a 1D SEM element (an interval with two endpoints).
     *
     * @param rPoints node positions of the element
     * @param localSpacing characteristic inter-node spacing
     * @param alpha alpha parameter (stored in the result)
     * @param expansionRadius outward expansion applied to the interval endpoints
     * @return generated surface data for the 1D element
     */
    static SemElementSurface<DIM> GenerateSurface1d(std::vector<c_vector<double, DIM> >& rPoints,
                                                    double localSpacing,
                                                    double alpha,
                                                    double expansionRadius)
    {
        if (rPoints.size() < 2u)
        {
            EXCEPTION("A 1D SEM element surface requires at least two distinct nodes");
        }

        double min_x = rPoints[0][0];
        double max_x = rPoints[0][0];
        for (const auto& r_point : rPoints)
        {
            min_x = std::min(min_x, r_point[0]);
            max_x = std::max(max_x, r_point[0]);
        }

        min_x -= expansionRadius;
        max_x += expansionRadius;
        if (max_x <= min_x)
        {
            EXCEPTION("A 1D SEM element surface requires a non-zero interval");
        }

        SemElementSurface<DIM> surface;
        surface.Points.resize(2u);
        surface.Points[0][0] = min_x;
        surface.Points[1][0] = max_x;
        surface.Lines.push_back({{0u, 1u}});
        surface.Measure = max_x - min_x;
        surface.LocalSpacing = localSpacing;
        surface.Alpha = alpha;
        surface.ExpansionRadius = expansionRadius;
        return surface;
    }

    /**
     * Check that the point cloud has enough points for alpha-shape surface generation.
     * Throws an exception if there are too few points (< 3 in 2D, < 4 in 3D).
     *
     * @param rPoints the point cloud to validate
     */
    static void EnsureSufficientPoints(const std::vector<c_vector<double, DIM> >& rPoints)
    {
        if constexpr (DIM == 2)
        {
            if (rPoints.size() < 3u)
            {
                EXCEPTION("A 2D SEM alpha-shape surface requires at least three distinct nodes");
            }
        }
        else if constexpr (DIM == 3)
        {
            if (rPoints.size() < 4u)
            {
                EXCEPTION("A 3D SEM alpha-shape surface requires at least four distinct nodes");
            }
        }
    }

#ifdef CHASTE_VTK
    static vtkSmartPointer<vtkPolyData> MakeInputPolyData(const std::vector<c_vector<double, DIM> >& rPoints)
    {
        vtkSmartPointer<vtkPoints> p_points = vtkSmartPointer<vtkPoints>::New();
        p_points->SetDataTypeToDouble();

        for (const auto& r_point : rPoints)
        {
            if constexpr (DIM == 2)
            {
                p_points->InsertNextPoint(r_point[0], r_point[1], 0.0);
            }
            else
            {
                p_points->InsertNextPoint(r_point[0], r_point[1], r_point[2]);
            }
        }

        vtkSmartPointer<vtkPolyData> p_poly_data = vtkSmartPointer<vtkPolyData>::New();
        p_poly_data->SetPoints(p_points);
        return p_poly_data;
    }

    static void CopyVtkPoints(vtkPoints* pVtkPoints, SemElementSurface<DIM>& rSurface)
    {
        if (pVtkPoints == nullptr || pVtkPoints->GetNumberOfPoints() == 0)
        {
            EXCEPTION("SEM alpha-shape generation produced no surface points");
        }

        rSurface.Points.reserve(pVtkPoints->GetNumberOfPoints());
        for (vtkIdType i = 0; i < pVtkPoints->GetNumberOfPoints(); ++i)
        {
            double point[3];
            pVtkPoints->GetPoint(i, point);
            c_vector<double, DIM> location;
            for (unsigned dim = 0u; dim < DIM; ++dim)
            {
                location[dim] = point[dim];
            }
            rSurface.Points.push_back(location);
        }
    }

    static SemElementSurface<DIM> GenerateSurface2d(const std::vector<c_vector<double, DIM> >& rPoints,
                                                    double localSpacing,
                                                    double alpha,
                                                    double expansionRadius)
    {
        vtkSmartPointer<vtkPolyData> p_input = MakeInputPolyData(rPoints);
        vtkSmartPointer<vtkDelaunay2D> p_delaunay = vtkSmartPointer<vtkDelaunay2D>::New();
        p_delaunay->SetAlpha(alpha);
        p_delaunay->BoundingTriangulationOff();
        p_delaunay->SetInputData(p_input);
        p_delaunay->Update();

        vtkSmartPointer<vtkTriangleFilter> p_triangle_filter = vtkSmartPointer<vtkTriangleFilter>::New();
        p_triangle_filter->SetInputData(p_delaunay->GetOutput());
        p_triangle_filter->Update();

        vtkPolyData* p_output = p_triangle_filter->GetOutput();

        SemElementSurface<DIM> surface;
        surface.LocalSpacing = localSpacing;
        surface.Alpha = alpha;
        surface.ExpansionRadius = expansionRadius;
        CopyVtkPoints(p_output->GetPoints(), surface);

        vtkSmartPointer<vtkIdList> p_cell_point_ids = vtkSmartPointer<vtkIdList>::New();
        vtkCellArray* p_polys = p_output->GetPolys();
        p_polys->InitTraversal();

        std::map<std::pair<unsigned, unsigned>, unsigned> edge_counts;
        while (p_polys->GetNextCell(p_cell_point_ids))
        {
            if (p_cell_point_ids->GetNumberOfIds() != 3)
            {
                continue;
            }

            const unsigned a = static_cast<unsigned>(p_cell_point_ids->GetId(0));
            const unsigned b = static_cast<unsigned>(p_cell_point_ids->GetId(1));
            const unsigned c = static_cast<unsigned>(p_cell_point_ids->GetId(2));

            const auto& r_a = surface.Points[a];
            const auto& r_b = surface.Points[b];
            const auto& r_c = surface.Points[c];
            surface.Measure += 0.5 * std::fabs((r_b[0] - r_a[0])*(r_c[1] - r_a[1])
                                                - (r_b[1] - r_a[1])*(r_c[0] - r_a[0]));

            const std::array<std::pair<unsigned, unsigned>, 3> edges =
            {{
                std::minmax(a, b),
                std::minmax(b, c),
                std::minmax(c, a)
            }};
            for (const auto& r_edge : edges)
            {
                edge_counts[r_edge]++;
            }
        }

        for (const auto& r_edge_count : edge_counts)
        {
            if (r_edge_count.second == 1u)
            {
                surface.Lines.push_back({{r_edge_count.first.first, r_edge_count.first.second}});
            }
        }

        if (surface.Lines.empty() || surface.Measure <= 0.0)
        {
            EXCEPTION("SEM alpha-shape generation failed to produce a 2D boundary");
        }

        return surface;
    }

    static SemElementSurface<DIM> GenerateSurface3d(const std::vector<c_vector<double, DIM> >& rPoints,
                                                    double localSpacing,
                                                    double alpha,
                                                    double expansionRadius)
    {
        vtkSmartPointer<vtkPolyData> p_input = MakeInputPolyData(rPoints);
        vtkSmartPointer<vtkDelaunay3D> p_delaunay = vtkSmartPointer<vtkDelaunay3D>::New();
        p_delaunay->SetAlpha(alpha);
        p_delaunay->AlphaTetsOn();
        p_delaunay->AlphaTrisOff();
        p_delaunay->AlphaLinesOff();
        p_delaunay->AlphaVertsOff();
        p_delaunay->SetInputData(p_input);
        p_delaunay->Update();

        vtkSmartPointer<vtkDataSetSurfaceFilter> p_surface_filter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
        p_surface_filter->SetInputData(p_delaunay->GetOutput());
        p_surface_filter->Update();

        vtkSmartPointer<vtkTriangleFilter> p_triangle_filter = vtkSmartPointer<vtkTriangleFilter>::New();
        p_triangle_filter->SetInputData(p_surface_filter->GetOutput());
        p_triangle_filter->Update();

        vtkPolyData* p_output = p_triangle_filter->GetOutput();

        SemElementSurface<DIM> surface;
        surface.LocalSpacing = localSpacing;
        surface.Alpha = alpha;
        surface.ExpansionRadius = expansionRadius;
        CopyVtkPoints(p_output->GetPoints(), surface);

        vtkSmartPointer<vtkIdList> p_cell_point_ids = vtkSmartPointer<vtkIdList>::New();
        vtkCellArray* p_polys = p_output->GetPolys();
        p_polys->InitTraversal();
        while (p_polys->GetNextCell(p_cell_point_ids))
        {
            if (p_cell_point_ids->GetNumberOfIds() == 3)
            {
                surface.Triangles.push_back({{
                    static_cast<unsigned>(p_cell_point_ids->GetId(0)),
                    static_cast<unsigned>(p_cell_point_ids->GetId(1)),
                    static_cast<unsigned>(p_cell_point_ids->GetId(2))
                }});
            }
        }

        if (surface.Triangles.empty())
        {
            EXCEPTION("SEM alpha-shape generation failed to produce a 3D surface");
        }

        vtkSmartPointer<vtkMassProperties> p_mass_properties = vtkSmartPointer<vtkMassProperties>::New();
        p_mass_properties->SetInputData(p_output);
        surface.Measure = p_mass_properties->GetVolume();

        if (surface.Measure <= 0.0)
        {
            EXCEPTION("SEM alpha-shape generation failed to produce a positive 3D volume");
        }

        return surface;
    }
#endif // CHASTE_VTK
};

#endif /*SEMELEMENTGEOMETRY_HPP_*/
