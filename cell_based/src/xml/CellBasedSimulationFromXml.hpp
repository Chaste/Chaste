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

#ifndef CELLBASEDSIMULATIONFROMXML_HPP_
#define CELLBASEDSIMULATIONFROMXML_HPP_

/**
 * @file CellBasedSimulationFromXml.hpp
 *
 * Loader that constructs and returns a ready-to-Solve() off-lattice
 * MeshBasedCellPopulation simulation from an XML specification file.
 *
 * ### XML format
 *
 * ```xml
 * <CellBasedSimulation version="1" elementDim="2" spaceDim="2">
 *   <Geometry>
 *     <HoneycombMeshGenerator numCellsAcross="5" numCellsUp="5" numGhostLayers="0"/>
 *   </Geometry>
 *   <InitialCells>
 *     <DefaultCellCycleModel type="UniformCellCycleModel">
 *       <MinCellCycleDuration>10</MinCellCycleDuration>
 *       <MaxCellCycleDuration>12</MaxCellCycleDuration>
 *     </DefaultCellCycleModel>
 *     <DefaultProliferativeType>TransitCellProliferativeType</DefaultProliferativeType>
 *     <DefaultMutationState>WildTypeCellMutationState</DefaultMutationState>
 *   </InitialCells>
 *   <Population type="MeshBasedCellPopulation">
 *     ...
 *   </Population>
 *   <Forces>
 *     <Force type="GeneralisedLinearSpringForce">...</Force>
 *   </Forces>
 *   <CellKillers/>
 *   <SimulationModifiers/>
 *   <BoundaryConditions/>
 *   <NumericalMethod type="ForwardEulerNumericalMethod">...</NumericalMethod>
 *   <Simulation>
 *     <EndTime>10</EndTime>
 *     ...
 *   </Simulation>
 * </CellBasedSimulation>
 * ```
 *
 * ### Design notes
 * - Parameter element names are identical to those written by
 *   `Output*Parameters()` methods (member variable name minus leading `m`),
 *   so the round-trip is exact once the output path is implemented.
 * - Birth times are drawn uniformly from `[-MaxCellCycleDuration, 0]` to
 *   randomise the initial cell phases (the "RandomCellGenerator" behaviour).
 * - Currently supports 2-D `MeshBasedCellPopulation` with
 *   `UniformCellCycleModel`, optionally with a `GeneralisedLinearSpringForce`.
 *   Support for additional population types and component families will be
 *   added incrementally.
 */

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Exception.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "CellsGenerator.hpp"
#include "UniformCellCycleModel.hpp"
#include "NullSrnModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "RandomNumberGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "ForwardEulerNumericalMethod.hpp"

/**
 * Minimal, dependency-free XML helper.
 *
 * The implementation uses only `std::string` operations.  It is intentionally
 * restricted to the simple, flat (non-recursive) structure of the Chaste
 * simulation-specification XML, and is not a general-purpose XML parser.
 */
namespace CellBasedXmlDetail
{

/**
 * Read the entire content of a file into a string.
 *
 * @param rFilePath  path to the file to read
 * @return file content as a string
 */
inline std::string ReadFile(const std::string& rFilePath)
{
    std::ifstream ifs(rFilePath.c_str());
    if (!ifs.is_open())
    {
        EXCEPTION("CellBasedSimulationFromXml: cannot open file '" + rFilePath + "'");
    }
    std::ostringstream oss;
    oss << ifs.rdbuf();
    return oss.str();
}

/**
 * Strip XML comments (`<!--...-->`) from a string.
 *
 * @param rXml  input XML string (modified in place)
 */
inline void StripComments(std::string& rXml)
{
    std::string result;
    result.reserve(rXml.size());
    std::size_t pos = 0;
    while (pos < rXml.size())
    {
        std::size_t comment_start = rXml.find("<!--", pos);
        if (comment_start == std::string::npos)
        {
            result += rXml.substr(pos);
            break;
        }
        result += rXml.substr(pos, comment_start - pos);
        std::size_t comment_end = rXml.find("-->", comment_start + 4);
        if (comment_end == std::string::npos)
        {
            break; // malformed comment — skip rest
        }
        pos = comment_end + 3;
    }
    rXml = result;
}

/**
 * Trim leading and trailing whitespace from a string.
 *
 * @param rStr  string to trim
 * @return trimmed string
 */
inline std::string Trim(const std::string& rStr)
{
    const std::string ws = " \t\n\r";
    std::size_t first = rStr.find_first_not_of(ws);
    if (first == std::string::npos)
    {
        return "";
    }
    std::size_t last = rStr.find_last_not_of(ws);
    return rStr.substr(first, last - first + 1);
}

/**
 * Extract the value of a named attribute from an opening tag string.
 *
 * E.g., from `<HoneycombMeshGenerator numCellsAcross="5" .../>` with
 * `attrName="numCellsAcross"` this returns `"5"`.
 *
 * @param rTag      the tag string (everything up to and including `>` or `/>`)
 * @param rAttrName name of the attribute to find
 * @param rDefault  value to return if the attribute is not present
 * @return attribute value as a string, or rDefault
 */
inline std::string GetAttr(const std::string& rTag,
                           const std::string& rAttrName,
                           const std::string& rDefault = "")
{
    std::string search = rAttrName + "=\"";
    std::size_t pos = rTag.find(search);
    if (pos == std::string::npos)
    {
        return rDefault;
    }
    pos += search.size();
    std::size_t end = rTag.find('"', pos);
    if (end == std::string::npos)
    {
        return rDefault;
    }
    return rTag.substr(pos, end - pos);
}

/**
 * Find the opening tag for a given element name inside `rXml`.
 *
 * Handles both self-closing (`<Foo .../>`) and normal (`<Foo ...>`) tags.
 *
 * @param rXml        source string to search in
 * @param rTagName    element name to find
 * @param outTagStart position of `<` in rXml, or npos if not found
 * @param outTagEnd   position one-past-end of the opening tag (`>` + 1),
 *                    or npos
 * @param outSelfClose true when the tag is self-closing (`/>`)
 */
inline void FindOpenTag(const std::string& rXml,
                        const std::string& rTagName,
                        std::size_t& outTagStart,
                        std::size_t& outTagEnd,
                        bool& outSelfClose,
                        std::size_t searchFrom = 0)
{
    outTagStart = std::string::npos;
    outTagEnd   = std::string::npos;
    outSelfClose = false;

    std::size_t pos = searchFrom;
    while (pos < rXml.size())
    {
        std::size_t candidate = rXml.find("<" + rTagName, pos);
        if (candidate == std::string::npos)
        {
            return;
        }
        // Verify that the character after the tag name is whitespace, '>' or '/>'
        std::size_t after = candidate + 1 + rTagName.size();
        if (after >= rXml.size())
        {
            return;
        }
        char ch = rXml[after];
        if (ch != ' ' && ch != '\t' && ch != '\n' && ch != '\r'
            && ch != '>' && ch != '/')
        {
            pos = candidate + 1;
            continue; // partial match (e.g. <FooBar when searching for <Foo)
        }

        outTagStart = candidate;
        // Find end of this opening tag
        std::size_t close = rXml.find('>', after);
        if (close == std::string::npos)
        {
            return;
        }
        outTagEnd = close + 1;
        outSelfClose = (rXml[close - 1] == '/');
        return;
    }
}

/**
 * Extract the inner content of the first occurrence of element `rTagName`
 * within `rXml`.
 *
 * For self-closing tags the inner content is an empty string.
 *
 * @param rXml       source string to search in
 * @param rTagName   element name
 * @param outFound   set to true if the element was found
 * @param outOpenTag the full opening tag string (useful for reading attributes)
 * @return inner XML content of the element (empty string if self-closing)
 */
inline std::string GetInner(const std::string& rXml,
                            const std::string& rTagName,
                            bool& outFound,
                            std::string& outOpenTag,
                            std::size_t searchFrom = 0)
{
    outFound = false;
    outOpenTag = "";

    std::size_t tagStart, tagEnd;
    bool selfClose;
    FindOpenTag(rXml, rTagName, tagStart, tagEnd, selfClose, searchFrom);
    if (tagStart == std::string::npos)
    {
        return "";
    }
    outOpenTag = rXml.substr(tagStart, tagEnd - tagStart);
    outFound = true;

    if (selfClose)
    {
        return "";
    }

    // Find matching close tag, accounting for nesting
    std::string openTagStr  = "<"  + rTagName;
    std::string closeTagStr = "</" + rTagName + ">";
    int depth = 1;
    std::size_t searchPos = tagEnd;
    while (depth > 0 && searchPos < rXml.size())
    {
        std::size_t nextOpen  = rXml.find(openTagStr, searchPos);
        std::size_t nextClose = rXml.find(closeTagStr, searchPos);

        if (nextClose == std::string::npos)
        {
            break;
        }
        // Check that nextOpen actually starts a new tag (not a close tag)
        if (nextOpen != std::string::npos && nextOpen < nextClose)
        {
            std::size_t charAfter = nextOpen + openTagStr.size();
            char ch = (charAfter < rXml.size()) ? rXml[charAfter] : 0;
            bool isRealOpen = (ch == ' ' || ch == '\t' || ch == '\n'
                               || ch == '\r' || ch == '>' || ch == '/');
            // Check it's not a close tag (starts with '</')
            bool isClose = (nextOpen + 1 < rXml.size() && rXml[nextOpen + 1] == '/');
            if (isRealOpen && !isClose)
            {
                ++depth;
                searchPos = nextOpen + openTagStr.size();
                continue;
            }
        }
        --depth;
        if (depth == 0)
        {
            return rXml.substr(tagEnd, nextClose - tagEnd);
        }
        searchPos = nextClose + closeTagStr.size();
    }
    return "";
}

/**
 * Convenience overload without outOpenTag.
 */
inline std::string GetInner(const std::string& rXml,
                            const std::string& rTagName,
                            bool& outFound)
{
    std::string dummy;
    return GetInner(rXml, rTagName, outFound, dummy);
}

/**
 * Get the text content of a direct child element, trimmed of whitespace.
 *
 * @param rXml     parent element inner content
 * @param rTag     child element name
 * @param rDefault value to return when element is absent
 * @return trimmed text content, or rDefault
 */
inline std::string ChildText(const std::string& rXml,
                             const std::string& rTag,
                             const std::string& rDefault = "")
{
    bool found;
    std::string inner = GetInner(rXml, rTag, found);
    if (!found)
    {
        return rDefault;
    }
    return Trim(inner);
}

/**
 * Get the text content of a child element as a `double`, returning
 * `defaultVal` if the element is absent.
 */
inline double ChildDouble(const std::string& rXml,
                          const std::string& rTag,
                          double defaultVal)
{
    std::string s = ChildText(rXml, rTag, "");
    if (s.empty())
    {
        return defaultVal;
    }
    return std::stod(s);
}

/**
 * Get the text content of a child element as an `unsigned int`, returning
 * `defaultVal` if the element is absent.
 */
inline unsigned ChildUnsigned(const std::string& rXml,
                              const std::string& rTag,
                              unsigned defaultVal)
{
    std::string s = ChildText(rXml, rTag, "");
    if (s.empty())
    {
        return defaultVal;
    }
    return static_cast<unsigned>(std::stoul(s));
}

/**
 * Get the text content of a child element as a `bool` (0/1 or false/true),
 * returning `defaultVal` if the element is absent.
 */
inline bool ChildBool(const std::string& rXml,
                      const std::string& rTag,
                      bool defaultVal)
{
    std::string s = ChildText(rXml, rTag, "");
    if (s.empty())
    {
        return defaultVal;
    }
    return (s == "1" || s == "true" || s == "True");
}

/**
 * Collect all direct child elements with name `rTag` as a vector of their
 * inner-content strings.
 *
 * @param rXml   parent element inner content
 * @param rTag   element name to collect
 * @return vector of inner-content strings, one per occurrence
 */
inline std::vector<std::string> AllChildren(const std::string& rXml,
                                            const std::string& rTag)
{
    std::vector<std::string> result;
    std::size_t pos = 0;
    while (pos < rXml.size())
    {
        std::size_t tagStart, tagEnd;
        bool selfClose;
        FindOpenTag(rXml, rTag, tagStart, tagEnd, selfClose, pos);
        if (tagStart == std::string::npos)
        {
            break;
        }
        if (selfClose)
        {
            result.push_back("");
            pos = tagEnd;
            continue;
        }
        std::string closeTagStr = "</" + rTag + ">";
        std::size_t closePos = rXml.find(closeTagStr, tagEnd);
        if (closePos == std::string::npos)
        {
            break;
        }
        result.push_back(rXml.substr(tagEnd, closePos - tagEnd));
        pos = closePos + closeTagStr.size();
    }
    return result;
}

} // namespace CellBasedXmlDetail

// ============================================================================

/**
 * Top-level loader: constructs a 2-D off-lattice MeshBasedCellPopulation
 * simulation from an XML specification file and returns a
 * `boost::shared_ptr<OffLatticeSimulation<2> >` ready to call `Solve()`.
 *
 * See file-level documentation for the expected XML format.
 */
class CellBasedSimulationFromXml
{
public:

    /**
     * Parse the XML file at `rXmlPath` and return a configured simulation.
     *
     * Only 2-D simulations (`elementDim="2" spaceDim="2"`) are currently
     * supported.  Other dimensions will throw an `Exception`.
     *
     * @param rXmlPath  absolute or relative path to the XML specification file
     * @return shared pointer to the assembled `OffLatticeSimulation<2>`
     */
    static boost::shared_ptr<OffLatticeSimulation<2> > Load(
        const std::string& rXmlPath)
    {
        using namespace CellBasedXmlDetail;

        // ── 1. Read and pre-process ───────────────────────────────────────
        std::string xml = ReadFile(rXmlPath);
        StripComments(xml);

        // Find root element
        bool foundRoot;
        std::string rootOpenTag;
        std::string rootInner = GetInner(xml, "CellBasedSimulation",
                                         foundRoot, rootOpenTag);
        if (!foundRoot)
        {
            EXCEPTION("CellBasedSimulationFromXml: no <CellBasedSimulation> "
                      "root element found in '" + rXmlPath + "'");
        }

        const std::string elemDimStr = GetAttr(rootOpenTag, "elementDim", "2");
        const std::string spaceDimStr = GetAttr(rootOpenTag, "spaceDim", "2");
        if (elemDimStr != "2" || spaceDimStr != "2")
        {
            EXCEPTION("CellBasedSimulationFromXml: only elementDim=2 / "
                      "spaceDim=2 is currently supported.");
        }

        // ── 2. Geometry ───────────────────────────────────────────────────
        bool foundGeom;
        std::string geomOpenTag;
        std::string geomInner = GetInner(rootInner, "Geometry",
                                         foundGeom, geomOpenTag);
        if (!foundGeom)
        {
            EXCEPTION("CellBasedSimulationFromXml: <Geometry> element "
                      "is required.");
        }

        // Currently: HoneycombMeshGenerator only
        bool foundGen;
        std::string genOpenTag;
        GetInner(geomInner, "HoneycombMeshGenerator", foundGen, genOpenTag);
        if (!foundGen)
        {
            EXCEPTION("CellBasedSimulationFromXml: only "
                      "<HoneycombMeshGenerator> is currently supported inside "
                      "<Geometry>.");
        }
        unsigned numX      = static_cast<unsigned>(
            std::stoul(GetAttr(genOpenTag, "numCellsAcross", "5")));
        unsigned numY      = static_cast<unsigned>(
            std::stoul(GetAttr(genOpenTag, "numCellsUp",     "5")));
        unsigned numGhosts = static_cast<unsigned>(
            std::stoul(GetAttr(genOpenTag, "numGhostLayers", "0")));

        HoneycombMeshGenerator generator(numX, numY, numGhosts);
        boost::shared_ptr<MutableMesh<2,2> > p_mesh = generator.GetMesh();

        // ── 3. Initial cells ──────────────────────────────────────────────
        bool foundIC;
        std::string icInner = GetInner(rootInner, "InitialCells", foundIC);

        // Default parameters for the cell-cycle model
        double minCCDuration = 10.0;
        double maxCCDuration = 12.0;

        if (foundIC)
        {
            bool foundCCM;
            std::string ccmInner = GetInner(icInner,
                                            "DefaultCellCycleModel", foundCCM);
            if (foundCCM)
            {
                minCCDuration = ChildDouble(ccmInner, "MinCellCycleDuration",
                                           minCCDuration);
                maxCCDuration = ChildDouble(ccmInner, "MaxCellCycleDuration",
                                           maxCCDuration);
            }
        }

        // Build the cell vector with randomised birth times
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        RandomNumberGenerator* const p_rng = RandomNumberGenerator::Instance();
        unsigned num_cells = p_mesh->GetNumNodes();
        cells.reserve(num_cells);
        for (unsigned i = 0; i < num_cells; ++i)
        {
            UniformCellCycleModel* p_model = new UniformCellCycleModel();
            p_model->SetMinCellCycleDuration(minCCDuration);
            p_model->SetMaxCellCycleDuration(maxCCDuration);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_transit_type);

            // Random birth time in [-maxCCDuration, 0]
            double birth_time = -maxCCDuration * p_rng->ranf();
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // ── 4. Cell population ────────────────────────────────────────────
        boost::shared_ptr<MeshBasedCellPopulation<2> > p_population(
            new MeshBasedCellPopulation<2>(*p_mesh, cells));

        // Apply population parameters from <Population> element (if present)
        bool foundPop;
        std::string popInner = GetInner(rootInner, "Population", foundPop);
        if (foundPop)
        {
            p_population->SetAreaBasedDampingConstant(
                ChildBool(popInner, "UseAreaBasedDampingConstant", false));
            p_population->SetAreaBasedDampingConstantParameter(
                ChildDouble(popInner, "AreaBasedDampingConstantParameter", 0.1));
            p_population->SetWriteVtkAsPoints(
                ChildBool(popInner, "WriteVtkAsPoints", false));
            p_population->SetBoundVoronoiTessellation(
                ChildBool(popInner, "BoundVoronoiTessellation", false));
            p_population->SetMeinekeDivisionSeparation(
                ChildDouble(popInner, "MeinekeDivisionSeparation", 0.3));
            p_population->SetDampingConstantNormal(
                ChildDouble(popInner, "DampingConstantNormal", 1.0));
            p_population->SetDampingConstantMutant(
                ChildDouble(popInner, "DampingConstantMutant", 1.0));
            p_population->SetOutputResultsForChasteVisualizer(
                ChildBool(popInner, "OutputResultsForChasteVisualizer", true));
        }

        // ── 5. Simulation ─────────────────────────────────────────────────
        bool foundSim;
        std::string simInner = GetInner(rootInner, "Simulation", foundSim);

        double end_time      = foundSim ? ChildDouble(simInner, "EndTime",   10.0) : 10.0;
        double dt            = foundSim ? ChildDouble(simInner, "Dt",        1.0/120.0) : 1.0/120.0;
        unsigned sample_mult = foundSim ? ChildUnsigned(simInner, "SamplingTimestepMultiple", 12u) : 12u;
        std::string output_dir = foundSim
            ? ChildText(simInner, "OutputDirectory", "CellBasedSimulationFromXml")
            : "CellBasedSimulationFromXml";

        boost::shared_ptr<OffLatticeSimulation<2> > p_simulator(
            new OffLatticeSimulation<2>(*p_population));
        p_simulator->SetOutputDirectory(output_dir);
        p_simulator->SetDt(dt);
        p_simulator->SetEndTime(end_time);
        p_simulator->SetSamplingTimestepMultiple(sample_mult);

        if (foundSim)
        {
            p_simulator->SetUpdateCellPopulationRule(
                ChildBool(simInner, "UpdateCellPopulation", true));
            p_simulator->SetOutputDivisionLocations(
                ChildBool(simInner, "OutputDivisionLocations", false));
            p_simulator->SetOutputCellVelocities(
                ChildBool(simInner, "OutputCellVelocities", false));
        }

        // ── 6. Numerical method ───────────────────────────────────────────
        // Default (and currently only supported): ForwardEulerNumericalMethod
        bool foundNM;
        std::string nmOpenTag;
        std::string nmInner = GetInner(rootInner, "NumericalMethod",
                                       foundNM, nmOpenTag);
        {
            MAKE_PTR(ForwardEulerNumericalMethod<2>, p_method);
            if (foundNM)
            {
                p_method->SetUseAdaptiveTimestep(
                    ChildBool(nmInner, "UseAdaptiveTimestep", false));
                p_method->SetUseUpdateNodeLocation(
                    ChildBool(nmInner, "UseUpdateNodeLocation", false));
            }
            p_simulator->SetNumericalMethod(p_method);
        }

        // ── 7. Forces ─────────────────────────────────────────────────────
        bool foundForces;
        std::string forcesInner = GetInner(rootInner, "Forces", foundForces);
        if (foundForces)
        {
            // Collect all <Force ...> child elements
            std::size_t pos = 0;
            while (pos < forcesInner.size())
            {
                std::size_t tagStart, tagEnd;
                bool selfClose;
                CellBasedXmlDetail::FindOpenTag(forcesInner, "Force",
                                               tagStart, tagEnd,
                                               selfClose, pos);
                if (tagStart == std::string::npos)
                {
                    break;
                }
                std::string forceOpenTag = forcesInner.substr(tagStart,
                                                              tagEnd - tagStart);
                std::string forceType = GetAttr(forceOpenTag, "type", "");
                std::string forceInner;
                if (!selfClose)
                {
                    const std::string closeTag = "</Force>";
                    std::size_t closePos = forcesInner.find(closeTag, tagEnd);
                    if (closePos != std::string::npos)
                    {
                        forceInner = forcesInner.substr(tagEnd,
                                                       closePos - tagEnd);
                        pos = closePos + closeTag.size();
                    }
                    else
                    {
                        pos = tagEnd;
                    }
                }
                else
                {
                    pos = tagEnd;
                }

                if (forceType == "GeneralisedLinearSpringForce")
                {
                    MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
                    p_force->SetMeinekeSpringStiffness(
                        ChildDouble(forceInner, "MeinekeSpringStiffness", 15.0));
                    p_force->SetMeinekeDivisionRestingSpringLength(
                        ChildDouble(forceInner, "MeinekeDivisionRestingSpringLength", 0.5));
                    p_force->SetMeinekeSpringGrowthDuration(
                        ChildDouble(forceInner, "MeinekeSpringGrowthDuration", 1.0));
                    bool useCutOff = ChildBool(forceInner, "UseCutOffLength", false);
                    if (useCutOff)
                    {
                        double cutOff = ChildDouble(forceInner, "CutOffLength", 1.5);
                        p_force->SetCutOffLength(cutOff);
                    }
                    p_simulator->AddForce(p_force);
                }
                else if (!forceType.empty())
                {
                    // Unknown force type — warn but continue
                    // (WARN_IF_PARALLEL is not available here; use std::cerr)
                }
            }
        }

        return p_simulator;
    }
};

#endif // CELLBASEDSIMULATIONFROMXML_HPP_
