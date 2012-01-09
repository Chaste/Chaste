/*

Copyright (C) University of Oxford, 2005-2012

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef XMLTOOLS_HPP_
#define XMLTOOLS_HPP_

#include <string>
#include <vector>

#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/util/XercesDefs.hpp> // XMLCh

#include <xsd/cxx/version.hxx>
#include <xsd/cxx/xml/string.hxx>
#include <xsd/cxx/xml/dom/auto-ptr.hxx>
#include <xsd/cxx/tree/error-handler.hxx>
#include <xsd/cxx/tree/elements.hxx>

/**
 * Convenience macro for transcoding C++ strings to Xerces' format.
 * @param str  the string to transcode
 */
#define X(str) xsd::cxx::xml::string(str).c_str()

/**
 * Convenience macro for transcoding an XML string to a C++ std::string.
 * @param str  the string to transcode
 */
#define X2C(str) xsd::cxx::xml::transcode<char>(str)

/**
 * A class of utility methods for processing XML files, using Xerces and CodeSynthesis XSD.
 */
class XmlTools
{
public:
    /**
     * Read an XML file into a DOM document, turning parsing errors into Chaste Exceptions.
     * Handles initialising the Xerces runtime.
     *
     * @param rFileName  the file to read
     * @param rProps  properties that specify fixed schema locations, if wanted
     * @param validate  whether to perform schema validation
     */
    static xsd::cxx::xml::dom::auto_ptr<xercesc::DOMDocument> ReadXmlFile(
        const std::string& rFileName,
        const ::xsd::cxx::tree::properties<char>& rProps,
        bool validate=true);

    /**
     * Must be called after you have finished working with a document returned by the ReadXmlFile methods.
     * An alternative is to instantiate
     * \code
     * XmlTools::Finalizer finalizer(false);
     * \endcode
     * just before calling ReadXmlFile, provided that you will do all your processing within that scope.
     * The finalizer object will call Finalize in its destructor.
     */
    static void Finalize();

    /**
     * A little class that automatically finalizes Xerces in its destructor.
     */
    class Finalizer
    {
    public:
        /**
         * Create the object.
         * @param init  whether to initialize the Xerces runtime also
         */
        Finalizer(bool init);

        /**
         * Finalize the Xerces runtime.
         */
        ~Finalizer();
    };

    /**
     * Read an XML file into a DOM document.
     * Useful for figuring out what version of the parameters file we're dealing with,
     * so we can construct the right version of the object model.
     *
     * Based on http://wiki.codesynthesis.com/Tree/FAQ#How_do_I_parse_an_XML_document_to_a_Xerces-C.2B.2B_DOM_document.3F
     *
     * Requires the Xerces runtime to have been initialised.
     *
     * @param rFileName  the file to read
     * @param rErrorHandler  handler for any parsing errors
     * @param rProps  properties that specify fixed schema locations, if wanted
     * @param validate  whether to perform schema validation
     */
    static xsd::cxx::xml::dom::auto_ptr<xercesc::DOMDocument> ReadFileToDomDocument(
        const std::string& rFileName,
        ::xsd::cxx::xml::error_handler<char>& rErrorHandler,
        const ::xsd::cxx::tree::properties<char>& rProps,
        bool validate=true);

    /**
     * Display key info about an XML node for debugging.
     *
     * @param rMsg  message to prepend to the report
     * @param pNode  the node to display
     * @param showChildren  whether to recursive display the node's children
     */
    static void PrintNode(const std::string& rMsg, xercesc::DOMNode* pNode, bool showChildren=false);

    /**
     * Fake having a namespace in older configuration files, by setting the namespace
     * on each element in a tree.
     *
     * Based on http://wiki.codesynthesis.com/Tree/FAQ#How_do_I_parse_an_XML_document_that_is_missing_namespace_information.3F
     *
     * @param pDocument  the DOM document containing the tree to be transformed
     * @param pElement  the root of the tree to be transformed
     * @param rNamespace  the namespace to put elements in
     */
    static xercesc::DOMElement* SetNamespace(xercesc::DOMDocument* pDocument,
                                             xercesc::DOMElement* pElement,
                                             const std::string& rNamespace);

    /**
     * Fake having a namespace in older configuration files, by setting the namespace
     * on each element in a tree.
     *
     * Based on http://wiki.codesynthesis.com/Tree/FAQ#How_do_I_parse_an_XML_document_that_is_missing_namespace_information.3F
     *
     * @param pDocument  the DOM document containing the tree to be transformed
     * @param pElement  the root of the tree to be transformed
     * @param pNamespace  the namespace to put elements in
     */
    static xercesc::DOMElement* SetNamespace(xercesc::DOMDocument* pDocument,
                                             xercesc::DOMElement* pElement,
                                             const XMLCh* pNamespace);

    /**
     * Wrap the content (children) of an element within a new element.  The
     * new element becomes the sole child of the original element.
     *
     * @note Doesn't transfer attributes.
     *
     * @param pDocument  the DOM document containing the tree to be transformed
     * @param pElement  the element whose content is to be wrapped
     * @param pNewElementLocalName  the local name (i.e. without namespace prefix) of the wrapping element
     *   (the namespace of pElement will be used).
     */
    static void WrapContentInElement(xercesc::DOMDocument* pDocument,
                                     xercesc::DOMElement* pElement,
                                     const XMLCh* pNewElementLocalName);

    /**
     * Get all the child elements of the given element.
     *
     * @param pElement  the parent element
     */
    static std::vector<xercesc::DOMElement*> GetChildElements(xercesc::DOMElement* pElement);

    /**
     * Find all elements matching the given path from this context element.
     *
     * @param pContextElement  the root element to search from
     * @param rPath  where to search.  This should be a '/'-separated path of element names.
     */
    static std::vector<xercesc::DOMElement*> FindElements(xercesc::DOMElement* pContextElement,
                                                          const std::string& rPath);

    /**
     * Find all elements matching the given path from this context element.
     *
     * @param pContextElement  the root element to search from
     * @param rNames  a list of element names, the first of which is looked for as children of
     *   pContextElement; the next as children of those, etc.
     * @param rResults  vector to be filled in with matching elements
     * @param depth  for managing recursion; should not be provided by users
     */
    static void FindElements(xercesc::DOMElement* pContextElement,
                             const std::vector<std::string>& rNames,
                             std::vector<xercesc::DOMElement*>& rResults,
                             unsigned depth=0);

    /**
     * Helper method for URL-escaping spaces in file paths, to avoid confusing Xerces
     * regarding schema locations.  Note that this is a very specific fix: it doesn't
     * do general URL-escaping.
     *
     * @param rPath  the path to escape
     */
    static std::string EscapeSpaces(const std::string& rPath);

};

#endif /* XMLTOOLS_HPP_ */
