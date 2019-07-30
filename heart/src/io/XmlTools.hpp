/*

Copyright (c) 2005-2019, University of Oxford.
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

#ifndef XMLTOOLS_HPP_
#define XMLTOOLS_HPP_

#define XSD_CXX11

#include <string>
#include <vector>

#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/util/XercesDefs.hpp> // XMLCh

#include <xsd/cxx/tree/elements.hxx>
#include <xsd/cxx/tree/error-handler.hxx>
#include <xsd/cxx/version.hxx>
#include <xsd/cxx/xml/dom/auto-ptr.hxx>
#include <xsd/cxx/xml/string.hxx>

#ifndef X //Also used in XSD code for XdmfMeshWriter
/**
 * Convenience macro for transcoding C++ strings to Xerces' format.
 * @param str  the string to transcode
 */
#define X(str) xsd::cxx::xml::string(str).c_str()
#endif //X

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
     * @return Xerces convenience object
     */
    static XSD_DOM_AUTO_PTR<xercesc::DOMDocument> ReadXmlFile(
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
     * @return Xerces convenience object
     */
    static XSD_DOM_AUTO_PTR<xercesc::DOMDocument> ReadFileToDomDocument(
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
     * @return Xerces element
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
     * @return Xerces element
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
     * @return all the child elements of the given element.
     *
     * @param pElement  the parent element
     */
    static std::vector<xercesc::DOMElement*> GetChildElements(const xercesc::DOMElement* pElement);

    /**
     * @return all elements matching the given path from this context element.
     *
     * @param pContextElement  the root element to search from
     * @param rPath  where to search.  This should be a '/'-separated path of element names.
     */
    static std::vector<xercesc::DOMElement*> FindElements(const xercesc::DOMElement* pContextElement,
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
    static void FindElements(const xercesc::DOMElement* pContextElement,
                             const std::vector<std::string>& rNames,
                             std::vector<xercesc::DOMElement*>& rResults,
                             unsigned depth=0);

    /**
     * Helper method for URL-escaping spaces in file paths, to avoid confusing Xerces
     * regarding schema locations.  Note that this is a very specific fix: it doesn't
     * do general URL-escaping.
     *
     * @param rPath  the path to escape
     * @return path with spaces escaped
     */
    static std::string EscapeSpaces(const std::string& rPath);
};

#endif /* XMLTOOLS_HPP_ */
