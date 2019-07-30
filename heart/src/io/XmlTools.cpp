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

#include "XmlTools.hpp"

#include <iostream>

#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/QName.hpp>
#include <xercesc/util/XMLUniDefs.hpp> // chLatin_*
#include <xercesc/framework/Wrapper4InputSource.hpp>
#include <xercesc/validators/common/Grammar.hpp>

#include <xsd/cxx/xml/sax/std-input-source.hxx>
#include <xsd/cxx/xml/dom/bits/error-handler-proxy.hxx>
#include <xsd/cxx/tree/exceptions.hxx>

#include "Exception.hpp"

XSD_DOM_AUTO_PTR<xercesc::DOMDocument> XmlTools::ReadXmlFile(
    const std::string& rFileName,
    const ::xsd::cxx::tree::properties<char>& rProps,
    bool validate)
{
    XSD_DOM_AUTO_PTR<xercesc::DOMDocument> p_doc;
    try
    {
        // Initialise Xerces
        xercesc::XMLPlatformUtils::Initialize();
        // Set up an error handler
        ::xsd::cxx::tree::error_handler<char> error_handler;
        // Parse XML to DOM
        p_doc = XmlTools::ReadFileToDomDocument(rFileName, error_handler, rProps, validate);
        // Any errors?
        error_handler.throw_if_failed< ::xsd::cxx::tree::parsing<char> >();
    }
    catch (const ::xsd::cxx::tree::parsing<char>& e)
    {
        Finalize();
        // Test for missing schema/xml file
#if (XSD_INT_VERSION >= 3000000L)
        const ::xsd::cxx::tree::diagnostics<char>& diags = e.diagnostics();
        const ::xsd::cxx::tree::error<char>& first_error = diags[0];
#else
        const ::xsd::cxx::tree::errors<char>& errors = e.errors();
        const ::xsd::cxx::tree::error<char>& first_error = errors[0];
#endif
        if (first_error.line() == 0u)
        {
            std::cerr << first_error << std::endl;
            EXCEPTION("Missing file parsing configuration file: " + rFileName);
        }
        else
        {
            std::cerr << e << std::endl;
            EXCEPTION("XML parsing error in configuration file: " + rFileName);
        }
    }
// LCOV_EXCL_START
    catch (...)
    { // This shouldn't happen, but just in case...
        Finalize();
        throw;
    }
// LCOV_EXCL_STOP
    return p_doc;
}


void XmlTools::Finalize()
{
    xercesc::XMLPlatformUtils::Terminate();
}

XmlTools::Finalizer::Finalizer(bool init)
{
    // The init=true case will very rarely be used, but a parameter to the constructor is needed
    // to stop some compilers complaining about an unused variable!
    if (init)
    {
// LCOV_EXCL_START
        xercesc::XMLPlatformUtils::Initialize();
// LCOV_EXCL_STOP
    }
}

XmlTools::Finalizer::~Finalizer()
{
    XmlTools::Finalize();
}

XSD_DOM_AUTO_PTR<xercesc::DOMDocument> XmlTools::ReadFileToDomDocument(
        const std::string& rFileName,
        ::xsd::cxx::xml::error_handler<char>& rErrorHandler,
        const ::xsd::cxx::tree::properties<char>& rProps,
        bool validate)
{
    using namespace xercesc;
    namespace xml = xsd::cxx::xml;

    // Get an implementation of the Load-Store (LS) interface.
    const XMLCh ls_id [] = {chLatin_L, chLatin_S, chNull};
    DOMImplementation* p_impl(DOMImplementationRegistry::getDOMImplementation(ls_id));

#if _XERCES_VERSION >= 30000
    // Xerces-C++ 3.0.0 and later.
    XSD_DOM_AUTO_PTR<DOMLSParser> p_parser(p_impl->createLSParser(DOMImplementationLS::MODE_SYNCHRONOUS, 0));
    DOMConfiguration* p_conf(p_parser->getDomConfig());

    // Discard comment nodes in the document.
    p_conf->setParameter(XMLUni::fgDOMComments, false);

    // Enable datatype normalization.
    p_conf->setParameter(XMLUni::fgDOMDatatypeNormalization, true);

    // Do not create EntityReference nodes in the DOM tree.  No
    // EntityReference nodes will be created, only the nodes
    // corresponding to their fully expanded substitution text
    // will be created.
    p_conf->setParameter(XMLUni::fgDOMEntities, false);

    // Perform namespace processing.
    p_conf->setParameter(XMLUni::fgDOMNamespaces, true);

    // Do not include ignorable whitespace in the DOM tree.
    p_conf->setParameter(XMLUni::fgDOMElementContentWhitespace, false);

    // Enable validation.
    if (validate)
    {
        p_conf->setParameter(XMLUni::fgDOMValidate, true);
        p_conf->setParameter(XMLUni::fgXercesSchema, true);
        p_conf->setParameter(XMLUni::fgXercesSchemaFullChecking, false);
        // Code taken from xsd/cxx/xml/dom/parsing-source.txx
        if (!rProps.schema_location().empty())
        {
            xml::string locn(rProps.schema_location());
            const void* p_locn(locn.c_str());
            p_conf->setParameter(XMLUni::fgXercesSchemaExternalSchemaLocation,
                                 const_cast<void*>(p_locn));
        }
        if (!rProps.no_namespace_schema_location().empty())
        {
            xml::string locn(rProps.no_namespace_schema_location());
            const void* p_locn(locn.c_str());

            p_conf->setParameter(XMLUni::fgXercesSchemaExternalNoNameSpaceSchemaLocation,
                                 const_cast<void*>(p_locn));
        }
    }
    else
    {
        // This branch is only used by projects
// LCOV_EXCL_START
        p_conf->setParameter(XMLUni::fgDOMValidate, false);
        p_conf->setParameter(XMLUni::fgXercesSchema, false);
        p_conf->setParameter(XMLUni::fgXercesSchemaFullChecking, false);
// LCOV_EXCL_STOP
    }

    // We will release the DOM document ourselves.
    p_conf->setParameter(XMLUni::fgXercesUserAdoptsDOMDocument, true);

    // Set error handler.
    xml::dom::bits::error_handler_proxy<char> ehp(rErrorHandler);
    p_conf->setParameter(XMLUni::fgDOMErrorHandler, &ehp);

#else // _XERCES_VERSION < 30000
    // Same as above but for Xerces-C++ 2 series.
    XSD_DOM_AUTO_PTR<DOMBuilder> p_parser(p_impl->createDOMBuilder(DOMImplementationLS::MODE_SYNCHRONOUS, 0));

    p_parser->setFeature(XMLUni::fgDOMComments, false);
    p_parser->setFeature(XMLUni::fgDOMDatatypeNormalization, true);
    p_parser->setFeature(XMLUni::fgDOMEntities, false);
    p_parser->setFeature(XMLUni::fgDOMNamespaces, true);
    p_parser->setFeature(XMLUni::fgDOMWhitespaceInElementContent, false);
    p_parser->setFeature(XMLUni::fgXercesUserAdoptsDOMDocument, true);

    // Code taken from xsd/cxx/xml/dom/parsing-source.txx
    if (validate)
    {
        p_parser->setFeature(XMLUni::fgDOMValidation, true);
        p_parser->setFeature(XMLUni::fgXercesSchema, true);
        p_parser->setFeature(XMLUni::fgXercesSchemaFullChecking, false);
        if (!rProps.schema_location().empty())
        {
            xml::string locn(rProps.schema_location());
            const void* p_locn(locn.c_str());
            p_parser->setProperty(XMLUni::fgXercesSchemaExternalSchemaLocation,
                                  const_cast<void*>(p_locn));
        }

        if (!rProps.no_namespace_schema_location().empty())
        {
            xml::string locn(rProps.no_namespace_schema_location());
            const void* p_locn(locn.c_str());

            p_parser->setProperty(XMLUni::fgXercesSchemaExternalNoNameSpaceSchemaLocation,
                                  const_cast<void*>(p_locn));
        }
    }
    else
    {
        // This branch is only used by projects
// LCOV_EXCL_START
        p_parser->setFeature(XMLUni::fgDOMValidation, false);
        p_parser->setFeature(XMLUni::fgXercesSchema, false);
        p_parser->setFeature(XMLUni::fgXercesSchemaFullChecking, false);
// LCOV_EXCL_STOP
    }

    xml::dom::bits::error_handler_proxy<char> ehp(rErrorHandler);
    p_parser->setErrorHandler(&ehp);

#endif // _XERCES_VERSION >= 30000

    // Do the parse
    XSD_DOM_AUTO_PTR<DOMDocument> p_doc(p_parser->parseURI(rFileName.c_str()));

    if (ehp.failed())
    {
        p_doc.reset();
    }

    return p_doc;
}

// LCOV_EXCL_START
void XmlTools::PrintNode(const std::string& rMsg, xercesc::DOMNode* pNode, bool showChildren)
{
    std::string prefix = X2C(pNode->getPrefix());
    std::string name = X2C(pNode->getLocalName());
    std::string nsuri = X2C(pNode->getNamespaceURI());
    std::cout << rMsg << " " << pNode << " " << prefix << ":" << name << " in " << nsuri << std::endl;
    if (showChildren)
    {
        for (xercesc::DOMNode* p_node = pNode->getFirstChild();
             p_node != NULL;
             p_node = p_node->getNextSibling())
        {
            std::cout << "     child type " << p_node->getNodeType();
            PrintNode("", p_node, false);
        }
        xercesc::DOMNamedNodeMap* p_attrs = pNode->getAttributes();
        if (p_attrs)
        {
            for (XMLSize_t i=0; i<p_attrs->getLength(); i++)
            {
                 xercesc::DOMNode* p_attr = p_attrs->item(i);
                 std::string value = X2C(p_attr->getNodeValue());
                 PrintNode("     attr (" + value + ")", p_attr, false);
            }
        }
    }
}
// LCOV_EXCL_STOP

xercesc::DOMElement* XmlTools::SetNamespace(xercesc::DOMDocument* pDocument,
                                            xercesc::DOMElement* pElement,
                                            const XMLCh* pNamespace)
{
    using namespace xercesc;

    //PrintNode("Renaming", pElement, true);
    DOMNamedNodeMap* p_orig_attrs = pElement->getAttributes();
    std::vector<std::string> attr_values;
    if (p_orig_attrs)
    {
        for (XMLSize_t i=0; i<p_orig_attrs->getLength(); i++)
        {
            DOMNode* p_attr = p_orig_attrs->item(i);
            attr_values.push_back(X2C(p_attr->getNodeValue()));
        }
    }
    DOMElement* p_new_elt = static_cast<DOMElement*>(
        pDocument->renameNode(pElement, pNamespace, pElement->getLocalName()));
    //PrintNode("   to", p_new_elt, true);
    // Fix attributes - some get broken by the rename!
    if (p_orig_attrs)
    {
        DOMNamedNodeMap* p_new_attrs = p_new_elt->getAttributes();
        assert(p_new_attrs);
        assert(p_new_attrs == p_orig_attrs);
        assert(p_new_attrs->getLength() == attr_values.size());
        for (XMLSize_t i=0; i<p_new_attrs->getLength(); i++)
        {
            DOMNode* p_attr = p_new_attrs->item(i);
            p_attr->setNodeValue(X(attr_values[i]));
        }
    }
    //PrintNode("   after attr fix", p_new_elt, true);

    std::vector<DOMElement*> children = GetChildElements(p_new_elt);
    for (std::vector<DOMElement*>::iterator it = children.begin(); it != children.end(); ++it)
    {
        SetNamespace(pDocument, *it, pNamespace);
    }

    return p_new_elt;
}

xercesc::DOMElement* XmlTools::SetNamespace(xercesc::DOMDocument* pDocument,
                                            xercesc::DOMElement* pElement,
                                            const std::string& rNamespace)
{
    return SetNamespace(pDocument, pElement, X(rNamespace));
}


std::vector<xercesc::DOMElement*> XmlTools::GetChildElements(const xercesc::DOMElement* pElement)
{
    std::vector<xercesc::DOMElement*> children;
    for (xercesc::DOMNode* p_node = pElement->getFirstChild();
         p_node != NULL;
         p_node = p_node->getNextSibling())
    {
        if (p_node->getNodeType() == xercesc::DOMNode::ELEMENT_NODE)
        {
            children.push_back(static_cast<xercesc::DOMElement*>(p_node));
        }
    }
    return children;
}


void XmlTools::FindElements(const xercesc::DOMElement* pContextElement,
                            const std::vector<std::string>& rNames,
                            std::vector<xercesc::DOMElement*>& rResults,
                            unsigned depth)
{
    for (xercesc::DOMNode* p_node = pContextElement->getFirstChild();
         p_node != NULL;
         p_node = p_node->getNextSibling())
    {
        if (p_node->getNodeType() == xercesc::DOMNode::ELEMENT_NODE &&
            X2C(p_node->getLocalName()) == rNames[depth])
        {
            xercesc::DOMElement* p_child_elt = static_cast<xercesc::DOMElement*>(p_node);
            if (depth == rNames.size() - 1)
            {
                rResults.push_back(p_child_elt);
            }
            else
            {
                FindElements(p_child_elt, rNames, rResults, depth+1);
            }
        }
    }
}

std::vector<xercesc::DOMElement*> XmlTools::FindElements(const xercesc::DOMElement* pContextElement,
                                                         const std::string& rPath)
{
    std::vector<xercesc::DOMElement*> results;
    std::vector<std::string> path;
    size_t start_pos = 0;
    size_t slash_pos = 0;
    while (slash_pos != std::string::npos)
    {
        slash_pos = rPath.find('/', start_pos);
        if (slash_pos == std::string::npos)
        {
            path.push_back(rPath.substr(start_pos));
        }
        else
        {
            path.push_back(rPath.substr(start_pos, slash_pos-start_pos));
        }
        start_pos = slash_pos + 1;
    }
    FindElements(pContextElement, path, results);
    return results;
}


void XmlTools::WrapContentInElement(xercesc::DOMDocument* pDocument,
                                    xercesc::DOMElement* pElement,
                                    const XMLCh* pNewElementLocalName)
{
    const XMLCh* p_namespace_uri = pElement->getNamespaceURI();
    const XMLCh* p_prefix = pElement->getPrefix();
    const XMLCh* p_qualified_name;
    if (p_prefix)
    {
// LCOV_EXCL_START
        // We can't actually cover this code, since versions of the parameters file which need this
        // transform didn't use a namespace, so can't have a namespace prefix!
        xercesc::QName qname(p_prefix, pNewElementLocalName, 0);
        p_qualified_name = qname.getRawName();
// LCOV_EXCL_STOP
    }
    else
    {
        p_qualified_name = pNewElementLocalName;
    }
    xercesc::DOMElement* p_wrapper_elt = pDocument->createElementNS(p_namespace_uri, p_qualified_name);
    // Move all child nodes of pElement to be children of p_wrapper_elt
    xercesc::DOMNodeList* p_children = pElement->getChildNodes();
    for (unsigned i=0; i<p_children->getLength(); i++)
    {
        xercesc::DOMNode* p_child = pElement->removeChild(p_children->item(i));
        p_wrapper_elt->appendChild(p_child);
    }
    // Add the wrapper as the sole child of pElement
    pElement->appendChild(p_wrapper_elt);
}


std::string XmlTools::EscapeSpaces(const std::string& rPath)
{
    std::string escaped_path;
    for (std::string::const_iterator it = rPath.begin(); it != rPath.end(); ++it)
    {
        if (*it == ' ')
        {
            escaped_path += "%20";
        }
        else
        {
            escaped_path += *it;
        }
    }
    return escaped_path;
}
