"""Copyright (c) 2005-2023, University of Oxford.
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
"""

import os
import re
import sys
import urllib
import urllib2
# For later? from urllib.request import urlopen

# Load our helper modules
our_path = os.path.dirname(os.path.realpath(__file__))
sys.path[0:0] = [our_path]
import TracWebUtils

#
# Configuration (extra settings at end of file)
#
# Where Trac is located
SERVER = "https://chaste.cs.ox.ac.uk"
BASE_URL = "/trac/"
# How to recognise page content by scraping the HTML
CONTENT_START = '<div id="content" class="wiki">'
CONTENT_END = '<div class="trac-modifiedby">'
# Header & footer for static pages
HEAD = """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" lang="en" xml:lang="en">
  <head>
    <title>%(title)s</title>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
    <link rel="stylesheet" href="/trac-static/css/trac.css" type="text/css" />
    <link rel="stylesheet" href="/trac-static/css/wiki.css" type="text/css" />
    <link rel="stylesheet" type="text/css" href="/trac-static/css/site.css" />
  </head>
  <body>
    <div id="banner">
      <div id="header">
        <p>
          <a id="logo" href="http://www.cs.ox.ac.uk/chaste"><img src="/logos/chaste-266x60.jpg" alt="Chaste logo" height="60" width="266" /></a>
          <em>Documentation for <a href="%(output_base_url)s">%(release_name)s</a>.</em>
        </p>
      </div>
    </div>
    <p>&nbsp;</p>
"""
FOOT = """
    </div>
  </body>
</html>
"""

# Log-in to Trac
TracWebUtils.TRAC_URL = SERVER + BASE_URL
#TracWebUtils.TracLogin()


class WikiConverter(object):
    """Takes care of converting a collection of Trac wiki pages to static versions."""
    def __init__(self, sections, releaseName, outputFolder, outputBaseUrl, indexPage=None):
        # Store args
        self.sections = sections
        self.releaseName = releaseName
        self.outputFolder = outputFolder
        if outputBaseUrl and outputBaseUrl[-1] != '/':
            outputBaseUrl += '/'
        self.outputBaseUrl = outputBaseUrl
        self.indexPage = indexPage
        # Queue of wiki pages to process
        self.pages = sections[:]
        # Which pages have already been done
        self.pagesDone = set()
        # Attachments to download
        self.attachments = set()
        # What links to other processable pages look like
        section_part = '|'.join(sections)
        common_re = r'(?P<attr>href|src)="%%s(?P<section>%s)(?P<page>[^#?&]*?)(?P<fragment>#.*?)?"' % section_part
        self.wikiHrefRegexp = re.compile(common_re % (BASE_URL+'wiki/'))
        self.attachmentHrefRegexp = re.compile(common_re % (BASE_URL+'(?:raw-)?attachment/wiki/'))
        self.doxygenHrefRegexp = re.compile(r'href="(?:%s)?/(?:public-)?docs/?(.*?)"' % SERVER)
        self.releaseDoxRegexp = re.compile(r'(<span style.*REL_DOXY.*</span>)')
        # Some bits of content we don't want
        strip = [r'\xc2', # Added in external web links!
                 r' - added by <em>.*?timeline.*?</a> ago.']
        self.stripRegexp = re.compile('|'.join(strip), re.DOTALL)
        self.replacements = {'\xe2\x80\x94': '&mdash;',
                             '\xa0': '&nbsp;',
                             '\xe2\x80\xa6': '&hellip;'}

    def _debug(self, *args):
        print(' '.join(map(str, args)))

    def PageUrl(self, wikiName):
        """Get the full URL to a wiki page."""
        return SERVER + BASE_URL + 'wiki/' + wikiName

    def AttachmentGetUrl(self, attachment):
        """Get the download URL for an attachment."""
        return SERVER + BASE_URL + 'raw-attachment/wiki/' + attachment

    def DocsUrl(self, docsPage):
        """Get the full URL to a static doxygen page."""
        return SERVER + '/chaste/docs/' + self.outputFolder + '/' + docsPage

    def GenerateHtml(self, wikiName, title):
        """Download a wiki page and return a static page with its interesting content."""
        current_page = urllib2.urlopen(self.PageUrl(wikiName))
        copy = False
        subst = {'title': title,
                 'output_base_url': self.outputBaseUrl,
                 'release_name': self.releaseName}
        output = [HEAD % subst]
        for line in current_page:
            if line.strip() == CONTENT_END:
                # Finished
                break
            if line.strip() == CONTENT_START:
                copy = True
            if copy:
                output.append(self.ConvertLinks(line))
        output.append(FOOT)
        output = self.stripRegexp.sub('', ''.join(output))
        for old, new in self.replacements.iteritems():
            output = output.replace(old, new)
        return output

    def _ReplFunction(self, match):
        """Used by ConvertLinks."""
        subst = match.groupdict('')
        subst['base'] = self.outputBaseUrl
        if match.re is self.attachmentHrefRegexp:
            subst['ext'] = ''
        else:
            subst['ext'] = '.html'
        return '%(attr)s="%(base)s%(section)s%(page)s%(ext)s%(fragment)s"' % subst

    def ConvertLinks(self, line):
        """Convert links to other processable wiki pages to refer to the static version.

        Also updates the queue of processable pages, and converts doxygen links.
        """
        # Wiki links
        for groups in self.wikiHrefRegexp.findall(line):
            page = ''.join(groups[1:3]) # Don't include the fragment!
            self._debug('Found page link to', page)
            self.pages.append(page)
        line = self.wikiHrefRegexp.sub(self._ReplFunction, line)
        # Attachment links
        for groups in self.attachmentHrefRegexp.findall(line):
            attachment = ''.join(groups[1:])
            self._debug('Found attachment link to', attachment)
            self.attachments.add(attachment)
        line = self.attachmentHrefRegexp.sub(self._ReplFunction, line)
        # Documentation links
        if line.find('<span style="display:none"') > 0:
            link = 'including <a href="%s"><b>documentation for this release</b></a>' % self.DocsUrl('')
            line = self.releaseDoxRegexp.sub(link, line)
        else:
            doxy_replace = 'href="%s"' % self.DocsUrl(r'\1')
            line = self.doxygenHrefRegexp.sub(doxy_replace, line)
        return line

    def EnsureFolderExists(self, page):
        """Ensure all the folders needed for the static copy of page exist."""
        parts = page.split('/')
        folders = parts[:-1]
        path = os.sep.join([self.outputFolder] + folders)
        if not os.path.exists(path):
            os.makedirs(path)
        return path, parts[-1]

    def ProcessPages(self):
        """Go through our queue, processing wiki pages."""
        if not self.indexPage:
            self.CreateIndexPage()
        while self.pages:
            page = self.pages.pop()
            if not page in self.pagesDone:
                self._debug('Processing page', page)
                # Ensure there's a folder to put it in
                folder_path, leaf_name = self.EnsureFolderExists(page)
                # Generate the HTML
                html = self.GenerateHtml(page, page + ' - Chaste')
                # Save the file
                fp = open(os.path.join(folder_path, leaf_name + '.html'), 'w')
                fp.write(html)
                fp.close()
                self.pagesDone.add(page)
        for attachment in self.attachments:
            self._debug('Downloading attachment', attachment)
            folder_path, leaf_name = self.EnsureFolderExists(attachment)
            local_file = open(os.path.join(folder_path, leaf_name), 'wb')
            for line in urllib2.urlopen(self.AttachmentGetUrl(attachment)):
                local_file.write(line)
            local_file.close()
        if self.indexPage:
            self.LinkIndexPage(self.indexPage)

    def LinkIndexPage(self, indexPage):
        """Create index.html as a symlink to the static version of the given wiki page."""
        self.EnsureFolderExists('')
        page_file = indexPage + '.html'
        if not os.path.exists(os.path.join(self.outputFolder, page_file)):
            raise ValueError("Unable to use %s as index since it doesn't exist" % indexPage)
        os.symlink(page_file, os.path.join(self.outputFolder, 'index.html'))

    def CreateIndexPage(self):
        """Create an index page linking to all the static pages created."""
        self.EnsureFolderExists('')
        fp = open(os.path.join(self.outputFolder, 'index.html'), 'w')
        subst = {'title': 'Chaste Tutorials and Manuals',
                 'output_base_url': self.outputBaseUrl,
                 'release_name': self.releaseName}
        output = [HEAD % subst, CONTENT_START]
        output.append("""
      <h1>%(title)s</h1>
      <ul>
""" % subst)
        for section in self.sections:
            output.append('        <li><a href="%s.html">%s</a></li>' % (section, section))
        output.extend(['      </ul>\n', FOOT])
        fp.write(''.join(output))
        fp.close()

if __name__ == '__main__':
    if len(sys.argv) > 1:
        sections = sys.argv[1:]
    else:
        sections = ['UserTutorials', 'ChasteGuides', 'InstallGuides', 'GettingStarted', 'ReleaseNotes']
    
    # Edit the output folder (which will dictate the URL of these files on the webserver) and the releaseName which will appear in the files.
    output_folder = 'release_2021.1' #Edit (and below)
    converter = WikiConverter(sections,
                              releaseName='Release 2021.1', #Edit (and above)
                              outputFolder=output_folder,
                              outputBaseUrl='/chaste/tutorials/'+output_folder,
                              indexPage='GettingStarted')
    converter.ProcessPages()
    print("You now need to copy the generated files to the correct location on the web server, using something like:")
    print("scp -r %s chaste@chaste.cs.ox.ac.uk:/var/www/html/chaste/tutorials/" % output_folder)
    print("I'm going to try doing this for you; you will be prompted three times for a password.")
    os.system('scp -r '+ output_folder + ' chaste@chaste.cs.ox.ac.uk:')
    os.system('ssh -t chaste@chaste.cs.ox.ac.uk "sudo mv '+ output_folder +' /var/www/html/chaste/tutorials/"')
