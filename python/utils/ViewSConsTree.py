
"""Copyright (c) 2005-2019, University of Oxford.
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

# This is a quick program that parses the output of a
#      scons --debug=tree
# run and displays it in a more friendly format, as a
# tree widget. It uses wxPython, so you'll need to have that
# installed. It's my first wxPython program, so feel free to point out
# anything I've done wrong.

# From http://www.scons.org/wiki/SconsTreeView

# Usage:
#  Run 'scons --debug=tree > tree.txt'
#  Run 'python python/ViewSConsTree.py'
#  Open the tree.txt file, and explore the dependency graph.

from wxPython.wx import *
import os
import sys

ID_OPEN=102

class SconsTreeView(wxApp):
  def OnInit(self):
    frame = wxFrame(NULL, -1, 'Scons tree view')
    self.frame=frame
    frame.Show(true)
    self.SetTopWindow(frame)

    filemenu = wxMenu()
    filemenu.Append(ID_OPEN, '&Open')
    menuBar = wxMenuBar()
    menuBar.Append(filemenu, '&File')
    frame.SetMenuBar(menuBar)
    EVT_MENU(self, ID_OPEN, self.OnOpen)

    self.tree = wxTreeCtrl(frame)
    sizer = wxBoxSizer(wxHORIZONTAL)
    sizer.Add(self.tree, 1, wxEXPAND)
    frame.SetSizer(sizer)
    frame.SetAutoLayout(1)
    sizer.Fit(frame)
    if (len(sys.argv) > 1):
      f = open(sys.argv[1], 'r')
      print f
      self.Parse(f)
    return true

  def OnOpen(self, evt):
    dlg = wxFileDialog(self.frame, "Choose a file", "", "", "*.*", wxOPEN)
    if dlg.ShowModal() == wxID_OK:
      filename = dlg.GetFilename()
      dirname = dlg.GetDirectory()
      f = open(os.path.join(dirname, filename), 'r')
      self.Parse(f)
      f.close()
    dlg.Destroy()

  def Parse(self, f):
    t = self.tree
    t.DeleteAllItems()
    for line in f:
      if line.startswith('+') and line.find('-') != -1:
        lastdepth = line.index('-')/2
        nodes = [ t.AddRoot('.') ]
        break

    lastnode = nodes[-1]
    for line in f:
      if line.find('-') == -1:
        break
      depth=line.index('-')/2
      content = line[line.index('-')+1:]
      if depth > lastdepth:
        nodes.append(lastnode)
        lastnode = t.AppendItem(nodes[-1], content)
      elif depth < lastdepth:
        for i in range(depth, lastdepth):
          nodes.pop()
        lastnode = t.AppendItem(nodes[-1], content)
      else:
        lastnode = t.AppendItem(nodes[-1], content)
      lastdepth = depth

app=SconsTreeView(0)
app.MainLoop()
