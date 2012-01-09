
"""Copyright (C) University of Oxford, 2005-2012

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
