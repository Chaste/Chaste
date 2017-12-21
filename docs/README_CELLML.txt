Using CellML with (Cardiac) Chaste
==================================

There is a companion tool to Chaste, named PyCml, which can generate
Chaste-compatible C++ code from cardiac ionic cell models described in
CellML.  This allows Chaste to make use of any such model.  It also applies
some optimisations to improve the speed of simulations.

Since PyCml generates C++ source code, the source code release of Chaste is
required to make use of CellML models.

Full documentation can be found on the Chaste wiki at
https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides/CodeGenerationFromCellML

Static versions of the guide accompanying Chaste releases are also available.  See
https://chaste.cs.ox.ac.uk/chaste/tutorials/release_2017.1/ChasteGuides/CodeGenerationFromCellML.html
for release 2017.1.
