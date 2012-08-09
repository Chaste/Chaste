Using CellML with (Cardiac) Chaste
==================================

There is a companion tool to Chaste, named PyCml, which can generate
Chaste-compatible C++ code from cardiac ionic cell models described in
CellML.  This allows Chaste to make use of any such model.  It also applies
some optimisations to improve the speed of simulations.

Since PyCml generates C++ source code, the source code release of Chaste is
required to make use of CellML models.

Full documentation can be found on the developers' wiki at
https://chaste.cs.ox.ac.uk/cgi-bin/trac.cgi/wiki/ChasteGuides/CodeGenerationFromCellML
For a guest login, use the username "anonymous", and your email address as
the password.

Static versions of the guide accompanying Chaste releases are also available.  See
https://chaste.cs.ox.ac.uk/chaste/tutorials/release_3.1/ChasteGuides/CodeGenerationFromCellML.html
for release 3.1.
