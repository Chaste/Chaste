A couple of notes...

You might expect these files to be in 'io' but they depend on mesh, so they can't be in there.

Note: there are two more HDF5 converters that also depend on heart as well.

Hdf5ToMeshalyzerConverter
Hdf5ToCmguiConverter

these are in heart/src/postprocessing as they currently use methods from HeartConfig. 

But eventually they should be migrated to here, and HeartConfig arguments passed into constructors, this started for 
Hdf5ToMeshalyzerConverter as part of #1660.