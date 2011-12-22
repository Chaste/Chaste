This folder contains saved Eclipse preferences that, at least on Ganymede,
can be used to set up Eclipse with normal Chaste developer settings.  To
use them, go to File > Import... and then select:

1) General > Preferences, and use the preferences.epf file.
2) Run/Debug > Launch Configurations, and use the whole eclipse_settings folder.

Among other things, the preferences contain Chaste code styles and templates
for header files, source files, and tests.  The latter must be selected explicitly
when creating a new file.

We haven't yet figured out how to include certain project-specific settings, found
by right-clicking on the Chaste project in the Navigator pane and selecting
Properties.

1) On the C/C++ General tab, tick "Enable project specific settings" and select
   'Doxygen' as the documentation tool.
