#!/bin/bash

# Add local chaste user
# Either use the USER_ID if passed in at runtime or fallback to 1000

USER_ID=${USER_ID:-1000}

echo "Starting with UID : $USER_ID"
useradd --shell /bin/bash -u $USER_ID -o -c "" --create-home chaste

# The -l makes bash behave like a login shell, so the USER & HOME
# environment variables will be set.  The -c at the end of the ENTRYPOINT
# means we can do things like "docker run -it ... scons" to build Chaste.
exec sudo -u chaste /bin/bash -l -c "$@"
