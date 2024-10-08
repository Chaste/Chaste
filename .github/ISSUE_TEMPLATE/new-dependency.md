---
name: New Dependency
about: Support and test against a new version of one of Chaste's dependencies
title: 'New dependency version: [library] [version]'
labels: 'component: infrastructure, enhancement'
assignees: ''

---

Library that has a new version: [e.g. PETSc / boost]

Link to release notes: [add here]

**TODO:**
- [ ]  Build a docker image containing the new dependency version - see the [Portability Testing Guide](https://github.com/Chaste/developer_wiki/wiki/Portability-Testing-Guide) for how to do this.
- [ ]  Create a branch off `develop` - you could name it, for example, `support-petsc-322`.
- [ ]  Add the new docker image to the [portability workflow](https://github.com/Chaste/Chaste/blob/develop/.github/workflows/portability.yml) on the new branch.
- [ ]  Create a PR with the new branch and see if it builds and passes all tests.
- [ ]  Do any necessary fixes.
- [ ]  Add permanently to testing rotation i.e. either replace an old docker image in the [portability workflow](https://github.com/Chaste/Chaste/blob/develop/.github/workflows/portability.yml) with the new image, or append the new image to the existing list of portability images. For example, if there's an old dependency version that is being retired, the new docker image can probably replace the old one.
- [ ]  Update the portability versions in the [Portability Testing Guide](https://github.com/Chaste/developer_wiki/wiki/Portability-Testing-Guide).
- [ ]  Update the [supported dependency versions](https://chaste.github.io/docs/installguides/dependency-versions/) page.
