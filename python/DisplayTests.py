
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


# Chaste tests display script.

# This module contains most of the functionality, and is loaded from the
# repository by a wrapper script.

import glob
import os
import operator
import time
import itertools
import re

# Compatibility with Python 2.3
try:
    set = set
except NameError:
    import sets
    set = sets.Set

_standalone = False
_db_module = None

#####################################################################
##                          Webpages                               ##
#####################################################################

def index(req):
    """The main test display page.
    
    This displays a summary of the most recent tests.
    """
    output = [_header()]
    output.append("""\
    <h1>Chaste Tests</h1>
    <p>
    This is the funky interface to Chaste's testing suite.
    </p>
    <ul>""")
    tests_types = os.listdir(_tests_dir)
    tests_types = filter(lambda s: s[0] != '.', tests_types)
    tests_types.sort()
    main_types = filter(lambda s: '-' not in s, tests_types)
    branch_types = filter(lambda s: '-' in s, tests_types)
    for tests_type in main_types:
        if tests_type.endswith('_old'):
            text = 'Old %s builds.' % tests_type[:-4]
        else:
            text = 'Recent %s builds.' % tests_type
        output.append('\n    <li><a href="%s/recent?type=%s">%s</a></li>' % (_our_url, tests_type, text))
        if tests_type == 'nightly':
            # Special case links to show only or all but project builds
            output.append('\n    <ul>')
            output.append('\n        <li><a href="%s/recent?type=nightly&targets=!*projects*">Only core builds</a></li>' % _our_url)
            output.append('\n        <li><a href="%s/recent?type=nightly&targets=*projects*">Only project builds</a></li>' % _our_url)
            output.append('\n    </ul>')
    output.append("""</ul>
    <p>Branch/project builds: (<a style='text-decoration: underline; color: blue;' onclick="toggle_visibility('branch-list');">toggle visibility</a>)</p>
    <ul id='branch-list' style='display:none;'>""")
    for tests_type in branch_types:
        output.append('\n    <li><a href="%s/recent?type=%s">Recent %s builds.</a></li>' % (_our_url, tests_type, tests_type))
    output.append("""</ul>
    
    <p>
    <a href="%s/profileHistory?n=100">Run time variation of profiled tests.</a>
    <br />
    <a href="%s/profileHistory?n=100&buildTypes=IntelProductionParallel4_onlytests_Weekly&buildTypes=IntelProduction_onlytests_Weekly">Run time variation of weekly tests.</a>
    </p>
    
    <h2>Latest continuous build</h2>
""" % (_our_url, _our_url))

    # Look for the latest revision present.
    type = 'continuous'
    revision_and_timestamps = os.listdir(os.path.join(_tests_dir, type))
    latest_revision = str(max(itertools.imap(lambda rts: int(rts[:rts.find('~')]) if '~' in rts else int(rts), revision_and_timestamps)))
    # Display summary of each build of this revision
    test_set_dirs = _getResultsParentDirs(type, latest_revision, timestamp=None)
    builds = []
    for test_set_dir in test_set_dirs:
        revision_and_timestamp = os.path.basename(test_set_dir)
        builds.extend(map(lambda machine_and_buildtype: (revision_and_timestamp, machine_and_buildtype), os.listdir(test_set_dir)))
    if len(builds) < 1:
        output.append(_error('No test set found for revision '+latest_revision+'. Probably the build is still in progress.'))
        output.append('<p><a href="/out/latest">Latest build log.</a></p>')
    else:
        for build in builds:
            revision, timestamp = _extractTildeSeparatedPair(build[0])
            machine, build_type = _extractDotSeparatedPair(build[1])
            output.append(_summary(req, type, latest_revision, machine, build_type, timestamp))

    output.append(_footer())
    return ''.join(output)


def testsuite(req, type, revision, machine, buildType, testsuite, status, runtime, timestamp=None):
    """Display the results for the given testsuite, by passing the file back to the user."""
    req.content_type = 'text/html'
    req.write(_header(), 0)
    test_set_dir = _testResultsDir(type, revision, machine, buildType, timestamp)
    buildTypesModule = _importBuildTypesModule(revision)
    build = _getBuildObject(buildTypesModule, buildType)
    for suite_name in _test_suite_name_aliases(testsuite):
        if _isWindows(build):
            ctest_results = glob.glob(os.path.join(test_set_dir, '*TestOutputs_*.txt*'))
            testsuite_file = ''.join(ctest_results) #TODO: Hack! (#2016)
        else:
            testsuite_file = build.ResultsFileName(test_set_dir, suite_name, status, float(runtime))
        if os.path.isfile(testsuite_file):
            req.write('\n<pre>\n', 0)
            fp = open(testsuite_file)
            for line in fp:
                req.write(line.replace('&', '&amp;').replace('<', '&lt;'))
            fp.close()
            req.write('\n</pre>\n', 0)
            break
    else:
        req.write(_error('The requested test suite was not found.'))
    req.write(_footer())


def graph(req, type, revision, machine, buildType, graphName, timestamp=None):
    """Send the given graph file (a .gif) back to the user."""
    test_set_dir = _testResultsDir(type, revision, machine, buildType, timestamp)
    graph_name = os.path.basename(graphName) # Foil hackers
    for graph_name in _test_suite_name_aliases(graph_name):
        graph_path = os.path.join(test_set_dir, graph_name)
        if os.path.isfile(graph_path):
            req.content_type = 'image/gif'
            req.sendfile(graph_path)
            break
    else:
        req.content_type = 'text/html'
        req.write(_header(), 0)
        req.write('\n<pre>%s\n%s</pre>\n' % (graphName, _test_suite_name_aliases(os.path.basename(graphName))))
        req.write(_error('The requested graph file was not found.'))
        req.write(_footer())


def recent(req, type='', start=0, n_per_page=30, **filters):
    "User-facing page. Content is generated by _recent."
    title = 'Recent '+type+' builds'
    page_body = """\
    <h1>%s</h1>
%s
""" % (title, _recent(req, type, int(start), int(n_per_page), **filters))
    return _header(title) + page_body + _footer()

def _recent(req, type='', start=0, n_per_page=30, **filters):
    """Display brief summaries of recent builds of the given type.
    
    Returns a string representing part of a webpage.
    """
    if not type:
        return _error('No type of test to summarise specified.')

    dir = os.path.join(_tests_dir, type)
    if not os.path.isdir(dir):
        return _error(type+' is not a valid type of test.')

    # What keys may be used to filter the builds list
    poss_filters = ['build_type', 'machine', 'targets']

    if _db_module and not req.form.getfirst('nocache', False):
        # Get information from the database
        db = _db_module.TestResultsDatabase(type, verbose=False)
        total_num_of_builds = db.CountResults()
        params = []
        where = []
        if filters:
            for filter_name, filter_value in filters.items():
                if filter_name in poss_filters:
                    if len(filter_value) > 0 and filter_value[0] == '!':
                        filter_value = filter_value[1:]
                        negated = True
                    else:
                        negated = False
                    op = [['=', '!='], ['GLOB', 'NOT GLOB']]['*' in filter_value][negated]
                    where.append(filter_name + ' ' + op + ' ?')
                    params.append(filter_value)
        if where:
            where = ' where ' + ' and '.join(where)
        else:
            where = ''
        params = tuple(params + [n_per_page, start])
        cur = db.conn.execute('select * from summary%s'
                              ' order by finished desc, revision desc, build_type, machine'
                              ' limit ? offset ?' % where, params)
        def gen_row(cur=cur):
            for row in cur:
                yield (row['revision'], row['machine'], row['build_type'], row['targets'],
                       db.FormatTimestamp(row['finished']), row['status'], row['colour'])
    else:
        # Parse the directory structure within dir into a list of builds
        builds = []
        for revision_and_timestamp in os.listdir(dir):
            for machine_and_build_type in os.listdir(os.path.join(dir, revision_and_timestamp)):
                revision, timestamp = _extractTildeSeparatedPair(revision_and_timestamp)
                if timestamp is None:
                    st = os.stat(os.path.join(dir, revision_and_timestamp, machine_and_build_type))
                    mod_time = st.st_mtime
                else:
                    mod_time = float(timestamp)
                machine, build_type = _extractDotSeparatedPair(machine_and_build_type)
                builds.append([mod_time, revision, build_type, machine])
        # Sort the list to be most recent first
        builds.sort()
        builds.reverse()
        # Just show a subset
        total_num_of_builds = len(builds)
        builds = builds[start:start+n_per_page] # About a screenful
        def gen_row(builds=builds):
            old_revision = -1
            for finished, revision, build_type, machine in builds:
                if revision != old_revision:
                    build_types_module = _importBuildTypesModule(revision)
                    old_revision = revision
                build = _getBuildObject(build_types_module, build_type)
                test_set_dir = _testResultsDir(type, revision, machine, build_type, finished)
                targets, overall_status, colour = _getTestSummary(test_set_dir, build)
                yield (revision, machine, build_type, targets, finished, overall_status, colour)

    def add_nav_links():
        if start > 0:
            output.append(_linkRecent('First page', type, start=0, **filters) + " ")
            output.append(_linkRecent('Previous page', type, start=start-n_per_page, **filters) + " ")
        if filters:
            output.append(_linkRecent('Remove filters', type, start=0) + " ")
        if total_num_of_builds > start+n_per_page:
            output.append(_linkRecent('Next page', type, start=start+n_per_page, **filters) + " ")
            output.append(_linkRecent('Last page', type, start=total_num_of_builds-n_per_page, **filters) + " ")

    output = []
    add_nav_links()
    output.append("""\
    <table border="1">
        <tr>
          <th>Date</th>
          <th>Revision</th>
          <th>Build Type</th>
          <th>Machine</th>
          <th>Targets</th>
          <th>Status</th>
        </tr>
""")

    bgcols = ["white", "#eedd82"]
    bgcol_index = 0
    old_revision = -1
    for revision, machine, build_type, targets, finished, overall_status, colour in gen_row():
        if type == 'nightly':
            date = time.strftime('%d/%m/%Y', time.localtime(float(finished)))
        else:
            date = time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(float(finished)))
        if revision != old_revision:
            bgcol_index = 1 - bgcol_index
            old_revision = revision
        subs = {'bgcol': bgcols[bgcol_index], 'status_col': colour,
                'date': date, 'machine': machine, 'targets': targets or 'default',
                'rev': _linkRevision(revision, changes=True),
                'build_type': _linkBuildType(build_type, revision, wrappableText=True),
                'status': _linkSummary(overall_status, type, revision, machine, build_type, finished)}
        for poss_filter in poss_filters:
            if poss_filter not in filters or '*' in filters[poss_filter]:
                new_filters = filters.copy()
                new_filters[poss_filter] = locals()[poss_filter]
                subs[poss_filter] += ' [' + _linkRecent('filter', type, 0, **new_filters) + ']'
        output.append("""\
        <tr>
          <td style='background-color: %(bgcol)s;'>%(date)s</td>
          <td style='background-color: %(bgcol)s;'>%(rev)s</td>
          <td style='background-color: %(bgcol)s;'>%(build_type)s</td>
          <td style='background-color: %(bgcol)s;'>%(machine)s</td>
          <td style='background-color: %(bgcol)s;'>%(targets)s</td>
          <td style='background-color: %(status_col)s;'>%(status)s</td>
        </tr>
""" % subs)
    output.append("  </table>\n")
    add_nav_links()
    if type == 'nightly':
        output.append('<p><a href="/out/latest-nightly">Latest nightly build log.</a></p>')
    elif type == 'continuous':
        output.append('<p><a href="/out/latest">Latest continuous build log.</a></p>')
    
    return ''.join(output)



def summary(req, type, revision, machine, buildType, timestamp=None):
    "User-facing page. Content is generated by _summary."
    page_body = """\
    <h1>Build Summary</h1>
""" +  _summary(req, type, revision, machine, buildType, timestamp)
    return _header('Test Summary') + page_body + _footer()
    
def _summary(req, type, revision, machine=None, buildType=None, timestamp=None):
    """Display a summary of a build.
    
    Returns a string representing part of a webpage.
    """
    output = []
    if not (type and revision):
        return _error('No test set to summarise specified.')
    if not (machine and buildType):
        return _error('No test set to summarise specified.')
    # Find the directory with appropriate test results
    if type == 'standalone':
        test_set_dir = _dir
    else:
        test_set_dir = _testResultsDir(type, revision, machine, buildType, timestamp)
    
    # Now test_set_dir should be the directory containing the test results
    # to summarise. Extract summary info from the filenames.
    if type == 'standalone':
        build = _build
    else:
        buildTypesModule = _importBuildTypesModule(revision)
        build = _getBuildObject(buildTypesModule, buildType)
    testsuite_status, overall_status, colour, runtime, graphs = _getTestStatus(test_set_dir, build)
    targets = _getBuildTargets(test_set_dir)
    # Store overall status for the standalone script case
    if _standalone:
        global _overall_status
        _overall_status = overall_status

    # Get the timestamp on the directory, if not specified
    if timestamp:
        mod_time = float(timestamp)
    else:
        st = os.stat(test_set_dir)
        mod_time = st.st_mtime
    date = time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(mod_time))

    # Work out the URL of the build log file
    # i.e. find out where test_set_dir/build.log points to
    build_log = "</p>"
    if not _standalone:
        logpath = os.path.realpath(test_set_dir + "/build.log")
        # Remove the '/var/www' part
        logurl = logpath[8:]
        if logurl:
            build_log = "Build log: <a href=\"%s\">%s</a></p>" % (logurl, logurl)
            timings = _parseBuildTimings(logpath).items()
            timings.sort(key=operator.itemgetter(1)) # Sort by time
            timings.reverse()
            build_log += "\nTimings (for entire build log): <table><tr><th>Activity</th><th>Time (minutes)</th></tr>\n"
            build_log += "\n".join(
                map(lambda row: '<tr><td>%s</td><td>%d:%f</td></tr>' % (row[0], row[1]//60, row[1]%60),
                    timings))
            build_log += "\n</table>\n"
    
    # Produce output HTML
    if graphs:
        extra_cols = '<th>Extras</th>'
    else:
        extra_cols = ''
    output.append("""\
    <p>
    Revision: %s<br />
    Date and time: %s<br />
    Overall status: %s<br />
    Build type: %s<br />
    Machine: %s<br />
    Targets: %s<br />
    %s
    <table border="1">
        <tr>
          <th>Test Suite</th>
          <th>Status</th>
          <th>Run Time</th>
          %s
        </tr>
""" % (_linkRevision(revision, changes=True), date, _colourText(overall_status, colour),
       _linkBuildType(buildType, revision), machine, targets, build_log, extra_cols))
    
    # Display the status of each test suite, in alphabetical order
    testsuites = testsuite_status.keys()
    testsuites.sort()
    bgcols = ["white", "#eedd82"]
    bgcol_index = 0
    for testsuite in testsuites:
        bgcol_index = 1 - bgcol_index
        subs = {'bgcol': bgcols[bgcol_index],
                'status_col': _statusColour(testsuite_status[testsuite], build),
                'testsuite': testsuite,
                'status': _linkTestSuite(type, revision, machine, buildType, timestamp,
                                         testsuite, testsuite_status[testsuite],
                                         runtime[testsuite], build),
                'runtime': _formatRunTime(runtime[testsuite]),
                'graph': _linkGraph(type, revision, machine, buildType, timestamp,
                                    graphs.get(testsuite,''))}
        output.append("""\
        <tr>
          <td style='background-color: %(bgcol)s;'>%(testsuite)s</td>
          <td style='background-color: %(status_col)s;'>%(status)s</td>
          <td style='background-color: %(bgcol)s;'>%(runtime)s</td>
          %(graph)s
        </tr>
""" % subs)

    output.append("  </table>\n")
    
    return ''.join(output)


def buildType(req, buildType, revision=None):
    """Display information on the compiler settings, etc. used to build a set of tests.
    buildType is the user-friendly name describing these settings, such as can be passed to scons build=buildType.
    revision is the code revision of the set of tests, in case the definition of buildType has changed since.
    """
    if revision is None:
        rev_text = ' at the latest revision'
    else:
        rev_text = ' at revision %s' % revision
    BuildTypes = _importBuildTypesModule(revision)
    build = _getBuildObject(BuildTypes, buildType)
    test_packs = ', '.join(build.TestPacks())
    # How test suites are run
    testsuite_exe = "testsuite.exe"
    testsuite_cmd = build.GetTestRunnerCommand(testsuite_exe)
    page_body = """\
    <h1>Explanation of build type '%s'%s</h1>
    <p>
    C++ compiler 'brand': %s<br />
    C++ extra compile flags: %s<br />
    Extra linker flags: %s<br />
    Test packs run: %s<br />
    Command to run '%s': %s<br />
    </p>
""" % (buildType, rev_text, build.CompilerType(),
       _getCcFlags(build), build.LinkFlags(),
       test_packs, testsuite_exe, testsuite_cmd)
    if buildType.startswith('acceptance') or buildType.startswith('longacceptance') :
        page_body += """
    <p>This build was actually used just to run the acceptance tests.</p>
"""
    elif buildType.startswith('_broke'):
        page_body += """
    <p>This build didn't actually succeed; check the build log for details.</p>
"""
    prefs = build.GetPreferedVersions()
    if prefs:
        page_body += """
    <p>Library/tool versions requested:</p>
    <ul>
"""
        for lib, version in prefs.iteritems():
            page_body += "    <li>%s: %s</li>\n" % (lib, version)
    page_body += "\n  </ul>\n"
    return _header() + page_body + _footer()


class FakeBuildType(object):
    """Fake build type.  Provides just enough for the profile history pages."""
    def GetInfoFromResultsFileName(self, leafname):
        """Extract the metadata held within the name of a results file.
        
        This returns a dictionary, with keys 'testsuite', 'status' and 'runtime'.
        testsuite is the name of the test suite.
        status is the encoded status string.
        runtime is the run time for the test suite in seconds.
        """
        # Components are separated by '.'
        i2 = leafname.rfind('.')
        i1 = leafname.rfind('.', 0, i2)
        if i1 == -1:
            # No runtime info available
            runtime = -1
            i1, i2 = i2, len(leafname)
        else:
            runtime = int(leafname[i2+1:])
        return {'testsuite': leafname[:i1],
                'status': leafname[i1+1:i2],
                'runtime': runtime}
    
    def ParseGraphFilename(self, leafname):
        """Return the test suite a graph file is associated with.
        
        Removes the string 'Runner.gif' from the end of the filename.
        """
        return leafname[:-10]
    
    def StatusColour(self, status):
        """
        Return a colour string indicating whether the given status string
        represents a 'successful' test suite under this build type.
        """
        # By default, 'OK' is ok and anything else isn't.
        if status == 'OK':
            return 'green'
        else:
            return 'red'


def profileHistory(req, n=20, buildTypes=None):
    """Show runtimes for the last n profile builds."""
    page_body = """
    <h1>Profile History</h1>
""" +  _profileHistory(req, int(n), buildTypes)
    return _header('Profile History') + page_body + _footer()

def _profileHistory(req, n=20, buildTypes=None):
    """Show runtimes for the last n profile builds."""
    if not _db_module:
        return _error('No database connection present; unable to produce performance plots')
    tests_dir = os.path.join(_tests_dir, 'nightly')
    if not os.path.isdir(tests_dir):
        return _error('No nightly tests found')
    output = []

    # These are the build types representing profile builds
    if not buildTypes:
        # Wildcard match on anything containing Profile
        buildTypes = ['*Profile*']
    if not isinstance(buildTypes, list):
        buildTypes = [buildTypes]
    if len(buildTypes) == 1 and '*' in buildTypes[0]:
        where = 'glob ?'
        glob_match = True
    else:
        qmarks = ','.join(['?'] * len(buildTypes))
        where = 'in (%s)' % qmarks
        glob_match = False

    # Find the last n revisions
    db = _db_module.TestResultsDatabase('nightly', verbose=False)
    cur = db.conn.execute('select distinct revision from summary '
                          'where build_type %s order by revision desc limit ?' % where,
                          tuple(buildTypes + [n]))
    revisions = [row[0] for row in cur]
    revisions.reverse()

    # Extract the build and run time information.  We create the following maps:
    #     builds maps (revision, build_type) -> [machine]
    #     timestamps maps (revision, build_type, machine) -> [timestamp]
    #     run_times maps suite_name -> {(revision, build_type, machine) -> (run_time, status)}
    # TODO #2518 The details table needs timestamp info!
    builds = {}
    timestamps = {}
    run_times = {}
    inf_test_names = ['Copyrights', 'DuplicateFileNames', 'OrphanedTests', 'Schemas']
    cur = db.conn.execute('select revision, machine, build_type, suite_name, suite_status, run_time from details'
                          ' where build_type %s and revision between ? and ?' % where,
                          tuple(buildTypes + [revisions[0], revisions[-1]]))
    for row in cur:
        # The builds dictionary
        k = (row['revision'], row['build_type'])
        machine = _canonical_machine_name(row['machine'])
        builds.setdefault(k, set()).add(machine)
        # The timestamps dictionary
        k_ts = k + (machine,)
        cur_ts = db.conn.execute('select finished from summary where revision=? and build_type=? and machine=?', k_ts)
        timestamps[k_ts] = sorted(db.FormatTimestamp(row['finished']) for row in cur_ts)
        # The run_times dictionary
        if row['suite_name'] not in inf_test_names:
            if not row['suite_name'] in run_times:
                run_times[row['suite_name']] = {}
            k = (row['revision'], row['build_type'], machine)
            #2518 TODO: Cope with multiple runs for these key combos
            run_times[row['suite_name']][k] = (row['run_time'], row['suite_status'])
    buildTypes = set()
    for k in builds:
        builds[k] = list(builds[k])
        builds[k].sort()
        buildTypes.add(k[1])
    buildTypes = list(buildTypes)
    output.append('\n\n<!-- Raw data:\n\n%s\n\n%s\n\n-->\n\n'
                  % (str(builds).replace(', (', ',\n ('),
                     str(run_times).replace('}, u', '},\n u')))

    # Display table headings
    buildTypes.sort()
    output.append('<table border="1">\n  <tr><th>Revision</th>\n')
    revbts = []
    for revision in revisions:
        cols = sum(map(lambda bt: len(builds.get((revision, bt), [])), buildTypes))
        if cols > 0:
            output.append('    <th colspan="%d">%s</th>\n'
                          % (cols, _linkChangeset(revision)))
        for bt in buildTypes:
            revbts.append((revision, bt))
    output.append('  </tr>\n  <tr><th>Build</th>\n')
    for rev, bt in revbts:
        if builds.has_key((rev, bt)):
            output.append('    <th colspan="%d">%s</th>\n'
                          % (len(builds[(rev, bt)]), _linkBuildType(bt, rev, wrappableText=True)))
    output.append('  </tr>\n  <tr><th>Machine</th>\n')
    for rev, bt in revbts:
        for machine in builds.get((rev, bt), []):
            output.append('    <th>%s</th>\n' %
                          _linkSummary(machine, 'nightly', rev, machine, bt, timestamps[(rev, bt, machine)][-1]))
    output.append('  </tr>\n')
    # Display the run times
    _handle_renamed_test_suites(run_times)
    test_suites = run_times.keys()
    test_suites.sort()
    for test_suite in test_suites:
        output.append('  <tr><th>%s</th>\n' % test_suite)
        for rev, bt in revbts:
            for machine in builds.get((rev, bt), []):
                k = (rev, bt, machine)
                if k in run_times[test_suite]:
                    run_time, status = run_times[test_suite][k]
                    link_text = _formatRunTime(run_time)
                    timestamp = timestamps[k][-1]
                    if bt.startswith('GoogleProfile'):
                        # _linkGraph includes the <td> tag
                        entry = _linkGraph('nightly', rev, machine, bt, timestamp,
                                           test_suite + 'Runner.gif', linkText=link_text)
                    else:
                        entry = _linkTestSuite('nightly', rev, machine, bt, timestamp, test_suite,
                                               status, run_time, None, linkText=link_text)
                        entry = '<td>%s</td>' % entry
                    output.append('    %s\n' % entry)
                else:
                    output.append('    <td></td>\n')
        output.append('  </tr>\n')
    output.append('</table>\n')
    
    # Graphs
    machines = set()
    for ms in builds.itervalues():
        machines.update(ms)
    gurl = '%s/profileHistoryGraph' % _our_url
    for machine in machines:
        for build_type in buildTypes:
            for test_suite in test_suites:
                # Graph plots run time against revision, and needs a list of (rev,time) pairs
                graph_data = []
                for r in revisions:
                    try:
                        graph_data.append((r, run_times[test_suite][(r, build_type, machine)][0]))
                    except KeyError:
                        pass
                graph_data_str = '|'.join(map(lambda p: ','.join(map(str, p)), graph_data))
                img_url = '%s?machine=%s&amp;buildType=%s&amp;testSuite=%s&amp;data=%s' % (gurl, machine, build_type, test_suite, graph_data_str)
                if len(img_url) > 2048:
                    # Internet Explorer will complain
                    img_url = '%s?machine=%s&amp;buildType=%s&amp;testSuite=%s&amp;n=%d' % (gurl, machine, build_type, test_suite, n)
                output.append('<h4>%s.%s on %s</h4>\n' % (build_type, test_suite, machine))
                output.append('<img src="%s" />\n' % img_url)
    
    return ''.join(output)

def _canonical_machine_name(machine):
    """Handle the .comlab -> .cs DNS change, mapping old results to the new names."""
    return machine.replace('.comlab.ox.ac.uk', '.cs.ox.ac.uk')

def _machine_name_aliases(machine):
    """Handle the .comlab -> .cs DNS change, returning a tuple of possible names for the given machine."""
    names = [machine]
    if machine.endswith('.cs.ox.ac.uk'):
        names.append(machine.replace('.cs.ox.ac.uk', '.comlab.ox.ac.uk'))
    return tuple(names)

def _handle_renamed_test_suites(runTimes, graphs={}):
    """Cope with the test result naming convention change in #2195.
    
    Test suite results files used to be named after just the leaf name of the test file.  They
    now use the full path.  In an attempt to make the history display still useful across this
    change, we change test suite names that have no path info to match one with the same leaf
    name but full path info.
    """
    name_map = {}
    for suite_name in runTimes.iterkeys():
        if '-' in suite_name:
            leaf_name = suite_name[1+suite_name.rfind('-'):]
            if not leaf_name in name_map and (leaf_name in runTimes or leaf_name in graphs):
                name_map[leaf_name] = suite_name
    for leaf_name, full_name in name_map.iteritems():
        if leaf_name in runTimes:
            runTimes[full_name].update(runTimes[leaf_name])
            del runTimes[leaf_name]
        if leaf_name in graphs:
            graphs[full_name] = graphs[leaf_name]
            del graphs[leaf_name]

def _test_suite_name_aliases(suiteName):
    """Cope with the test result naming convention change in #2195.
    
    If the given suiteName has path info, return it and its leaf name.  Otherwise just return
    suiteName.
    """
    if '-' in suiteName:
        return (suiteName, suiteName[1+suiteName.rfind('-'):])
    return (suiteName,)

def profileHistoryGraph(req, buildType, machine, testSuite, data='', n=''):
    """Show runtime graph for a specific testSuite in a profile build."""
    # Extract run-time data
    run_times = []
    if not data:
        if not n or not _db_module:
            return _error('Not enough data to plot, or no database available.')
        db = _db_module.TestResultsDatabase('nightly', verbose=False)
        aliases = _test_suite_name_aliases(testSuite)
        suite_test = 'suite_name in (%s)' % ','.join(['?'] * len(aliases))
        machines = _machine_name_aliases(machine)
        machine_test = 'machine in (%s)' % ','.join(['?'] * len(machines))
        cur = db.conn.execute('select revision, run_time from details where build_type=? and ' + machine_test
                              + ' and ' + suite_test + ' order by revision desc limit ?',
                              (buildType,) + machines + aliases + (n,))
        run_times = [(row['revision'], row['run_time']) for row in cur]
        run_times.reverse()
    else:
        for item in data.split('|'):
            run_times.append(map(float, item.split(',')))
    
    # Draw graph and send to browser
    req.content_type = 'image/png'
    from pychart import theme, axis, area, line_style, line_plot, canvas, category_coord, color
    theme.use_color = True
    theme.scale_factor = 2
    theme.title = 'Run time graph for ' + testSuite
    theme.reinitialize()
    xaxis = axis.X(format="/a-90/vM%d", label='Revision')
    yaxis = axis.Y(label='Run time (s)')
    ar = area.T(x_axis=xaxis, y_axis=yaxis, legend=None, size=(len(run_times)*8,120),
                x_coord=category_coord.T(run_times, 0))
    style = line_style.T(color=color.red)
    plot = line_plot.T(data=run_times, line_style=style)
    ar.add_plot(plot)
    c = canvas.init(req, format='png')
    ar.draw(c)
    c.close()
 
#####################################################################
##                    Helper functions.                            ##
#####################################################################

def _importModuleFromSvn(module_name, module_filepath, revision=None):
    """
    Use svn and imp to import the requested revision of the given
        module from the repository.
    module_name is the name to give the module.
    module_filepath is the path to the module file within the trunk
        directory of the repository.
    By default import the latest version.
    Return the module object.
    """
    filepath = _svn_repos + module_filepath
    command = ["svn", "cat"]
    if revision is not None:
        command.extend(["-r", str(revision)])
    command.extend(["--config-dir", "/home/svn/.subversion", filepath])
    stdin, stdout, stderr = os.popen3(command)
    module_text = ''.join(stdout.readlines())
    stdin.close()
    stderr.close()
    stdout.close()
    return _importCode(module_text, module_name)

def _importBuildTypesModule(revision=None):
    """
    Use svn and imp to import the requested revision of the BuildTypes.py
    module from the repository.
    By default import the latest version.
    """
    return _importModuleFromSvn('BuildTypes', '/python/BuildTypes.py', revision)

def _importCode(code, name, add_to_sys_modules=0):
    """
    Import dynamically generated code as a module. code is the
    object containing the code (a string, a file handle or an
    actual compiled code object, same types as accepted by an
    exec statement). The name is the name to give to the module,
    and the final argument says wheter to add it to sys.modules
    or not. If it is added, a subsequent import statement using
    name will return this module. If it is not added to sys.modules
    import will try to load it in the normal fashion.
    Code from the Python Cookbook.
    
    import foo
    
    is equivalent to
    
    foofile = open("/path/to/foo.py")
    foo = importCode(foofile,"foo",1)
    
    Returns a newly generated module.
    """
    import sys, imp
    
    module = imp.new_module(name)
    
    exec code in module.__dict__
    if add_to_sys_modules:
        sys.modules[name] = module
        
    return module

def _getBuildObject(buildTypesModule, buildType):
    """Get the build object for the given build type string."""
    if buildType.startswith('acceptance'):
        build = buildTypesModule.GetBuildType('default' + buildType[10:])
    elif buildType.startswith('longacceptance'):
        build = buildTypesModule.GetBuildType('default' + buildType[14:])
    elif buildType == 'failing':
        build = buildTypesModule.GetBuildType('default')
    else:
        build = buildTypesModule.GetBuildType(buildType)
    return build

def _extractDotSeparatedPair(string):
    """
    Extract both parts from a string of the form part1.part2.
    The '.' used is the last in the string.
    Returns a pair (part1, part2).
    Useful for parsing machine.buildType filenames.
    """
    i = string.rfind('.')
    return string[:i], string[i+1:]

def _extractTildeSeparatedPair(string):
    """Extract both parts from a string of the form revision~timestamp, where the ~timestamp component is optional.
    
    If no timestamp is present, returns (revision, None).
    """
    i = string.find('~')
    if i == -1:
        return string, None
    else:
        return string[:i], string[i+1:]

def _getResultsParentDirs(type, revision, timestamp=None):
    """Get the folder(s) holding test results for this test type, revision, and optionally timestamp.
    
    This deals with legacy results where the timestamp information is not included in the path, and
    the new style where we have revision~timestamp folders.
    If timestamp is None, a list of all folders for the given revision will be returned, newest first.
    If timestamp is specified then we check first for a folder with that timestamp, and secondly for a
    folder with no timestamp, returning a singleton list in either case.
    """
    base_path = os.path.join(_tests_dir, type, str(revision))
    if timestamp is None:
        results = glob.glob(base_path + '*')
        results.sort(key=lambda n: -float(n[n.find('~')+1:]) if '~' in n else 0)
    else:
        if os.path.isdir(base_path+'~'+str(timestamp)):
            results = [base_path+'~'+str(timestamp)]
        else:
            results = [base_path]
    return results

def _testResultsDir(type, revision, machine, buildType, timestamp=None):
    """Get the full path to a test results folder.
    This can cope with build machine renamings, using _machine_name_aliases.
    """
    for machine in reversed(_machine_name_aliases(machine)):
        for base_path in _getResultsParentDirs(type, revision, timestamp):
            results_path = os.path.join(base_path, machine+'.'+buildType)
            if os.path.isdir(results_path):
                break
    return results_path

_testSummaryRegexp = re.compile(r' *Overall status: <span style="color: (\w+);">(.*)</span>')
def _getTestSummary(test_set_dir, build):
    """
    Return a summary of the status of tests in the given directory,
    as a tuple of strings (targets, overall_status, colour).

    Does this by parsing the index.html page in the directory, looking
    for the overall status line.
    If this file doesn't exist or parsing fails, will fall back to
    using _getTestStatus.
    Target information is retrieved from the info.log file, if present.
    """
    index_path = os.path.join(test_set_dir, 'index.html')
    parsed_ok = False
    if os.path.isfile(index_path):
        # Load & parse file
        index_file = file(index_path, 'r')
        if index_file:
            for line in index_file:
                m = _testSummaryRegexp.match(line)
                if m:
                    overall_status = m.group(2)
                    colour = m.group(1)
                    parsed_ok = True
                    break
            index_file.close()
        if parsed_ok:
            # Check for build failure
            overall_status, colour = _checkBuildFailure(test_set_dir, overall_status, colour)
    if not parsed_ok:
        overall_status, colour = _getTestStatus(test_set_dir, build, True)
    targets = _getBuildTargets(test_set_dir)
    return targets, overall_status, colour

def _getBuildTargets(resultsDir):
    """Get the targets that were requested on the build command line, if known.
    
    It parses the info.log file, if present, to obtain them.
    """
    targets = ''
    info_path = os.path.join(resultsDir, 'info.log')
    try:
        info_file = open(info_path)
        for line in info_file:
            if line.startswith('Targets:'):
                targets = line[8:].strip()
                break
        info_file.close()
    except:
        pass
    return targets


_sconstruct_traceback_re = re.compile(r'  File ".*SConstruct", line ')
def _checkBuildFailure(test_set_dir, overall_status, colour):
    """Check whether the build failed, and return a new status if it did."""
    found_semget = False
    found_done_building = False
    try:
        log = file(os.path.join(test_set_dir, 'build.log'), 'r')
        if 'Windows' in test_set_dir:
            return _checkWinBuildFailure(log, overall_status, colour)
        for line in log:
            if (line.startswith('scons: building terminated because of errors.')
                or line.strip().endswith('(errors occurred during build).')
                or _sconstruct_traceback_re.match(line)):
                overall_status = 'Build failed (check build log for ": ***").  ' + overall_status
                colour = 'red'
                break
            if line.startswith('Error validating server certificate') and colour == 'green':
                overall_status = 'Probably failed svn update.  ' + overall_status
                colour = 'orange'
            if not found_semget and 'semget failed for setnum' in line:
                overall_status = 'Semaphore error.  ' + overall_status
                if colour == 'green':
                    colour = 'orange'
                found_semget = True
            if not found_done_building and line.startswith('scons: done building targets.'):
                found_done_building = True
        else:
            if not found_done_building:
                # Something went wrong that we didn't spot
                overall_status = "Build failed (didn't complete).  " + overall_status
                colour = 'red'
        log.close()
    except:
        # Build log may not exists for old builds
        pass
    return overall_status, colour

def _checkWinBuildFailure(log, overallStatus, colour):
    """Variant of _checkBuildFailure for Windows builds."""
    build_summary_re = re.compile(r'=+ Build: (\d+) succeeded, (\d+) failed, (\d+) up-to-date, (\d+) skipped =+')
    test_summary_re = re.compile(r'\d+% tests passed, (\d+) tests failed out of (\d+)')
    found_summaries = 0
    for line in log:
        m = build_summary_re.match(line)
        if m:
            found_summaries += 1
            if int(m.group(2)) > 0:
                overallStatus = 'Build failed.  ' + overallStatus
                colour = 'red'
            continue
        m = test_summary_re.match(line)
        if m:
            found_summaries += 1
            break
    if found_summaries < 2:
        overallStatus = "Build didn't complete.  " + overallStatus
        colour = 'red'
    return overallStatus, colour

def _isWindows(build):
    """Return build.is_windows in a backwards-compatible way."""
    try:
        return build.is_windows
    except AttributeError:
        return False

def _getTestStatus(test_set_dir, build, summary=False):
    """
    Return the status for all tests in the given directory, and compute
    a summary status given the build type.
    Return a tuple (dict, string, string, dict, dict) where the first entry maps
    test suite names to a string describing their status, the second is
    the overall status, the third is the colour in which to display
    the overall status, the fourth is a dictionary of run times, and the fifth
    is a dictionary of graphical files (used for graphical profiling).
    
    If summary is given as True, don't generate or return the dictionaries.
    """
    if _isWindows(build):
        return _getWinTestStatus(test_set_dir, build, summary)
    ignores = ['.html', '.log']
    result_files = os.listdir(test_set_dir)
    testsuite_status, runtime, graphs = {}, {}, {}
    for filename in result_files:
        ext = os.path.splitext(filename)[1]
        if ext == '.gif':
            if not summary:
                testsuite = build.ParseGraphFilename(filename)
                graphs[testsuite] = filename
        elif not (ext in ignores or filename[0] == '.'):
            d = build.GetInfoFromResultsFileName(filename)
            testsuite = d['testsuite']
            testsuite_status[testsuite] = d['status']
            if not summary:
                runtime[testsuite] = d['runtime']
    if graphs:
        _handle_renamed_test_suites(runtime, graphs)
    overall_status, colour = _overallStatus(testsuite_status, build)

    # Check for build failure
    overall_status, colour = _checkBuildFailure(test_set_dir, overall_status, colour)
    
    if summary:
        return overall_status, colour
    else:
        return testsuite_status, overall_status, colour, runtime, graphs

def _getWinTestStatus(testSetDir, build, summary=False):
    """Equivalent of _getTestStatus for Windows CMake-based builds.
    
    In CMake setups, we expect only three files in the results folder: info.log, build.log, and a ctest output text file.
    The latter contains all the test results, including the test output for those tests which didn't pass.
    It will be named like: ${TestPack}TestOutputs_${DateTime}.txt{,.tmp}
    Lines look like:
        Start   1: TestArchivingRunner
  1/248 Test   #1: TestArchivingRunner ..................................................................   Passed    0.55 sec
        Start   9: TestExceptionRunner
  9/248 Test   #9: TestExceptionRunner ..................................................................***Failed    0.42 sec
[test output]
        Start  10: TestExecutableSupportRunner
 10/248 Test  #10: TestExecutableSupportRunner ..........................................................   Passed    1.19 sec
        Start  59: TestMixedDimensionMeshRunner
 59/248 Test  #59: TestMixedDimensionMeshRunner .........................................................***Exception: Other  1.48 sec

69% tests passed, 77 tests failed out of 248

Total Test time (real) = 7724.14 sec

The following tests FAILED:
          9 - TestExceptionRunner (Failed)
    """
    ctest_results = glob.glob(os.path.join(testSetDir, '*TestOutputs_*.txt*'))
    # Regular expressions matching the key lines in the output 
    start_re = re.compile(r'\s+Start\s+\d+: ')
    test_re = re.compile(r'\s*(\d+)/(\d+) Test\s+#\d+: (\S+) \.*\s*(Passed|\*\*\*\D+)\s*([0-9.]+) sec')
    summary_re = re.compile(r'\d+% tests passed, (\d+) tests failed out of (\d+)')
    # Parse the output(s)
    testsuite_status, runtime, graphs = {}, {}, {}
    tests_completed = True
    expected_num_tests = -1 # In case the summary line isn't present in the output
    for results in ctest_results:
        if results.endswith('.tmp'):
            tests_completed = False
        for line in file(results):
            m = test_re.match(line)
            if m:
                testsuite = m.group(3)
                if m.group(4) == 'Passed':
                    testsuite_status[testsuite] = 'OK'
                elif m.group(4).startswith('Exception'):
                    testsuite_status[testsuite] = 'Unknown'
                else:
                    testsuite_status[testsuite] = '1_1'
                #TODO instead: testsuite_status = build.EncodeStatus(test_output)  (#2016)
                runtime[testsuite] = float(m.group(5))
            # Sanity check
            m = summary_re.match(line)
            if m:
                expected_num_tests = int(m.group(2))
    overall_status, colour = _overallStatus(testsuite_status, build)
    if not tests_completed:
        overall_status = "Tests didn't complete.  " + overall_status
        colour = 'red'
    if expected_num_tests != len(testsuite_status):
        overall_status = "Expected %d tests but found %d.  %s" % (expected_num_tests, len(testsuite_status), overall_status)
    # Check for build failure
    overall_status, colour = _checkBuildFailure(testSetDir, overall_status, colour)
    if summary:
        return overall_status, colour
    else:
        return testsuite_status, overall_status, colour, runtime, graphs


def _overallStatus(statuses, build):
    """
    Given a dict mapping test suite name to its status,
    and the type of build performed, return the overall status.
    Return value is a pair, the first item of which is a string giving
    the number of failing test suites, and the second a colour name.
    """
    total = len(statuses)
    failed, warnings = 0, 0
    components = set()
    for testsuite, status in statuses.iteritems():
        colour = _statusColour(status, build)
        if colour == 'red':
            failed += 1
            try:
                component = testsuite.split('-')[0]
                components.add(component)
            except IndexError:
                pass
        elif colour == 'orange':
            warnings += 1
    if failed > 0:
        if warnings:
            warnstr = " (with %d warnings)" % warnings
        else:
            warnstr = ""
        if components:
            warnstr += " (" + ', '.join(components) + ")"
        result = "Failed %d out of %d test suites%s" % (failed, total, warnstr)
        colour = "red"
    elif warnings > 0:
        result = "Warnings on %d out of %d test suites" % (warnings, total)
        colour = "orange"
    else:
        if total == 0:
            result = "No test results found"
            colour = "red"
        else:
            result = "All %d tests run passed" % total
            colour = "green"
    return result, colour


def _statusColour(status, build):
    """
    Return the name of the colour in which this status string should be
    displayed, given that the build type was build.
    """
    try:
        return build.StatusColour(status)
    except AttributeError:
        # Backwards compatibility
        if build.IsGoodStatus(status):
            return 'green'
        else:
            return 'red'

# Regexprs and state strings for _parseBuildTimings
#Command execution time: heart/build/debug/src/io/ChasteParameters_3_3.hpp: 0.737745 seconds
_time_re = re.compile(r"Command execution time:(?: ([^:]+):)? ([0-9.]+) seconds")
_states = ['Other', 'Compile', 'Object dependency analysis', 'CxxTest generation', 'PyCml execution', 'Test running']
# For older scons versions we need to parse the line before a timing output to determine what happened.
# Note that the order must match the _states list, except that we skip the 'Other' state.
_state_res = map(re.compile,
                 [r"[^ ]*mpicxx ",
                  r"BuildTest\(\[",
                  r"cxxtest/cxxtestgen.py",
                  r"RunPyCml\(\[",
                  r"(r|R)unning '(.*/build/.*/Test.*Runner|python/test/.*\.py)'"])
# For newer scons versions the timing line includes the target that was created, so we parse that instead
_target_state_map = [lambda t, ext: ext in ['.so', '.o', '.os'] or t.endswith('Runner'), # Compile
                     lambda t, ext: ext == '.dummy',                                     # Obj dep analysis
                     lambda t, ext: t.endswith('Runner.cpp'),                            # CxxTest
                     lambda t, ext: ext in ['.hpp', '.cpp'] and 'cellml' in t,           # PyCml
                     lambda t, ext: ext == '.log'                                        # Test running
                     ]

def _parseBuildTimings(logfilepath):
    """Parse a build log file to determine timings.
    
    Returns a dictionary mapping activity to time (in seconds).
    """
    times = [0] * len(_states)
    new_style = False
    try:
        logfile = open(logfilepath, 'r')
        if 'Windows' in logfilepath:
            return _parseWinBuildTimings(logfile)

        state = 0
        for line in logfile:
            m = _time_re.match(line)
            if m:
                if m.group(1):
                    # New-style with target info
                    new_style = True
                    ext = os.path.splitext(m.group(1))[1]
                    for i, func in enumerate(_target_state_map):
                        if func(m.group(1), ext):
                            state = i+1
                            break
                times[state] += float(m.group(2))
                state = 0
            elif not new_style and state == 0:
                for i, regexp in enumerate(_state_res):
                    m = regexp.match(line)
                    if m:
                        state = i+1
                        break
        logfile.close()

        result = dict(zip(_states, times))
    except IOError:
        result = dict.fromkeys(_states, -1.0)
    result['Total'] = sum(times)
    return result

def _parseWinBuildTimings(logfile):
    """Variant of _parseBuildTimings for Windows builds."""
    res = {'Compile': re.compile(r'\d+>Time Elapsed (\d+):(\d+):([0-9.]+)'),
           'Test running': re.compile(r'.*?\.+.*?([0-9.]+) sec')}
    times = dict([(k, 0.0) for k in res])
    for line in logfile:
        for key, regexp in res.iteritems():
            m = regexp.match(line)
            if m:
                multiplier = 1
                for time_part in reversed(m.groups()):
                    times[key] += float(time_part) * multiplier
                    multiplier *= 60
                break
    times['Total'] = sum(times.values())
    return times


def _getCcFlags(build):
    """A wrapper for BuildType.CcFlags that copes with ValueError being raised.
    
    This happens with the GccOptNative build type, since the compile flags are machine-specific.
    """
    try:
        return build.CcFlags()
    except ValueError:
        return "[machine specific optimisations are applied]"


#####################################################################
##                   HTML helper functions.                        ##
#####################################################################

def _buildLinkQuery(**parameters):
    """Build the query string for a hyperlink using the given key-value parameters."""
    return '&amp;'.join('%s=%s' % pair for pair in sorted(parameters.iteritems()))

def _linkRecent(text, type, start, **filters):
    """Return a link tag to the recent tests page, starting at the given position."""
    query = _buildLinkQuery(type=type, start=start, **filters)
    return '<a href="%s/recent?%s">%s</a>' % (_our_url, query, text)

def _linkRevision(revision, changes=False):
    """Return a link tag to the source browser for this revision."""
    try:
        revision = int(revision)
    except:
        return revision
    link = '<a href="%s?rev=%d">%d</a>' % (_source_browser_url, revision, revision)
    if changes:
        link += ' [%s]' % _linkChangeset(revision, 'changes')
    return link

def _linkChangeset(revision, text=None):
    """Return a link tag to the changes in this revision."""
    try:
        revision = int(revision)
    except:
        return revision
    if text is None:
        text = str(revision)
    return '<a href="%schangeset/%d">%s</a>' % (_trac_url, revision, text)

def _linkBuildType(buildType, revision, wrappableText=False):
    """Return a link tag to the detailed info page for this build type."""
    try:
        revision = int(revision)
    except:
        return buildType
    link_text = buildType
    if wrappableText:
        link_text = link_text.replace(',', ', ')
    query = _buildLinkQuery(buildType=buildType, revision=revision)
    return '<a href="%s/buildType?%s">%s</a>' % (_our_url, query, link_text)

def _linkSummary(text, type, revision, machine, buildType, timestamp):
    """Return a link tag to the summary page for this set of tests. text is the text of the link."""
    revision = int(revision)
    query = _buildLinkQuery(type=type, revision=revision, machine=machine, buildType=buildType, timestamp=timestamp)
    return '<a href="%s/summary?%s">%s</a>' % (_our_url, query, text)

def _linkTestSuite(type, revision, machine, buildType, timestamp, testsuite,
                   status, runtime, build, linkText=None):
    """Return a link tag to a page displaying the output from a single test suite."""
    if type == 'standalone':
        filename = build.ResultsFileName(os.curdir, testsuite, status, runtime)
        link = '<a href="%s">%s</a>' % (filename, build.DisplayStatus(status))
    else:
        query = _buildLinkQuery(type=type, revision=revision, machine=machine, buildType=buildType, timestamp=timestamp,
                                testsuite=testsuite, status=status, runtime=runtime)
        if linkText is None:
            linkText = build.DisplayStatus(status)
        link = '<a href="%s/testsuite?%s">%s</a>' % (_our_url, query, linkText)
    return link

def _linkGraph(type, revision, machine, buildType, timestamp, graphFilename, linkText='Profile callgraph'):
    """Return a link tag in a <td> to the graphics file which contains the graph."""
    if graphFilename == '':
        link = ''
    elif type == 'standalone':
        link = '<td><a href="%s">%s</a></td>' % (graphFilename, linkText)
    else:
        revision = int(revision)
        query = _buildLinkQuery(type=type, revision=revision, machine=machine, buildType=buildType, timestamp=timestamp, graphName=graphFilename)
        link = '<td><a href="%s/graph?%s">%s</a></td>' % (_our_url, query, linkText)
    return link

def _formatRunTime(runtime):
    "Return a human-readable version of the given runtime (which is in s)."
    if runtime < 0:
        s = 'Unknown'
    elif runtime < 60:
        s = str(runtime) + 's'
    else:
        # Minutes & seconds
        s = "%d:%02d" % (runtime // 60, runtime % 60)
    return s

def _colourText(text, colour):
    "Return text in the given colour."
    return '<span style="color: %s;">%s</span>' % (colour, text)

def _error(msg):
    "Encapsulate an error message."
    return '<p class="error">%s</p>' % msg

def _header(title=""):
    """HTML page header."""
    if title:
        title = " - " + title
    header = """\
<html>
    <head>
        <title>Chaste Tests%s</title>
        <link rel="stylesheet" href="/style.css" type="text/css">
    </head>
    <body>
    <script type="text/javascript">
<!--
        function toggle_visibility(id) {
           var e = document.getElementById(id);
           if(e.style.display == 'block')
              e.style.display = 'none';
           else
              e.style.display = 'block';
        }
//-->
</script>
""" % title
    return header

def _footer():
    """HTML page footer."""
    footer = """\
    <hr />
    <a href="%s">Tests index page</a><br />
    <a href="%s">Chaste project website</a>
    </body>
</html>""" % (_our_url, _trac_url)
    return footer


#####################################################################
##                     Standalone version                          ##
#####################################################################

if __name__ == '__main__':
    # We're being run from the command line rather than from within mod_python
    # Arguments should be:
    #  1: the directory containing test output files
    #  2: the type of build (string, defaults to 'default')
    # We write an index.html file for the local run of tests in the given directory.
    _standalone = True
    import sys, BuildTypes, socket
    if len(sys.argv) < 2:
        print "Syntax error."
        print "Usage:",sys.argv[0],"<test output dir> [<build type>]"
        sys.exit(1)
    _dir = sys.argv[1]
    if len(sys.argv) > 2:
        _build_type = sys.argv[2]
    else:
        _build_type = 'default'
    _build = BuildTypes.GetBuildType(_build_type)
    _machine = socket.getfqdn()

    # Alter the configuration slightly
    _tests_dir = '.'
    _source_browser_url = 'https://chaste.cs.ox.ac.uk/trac/browser/'
    _trac_url = 'https://chaste.cs.ox.ac.uk/trac/'
    _our_url = 'https://chaste.cs.ox.ac.uk/tests.py'

    _fp = file(os.path.join(_dir, 'index.html'), 'w')

    print >>_fp,_header('Test Summary For Local ' + _build_type + ' Build')
    print >>_fp,'<h1>Test Summary For Local Build</h1>'
    print >>_fp,'<p>Displaying info for tests with output stored in',_dir,'</p>'
    print >>_fp,_summary(None, 'standalone', 'working copy', _machine, _build_type)
    print >>_fp,_footer()
    _fp.close()

    print "Test summary generated in", 'file://' + os.path.abspath(os.path.join(_dir, 'index.html'))
    print "Overall test status:", _overall_status

