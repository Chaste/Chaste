
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


# Chaste tests display script.

# This module contains most of the functionality, and is loaded from the
# repository by a wrapper script.

import os
import operator
import time
import itertools
import re

# Compatability with Python 2.3
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
        output.append('\n    <li><a href="%s/recent?type=%s">%s</a></li>'
                        % (_our_url, tests_type, text))
    output.append("""</ul>
    <p>Branch/project builds: (<a style='text-decoration: underline; color: blue;' onclick="toggle_visibility('branch-list');">toggle visibility</a>)</p>
    <ul id='branch-list' style='display:none;'>""")
    for tests_type in branch_types:
        output.append('\n    <li><a href="%s/recent?type=%s">Recent %s builds.</a></li>' %
                      (_our_url, tests_type, tests_type))
    output.append("""</ul>
    
    <a href="%s/profileHistory">Run time variation of profiled tests.</a>
    
    <h2>Latest continuous build</h2>
""" % _our_url)

    # Look for the latest revision present.
    type = 'continuous'
    revisions = os.listdir(os.path.join(_tests_dir, type))
    revision = str(max(itertools.imap(int, revisions)))
    # Display summary of each machine & build type combination for this revision
    test_set_dir = os.path.join(_tests_dir, type, revision)
    builds = os.listdir(test_set_dir)
    if len(builds) < 1:
        output.append(_error('No test set found for revision '+revision+
                             '. Probably the build is still in progress.'))
        output.append('<p><a href="/out/latest">Latest build log.</a></p>')
    else:
        for build in builds:
            machine, buildType = _extractDotSeparatedPair(build)
            output.append(_summary(req, type, revision, machine, buildType))

    output.append(_footer())

    return ''.join(output)


def testsuite(req, type, revision, machine, buildType, testsuite, status, runtime):
    """
    Display the results for the given testsuite, by passing the file back
    to the user.
    """
    req.content_type = 'text/html'
    req.write(_header(), 0)
    test_set_dir = _testResultsDir(type, revision, machine, buildType)
    buildTypesModule = _importBuildTypesModule(revision)
    build = _getBuildObject(buildTypesModule, buildType)
    testsuite_file = build.ResultsFileName(test_set_dir, testsuite, status, runtime)
    if os.path.isfile(testsuite_file):
        req.write('\n<pre>\n', 0)
        fp = open(testsuite_file)
        for line in fp:
            req.write(line.replace('&', '&amp;').replace('<', '&lt;'))
        fp.close()
        req.write('\n</pre>\n', 0)
    else:
        req.write(_error('The requested test suite was not found.'))
    req.write(_footer())


def graph(req, type, revision, machine, buildType, graphName):
    """
    Send the given graph file (a .gif) back to the user.
    """
    test_set_dir = _testResultsDir(type, revision, machine, buildType)
    graphName = os.path.basename(graphName) # Foil hackers
    graph_path = os.path.join(test_set_dir, graphName)
    if os.path.isfile(graph_path):
        req.content_type = 'image/gif'
        req.sendfile(graph_path)
    else:
        req.content_type = 'text/html'
        req.write(_header(), 0)
        req.write(_error('The requested graph file was not found.'))
        req.write(_footer())


def recent(req, type='', start=0, buildType=''):
    "User-facing page. Content is generated by _recent."
    title = 'Recent '+type+' builds'
    page_body = """\
    <h1>%s</h1>
%s
""" % (title, _recent(req, type, int(start), buildType))
    return _header(title) + page_body + _footer()

def _recent(req, type='', start=0, buildType=''):
    """Display brief summaries of recent builds of the given type.
    
    Returns a string representing part of a webpage.
    """
    if not type:
        return _error('No type of test to summarise specified.')

    dir = os.path.join(_tests_dir, type)
    if not os.path.isdir(dir):
        return _error(type+' is not a valid type of test.')
    
    n_per_page = 20 # Number of builds to display per page
    if _db_module and not req.form.getfirst('nocache', False):
        # Get information from the database
        db = _db_module.TestResultsDatabase(type, verbose=False)
        db.FastUpdate()
        total_num_of_builds = db.CountResults()
        if not buildType:
            where = ''
            params = (n_per_page, start)
        else:
            where = ' where build_type=?'
            params = (buildType, n_per_page, start)
        cur = db.conn.execute('select * from summary%s'
                              ' order by finished desc, revision desc, build_type, machine'
                              ' limit ? offset ?' % where, params)
        def gen_row(cur=cur):
            for row in cur:
                yield (row['revision'], row['machine'], row['build_type'],
                       time.mktime(row['finished'].timetuple()), row['status'], row['colour'])
    else:
        # Parse the directory structure within dir into a list of builds
        builds = []
        for revision in os.listdir(dir):
            for machine_and_build_type in os.listdir(os.path.join(dir, revision)):
                st = os.stat(os.path.join(dir, revision, machine_and_build_type))
                mod_time = st.st_mtime
                machine, build_type = _extractDotSeparatedPair(machine_and_build_type)
                if not buildType or build_type == buildType:
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
                test_set_dir = _testResultsDir(type, revision, machine, build_type)
                overall_status, colour = _getTestSummary(test_set_dir, build)
                yield (revision, machine, build_type, finished, overall_status, colour)

    output = []
    if start > 0:
        output.append(_linkRecent('Previous page', type, start=start-n_per_page, buildType=buildType) + " ")
    if buildType:
        output.append(_linkRecent('All build types', type, start=start) + " ")
    if total_num_of_builds > start+n_per_page:
        output.append(_linkRecent('Next page', type, start=start+n_per_page, buildType=buildType))
    output.append("""\
    <table border="1">
        <tr>
          <th>Date</th>
          <th>Revision</th>
          <th>Build Type</th>
          <th>Machine</th>
          <th>Status</th>
        </tr>
""")

    bgcols = ["white", "#eedd82"]
    bgcol_index = 0
    old_revision = -1
    for revision, machine, build_type, finished, overall_status, colour in gen_row():
        if type == 'nightly':
            date = time.strftime('%d/%m/%Y', time.localtime(finished))
        else:
            date = time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(finished))
        if revision != old_revision:
            bgcol_index = 1 - bgcol_index
            old_revision = revision
        subs = {'bgcol': bgcols[bgcol_index], 'status_col': colour,
                'date': date, 'machine': machine,
                'rev': _linkRevision(revision),
                'build_type': _linkBuildType(build_type, revision),
                'status': _linkSummary(overall_status, type, revision, machine, build_type)}
        if not buildType:
            subs['build_type'] += ' [' + _linkRecent('filter', type, start, build_type) + ']'
        output.append("""\
        <tr>
          <td style='background-color: %(bgcol)s;'>%(date)s</td>
          <td style='background-color: %(bgcol)s;'>%(rev)s</td>
          <td style='background-color: %(bgcol)s;'>%(build_type)s</td>
          <td style='background-color: %(bgcol)s;'>%(machine)s</td>
          <td style='background-color: %(status_col)s;'>%(status)s</td>
        </tr>
""" % subs)
    output.append("  </table>\n")

    if start > 0:
        output.append(_linkRecent('Previous page', type, start=start-n_per_page, buildType=buildType) + " ")
    if buildType:
        output.append(_linkRecent('All build types', type, start=start) + " ")
    if total_num_of_builds > start+n_per_page:
        output.append(_linkRecent('Next page', type, start=start+n_per_page, buildType=buildType))

    if type == 'nightly':
        output.append('<p><a href="/out/latest-nightly">Latest nightly build log.</a></p>')
    elif type == 'continuous':
        output.append('<p><a href="/out/latest">Latest continuous build log.</a></p>')
    
    return ''.join(output)



def summary(req, type, revision, machine, buildType):
    "User-facing page. Content is generated by _summary."
    page_body = """\
    <h1>Build Summary</h1>
""" +  _summary(req, type, revision, machine, buildType)
    return _header('Test Summary') + page_body + _footer()
    
def _summary(req, type, revision, machine=None, buildType=None):
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
        test_set_dir = _testResultsDir(type, revision, machine, buildType)
    
    # Now test_set_dir should be the directory containing the test results
    # to summarise. Extract summary info from the filenames.
    if type == 'standalone':
        build = _build
    else:
        buildTypesModule = _importBuildTypesModule(revision)
        build = _getBuildObject(buildTypesModule, buildType)
    testsuite_status, overall_status, colour, runtime, graphs = _getTestStatus(test_set_dir, build)
    # Store overall status for the standalone script case
    if _standalone:
        global _overall_status
        _overall_status = overall_status

    # Get the timestamp on the directory
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
    %s
    <table border="1">
        <tr>
          <th>Test Suite</th>
          <th>Status</th>
          <th>Run Time</th>
          %s
        </tr>
""" % (_linkRevision(revision), date, _colourText(overall_status, colour),
           _linkBuildType(buildType, revision), machine, build_log, extra_cols))
    
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
                'status': _linkTestSuite(type, revision, machine, buildType,
                                         testsuite, testsuite_status[testsuite],
                                         runtime[testsuite], build),
                'runtime': _formatRunTime(runtime[testsuite]),
                'graph': _linkGraph(type, revision, machine, buildType,
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
    """
    Display information on the compiler settings, etc. used to build a set
    of tests.
    buildType is the user-friendly name describing these settings, such as
    can be passed to scons build=buildType.
    revision is the code revision of the set of tests, in case the 
    definition of buildType has changed since.
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
           build.CcFlags(), build.LinkFlags(),
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


def profileHistory(req, n=20):
    """Show runtimes for the last n profile builds."""
    page_body = """
    <h1>Profile History</h1>
""" +  _profileHistory(req, int(n))
    return _header('Profile History') + page_body + _footer()

def _profileHistory(req, n=20):
    """Show runtimes for the last n profile builds."""
    tests_dir = os.path.join(_tests_dir, 'nightly')
    if not os.path.isdir(tests_dir):
        return _error('No nightly tests found')

    # These are the build types representing profile builds
    build_types = ['Profile_ndebug', 'GoogleProfile_ndebug']

    # Find the last n revisions
    if _db_module and not req.form.getfirst('nocache', False):
        db = _db_module.TestResultsDatabase('nightly', verbose=False)
        db.FastUpdate()
        cur = db.conn.execute('select distinct revision from summary '
                              'where build_type in (?,?) order by revision desc limit ?',
                              tuple(build_types + [n]))
        revisions = [row[0] for row in cur]
        revisions.reverse()
    else:
        revisions = map(int, os.listdir(tests_dir))
        revisions.sort()
        revisions = revisions[-n:]
    
    # Extract the build and run information.  We create two maps: builds and run_times.
    # builds former maps (revision, build_type) -> [machine]
    # run_times maps suite_name -> {(revision, build_type, machine) -> (run_time, status)}
    # Find the appropriate builds for these revisions: map from (revision, build_type) -> [machine]
    builds = {}
    run_times = {}
    if _db_module:
        cur = db.conn.execute('select revision, machine, build_type, suite_name, suite_status, run_time from details'
                              ' where build_type in (?,?) and revision between ? and ?',
                              tuple(build_types + [revisions[0], revisions[-1]]))
        for row in cur:
            # The builds dictionary
            k = (row['revision'], row['build_type'])
            if not k in builds:
                builds[k] = set()
            builds[k].add(row['machine'])
            # The run_times dictionary
            if row['suite_name'][:4] == 'Test':
                if not row['suite_name'] in run_times:
                    run_times[row['suite_name']] = {}
                k = (row['revision'], row['build_type'], row['machine'])
                run_times[row['suite_name']][k] = (row['run_time'], row['suite_status'])
        for k in builds:
            builds[k] = list(builds[k])
            builds[k].sort()
    else:
        for revision in revisions:
            rev_dir = os.path.join(tests_dir, str(revision))
            for machine_and_build_type in os.listdir(rev_dir):
                machine, build_type = _extractDotSeparatedPair(machine_and_build_type)
                if build_type in build_types:
                    k = (revision, build_type)
                    if not builds.has_key(k):
                        builds[k] = []
                    builds[k].append(machine)
    
        build = FakeBuildType()
        for revision in revisions:
            for build_type in build_types:
                for machine in builds.get((revision, build_type), []):
                    k = (revision, build_type, machine)
                    d = _testResultsDir('nightly', revision, machine, build_type)
                    statuses, _, _, runtimes, _ = _getTestStatus(d, build)
                    for test_suite in statuses.keys():
                        if test_suite[:4] == 'Test':
                            if not run_times.has_key(test_suite):
                                run_times[test_suite] = {}
                            run_times[test_suite][k] = (runtimes[test_suite],
                                                        statuses[test_suite])

    # Display table headings
    output = ['<table border="1">\n  <tr><th>Revision</th>\n']
    revbts = []
    for revision in revisions:
        cols = sum(map(lambda bt: len(builds.get((revision, bt), [])), build_types))
        if cols > 0:
            output.append('    <th colspan="%d">%s</th>\n'
                          % (cols, _linkChangeset(revision)))
        for bt in build_types:
            revbts.append((revision, bt))
    output.append('  </tr>\n  <tr><th>Build</th>\n')
    for rev, bt in revbts:
        if builds.has_key((rev, bt)):
            output.append('    <th colspan="%d">%s</th>\n'
                          % (len(builds[(rev, bt)]), _linkBuildType(bt, rev)))
    output.append('  </tr>\n  <tr><th>Machine</th>\n')
    for rev, bt in revbts:
        for machine in builds.get((rev, bt), []):
            output.append('    <th>%s</th>\n' %
                          _linkSummary(machine, 'nightly', rev, machine, bt))
    output.append('  </tr>\n')
    # Display the run times
    test_suites = run_times.keys()
    test_suites.sort()
    for test_suite in test_suites:
        output.append('  <tr><th>%s</th>\n' % test_suite)
        for rev, bt in revbts:
            for machine in builds.get((rev, bt), []):
                k = (rev, bt, machine)
                if run_times[test_suite].has_key(k):
                    run_time, status = run_times[test_suite][k]
                    link_text = _formatRunTime(run_time)
                    if bt.startswith('GoogleProfile'):
                        # _linkGraph includes the <td> tag
                        entry = _linkGraph('nightly', rev, machine, bt,
                                           test_suite + 'Runner.gif', linkText=link_text)
                    else:
                        entry = _linkTestSuite('nightly', rev, machine, bt, test_suite,
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
        for build_type in build_types:
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

def profileHistoryGraph(req, buildType, machine, testSuite, data='', n=''):
    """Show runtime graph for a specific testSuite in a profile build."""
    # Extract run-time data
    run_times = []
    if not data:
        if not n or not _db_module:
            return _error('Not enough data to plot, or no database available.')
        db = _db_module.TestResultsDatabase('nightly', verbose=False)
        db.FastUpdate()
        cur = db.conn.execute('select revision, run_time from details where build_type=? and machine=? and suite_name=?'
                              ' order by revision desc limit ?',
                              (buildType, machine, testSuite, n))
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

def _importModuleFromSvn(module_name, module_filepath,
                             revision=None):
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

def _testResultsDir(type, revision, machine, buildType):
    """
    Return the directory in which test results are stored for this
    test type, code revision, build machine and build type.
    """
    return os.path.join(_tests_dir, type, str(revision), machine+'.'+buildType)

_testSummaryRegexp = re.compile(r' *Overall status: <span style="color: (\w+);">(.*)</span>')
def _getTestSummary(test_set_dir, build):
    """
    Return a summary of the status of tests in the given directory,
    as a tuple of strings (overall_status, colour).

    Does this by parsing the index.html page in the directory, looking
    for the overall status line.
    If this file doesn't exist or parsing fails, will fall back to
    using _getTestStatus.
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
    return overall_status, colour

def _checkBuildFailure(test_set_dir, overall_status, colour):
    """Check whether the build failed, and return a new status if it did."""
    found_semget = False
    try:
        log = file(os.path.join(test_set_dir, 'build.log'), 'r')
        for line in log:
            if (line.startswith('scons: building terminated because of errors.')
                or line.strip().endswith('(errors occurred during build).')
                or line.startswith('  File "SConstruct", line ')):
                overall_status = 'Build failed (check build log for ": ***").  ' + overall_status
                colour = 'red'
                break
            if not found_semget and 'semget failed for setnum' in line:
                overall_status = 'Semaphore error.  ' + overall_status
                if colour == 'green':
                    colour = 'orange'
                found_semget = True
        log.close()
    except:
        # Build log may not exists for old builds
        pass
    return overall_status, colour

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
    ignores = ['index.html', 'build.log']
    result_files = os.listdir(test_set_dir)
    testsuite_status, runtime, graphs = {}, {}, {}
    for filename in result_files:
        if filename[-3:]=='gif':
            testsuite = build.ParseGraphFilename(filename)
            graphs[testsuite] = filename
        elif not (filename in ignores or filename[0] == '.'):
            d = build.GetInfoFromResultsFileName(filename)
            testsuite = d['testsuite']
            testsuite_status[testsuite] = d['status']
            if not summary:
                runtime[testsuite] = d['runtime']
    overall_status, colour = _overallStatus(testsuite_status,
                                              build)

    # Check for build failure
    overall_status, colour = _checkBuildFailure(test_set_dir, overall_status, colour)
    
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
        try:
            colour = build.StatusColour(status)
        except AttributeError:
            # Backwards compatibility
            if build.IsGoodStatus(status):
                colour = 'green'
            else:
                colour = 'red'
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
            warnstr += " (" + ','.join(components) + ")"
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
_time_re = re.compile(r"Command execution time: ([0-9.]+) seconds")
_states = ['Other', 'Compile', 'Object dependency analysis', 'CxxTest generation', 'PyCml execution', 'Test running']
_state_res = map(re.compile,
                 [r"mpicxx ",
                  r"BuildTest\(\[",
                  r"cxxtest/cxxtestgen.py",
                  r"RunPyCml\(\[",
                  r"(r|R)unning '(.*/build/.*/Test.*Runner|python/test/.*\.py)'"])

def _parseBuildTimings(logfilename):
    """Parse a build log file to determine timings.
    
    Returns a dictionary mapping activity to time (in seconds).
    """
    times = [0] * len(_states)
    try:
        logfile = open(logfilename, 'r')
        state = 0
    
        for line in logfile:
            m = _time_re.match(line)
            if m:
                times[state] += float(m.group(1))
                state = 0
            elif state == 0:
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



#####################################################################
##                   HTML helper functions.                        ##
#####################################################################

def _linkRecent(text, type, start, buildType=''):
    "Return a link tag to the recent tests page, starting at the given position."
    query = 'recent?type=%s&amp;start=%d' % (type, start)
    if buildType:
        query += '&amp;buildType=' + buildType
    return '<a href="%s/%s">%s</a>' % (_our_url, query, text)

def _linkRevision(revision):
    "Return a link tag to the source browser for this revision."
    if revision == 'working copy':
        return revision
    return '<a href="%s?rev=%s">%s</a>' % (_source_browser_url,
                                             revision, revision)

def _linkChangeset(revision):
    """Return a link tag to the changes in this revision."""
    try:
        revision = int(revision)
    except:
        return revision
    return '<a href="%schangeset/%d">%d</a>' % (_trac_url, revision, revision)

def _linkBuildType(buildType, revision):
    "Return a link tag to the detailed info page for this build type."
    try:
        revision = int(revision)
    except:
        return buildType
    query = 'buildType?buildType=%s&revision=%d' % (buildType, revision)
    return '<a href="%s/%s">%s</a>' % (_our_url, query, buildType)

def _linkSummary(text, type, revision, machine, buildType):
    """
    Return a link tag to the summary page for this set of tests.
    text is the text of the link.
    """
    revision = int(revision)
    query = 'type=%s&revision=%d&machine=%s&buildType=%s' % (type, revision,
                                                               machine, buildType)
    return '<a href="%s/summary?%s">%s</a>' % (_our_url, query, text)

def _linkTestSuite(type, revision, machine, buildType, testsuite,
                       status, runtime, build, linkText=None):
    """
    Return a link tag to a page displaying the output from a single
    test suite.
    """
    if type == 'standalone':
        filename = build.ResultsFileName(os.curdir, testsuite, status, runtime)
        link = '<a href="%s">%s</a>' % (filename, build.DisplayStatus(status))
    else:
        query = 'type=%s&amp;revision=%s&amp;machine=%s&amp;buildType=%s'\
          % (type, revision, machine, buildType)
        query += '&amp;testsuite=%s&amp;status=%s&amp;runtime=%d'\
          % (testsuite, status, runtime)
        if linkText is None:
            linkText = build.DisplayStatus(status)
        link = '<a href="%s/testsuite?%s">%s</a>' % (_our_url, query, linkText)
    return link

def _linkGraph(type, revision, machine, buildType, graphFilename,
                   linkText='Profile callgraph'):
    """
    Return a link tag in a <td> to the graphics file which contains the graph.
    """
    if graphFilename == '':
        link = ''
    elif type == 'standalone':
        link = '<td><a href="%s">%s</a></td>' % (graphFilename, linkText)
    else:
        revision = int(revision)
        query = 'type=%s&amp;revision=%d&amp;machine=%s&amp;buildType=%s' % \
                (type, revision, machine, buildType)
        query = query + '&amp;graphName=%s' % (graphFilename)
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
    _source_browser_url = 'https://chaste.cs.ox.ac.uk/cgi-bin/trac.cgi/browser/'
    _trac_url = 'https://chaste.cs.ox.ac.uk/cgi-bin/trac.cgi/'
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

