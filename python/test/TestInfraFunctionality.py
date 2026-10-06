
"""Copyright (c) 2005-2026, University of Oxford.
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

import os
import unittest
import sys
import filecmp
import importlib.util


def _load_infra_module(name):
    """Import a python/infra script by path, regardless of cwd."""
    path = os.path.join(os.path.dirname(__file__), '..', 'infra', name + '.py')
    spec = importlib.util.spec_from_file_location(name, os.path.abspath(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestInfraFunctionality(unittest.TestCase):

    def TestRoundResultsFiles(self):
        original = 'python/test/data/input/rounding_input.txt'
        output = os.path.join(CHASTE_TEST_OUTPUT, 'rounding_output.txt')
        # 8 decimal places
        rc = os.system(sys.executable + ' python/infra/RoundResultsFiles.py 8 <' + original + ' > '+ output)
        self.assertEqual(rc, 0)
        self.assertTrue(filecmp.cmp(output,'python/test/data/output/rounding_output8.txt'))
        # 6 decimal places
        rc = os.system(sys.executable + ' python/infra/RoundResultsFiles.py 6 <' + original + ' > '+ output)
        self.assertEqual(rc, 0)
        self.assertTrue(filecmp.cmp(output,'python/test/data/output/rounding_output6.txt'))


class TestCheckReleaseNotes(unittest.TestCase):
    """Unit tests for the matching logic in CheckReleaseNotes.py."""

    @classmethod
    def setUpClass(cls):
        cls.m = _load_infra_module('CheckReleaseNotes')

    def TestTitleIssues(self):
        f = self.m.title_issues
        self.assertEqual(f('561 add gperftools profiling'), {561})
        self.assertEqual(f('#602 Suppress the PMIx leak'), {602})
        self.assertEqual(f('561: add a thing'), {561})
        self.assertEqual(f('fix #42 and #43 in solver'), {42, 43})
        self.assertEqual(f('Support PETSc 3.25'), set())
        self.assertEqual(f('Cleanup: remove unused #include statements'), set())
        self.assertEqual(f('Refactor cell-based force hierarchy'), set())

    def TestReferenced(self):
        page = (
            '- [#500](https://github.com/Chaste/Chaste/issues/500) headline\n'
            '- [#453](https://github.com/Chaste/Chaste/pull/453) boost\n'
            'Supports VTK 9.1 on Ubuntu 2026; see #607 too.\n'
            'Fragment https://github.com/Chaste/Chaste/pull/9#issuecomment-1 stays 9.\n'
        )
        self.assertEqual(self.m.referenced(page), {453, 500, 607, 9})

    def _pr(self, number, title='Some work', closes=frozenset(), labels=(),
            author='alice', merged_by='bob', participants=()):
        return {
            'number': number, 'title': title, 'closes': set(closes),
            'title_issues': self.m.title_issues(title),
            'labels': [{'name': name} for name in labels],
            'author': None if author is None else {'login': author},
            'mergedBy': None if merged_by is None else {'login': merged_by},
            'participants': set(participants),
        }

    def TestRecordedReason(self):
        refs = {500, 607, 624}
        self.assertIsNone(self.m.recorded_reason(self._pr(888, 'Unrelated'), refs))
        self.assertEqual(self.m.recorded_reason(self._pr(777, labels=[self.m.EXEMPT_LABEL]), refs),
                         'labelled')
        self.assertIn('PR #624', self.m.recorded_reason(self._pr(624), refs))
        self.assertIn('#500', self.m.recorded_reason(self._pr(999, closes={500}), refs))
        self.assertIn('#607', self.m.recorded_reason(self._pr(621, '#607 Join statements'), refs))

    def TestIsBot(self):
        self.assertTrue(self.m.is_bot(None))
        self.assertTrue(self.m.is_bot({'login': 'app/copilot-swe-agent', 'is_bot': True}))
        self.assertTrue(self.m.is_bot({'login': 'dependabot[bot]'}))
        self.assertTrue(self.m.is_bot({'login': 'Copilot'}))
        self.assertFalse(self.m.is_bot({'login': 'alice'}))

    def TestTaggees(self):
        org = {'alice', 'bob', 'carol'}
        # Author when a human org member; else the merger.
        self.assertEqual(self.m.taggees(self._pr(1, author='alice'), org), ['alice'])
        self.assertEqual(self.m.taggees(self._pr(2, author='outsider', merged_by='carol'), org),
                         ['carol'])
        bot = self._pr(3, author='app/copilot-swe-agent', merged_by='bob')
        bot['author']['is_bot'] = True
        self.assertEqual(self.m.taggees(bot, org), ['bob'])
        # Else org-member participants; else the backup list.
        self.assertEqual(self.m.taggees(
            self._pr(4, author='outsider', merged_by=None,
                     participants={'alice', 'outsider', 'Copilot'}), org), ['alice'])
        self.assertEqual(self.m.taggees(
            self._pr(5, author='outsider', merged_by=None, participants={'other'}), org),
            self.m.BACKUP_TAGGEES)
        # No member list: any human author is trusted.
        self.assertEqual(self.m.taggees(self._pr(6, author='outsider', merged_by='carol'), None),
                         ['outsider'])

    def TestPrsInBody(self):
        self.assertEqual(self.m.prs_in_body('| #12 | a | @x |\n| #34 | b | @y |\n'), {12, 34})
        self.assertEqual(self.m.prs_in_body(''), set())

    def TestSplice(self):
        s, e = self.m.TABLE_START, self.m.TABLE_END
        out = self.m.splice(f"Intro.\n\n{s}\nold\n{e}\n\nFooter.\n", "new \\1 table")
        self.assertIn("Intro.", out)
        self.assertIn("Footer.", out)
        self.assertIn(f"{s}\nnew \\1 table\n{e}", out)
        self.assertNotIn("old", out)
        appended = self.m.splice("Just intro.", "the table")
        self.assertTrue(appended.rstrip().endswith(e))
        self.assertIn("the table", appended)
