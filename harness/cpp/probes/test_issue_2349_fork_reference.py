# SPDX-License-Identifier: Apache-2.0
"""Wire-format and diagnostic regressions for the independent fork probe."""
import unittest
import tempfile
import json
import hashlib
from pathlib import Path

from issue_2349_fork_reference import decode_record, differences, validate_case_names, verify_archive_stamp, analyze_report, FORK, FIXTURE


class RecordTests(unittest.TestCase):
    def test_full_record_preserves_every_field(self):
        record = [1, 42, 0, 2, 3, 4, 2, 7, 8, 1, 0.25, 0, *range(12)]
        fields = decode_record(record, 'full')
        self.assertEqual(fields['result_topology'], [7, 8])
        self.assertEqual(fields['result_geometry'], [0.25])
        self.assertEqual(fields['mass_properties'], list(range(12)))
        self.assertEqual(sum(map(len, fields.values())), len(record) - 3)

    def test_rejects_trailing_and_short_records(self):
        record = [0, 1001, 0, 0, 0, 0, 0, 1001]
        self.assertEqual(decode_record(record, 'full')['mass_properties'], [])
        for invalid in [record + [1], record[:-1], [-1, *record], [1.5, *record]]:
            with self.subTest(record=invalid), self.assertRaises(ValueError):
                decode_record(invalid, 'full')

    def test_named_differences_keep_exact_integer_checks_and_strict_values(self):
        lhs = {'prefix_state': [0, 1], 'generated_geometry_samples': [1.0, 0.0]}
        rhs = {'prefix_state': [0, 2], 'generated_geometry_samples': [1.0 + 2e-9, 1e-30]}
        result = differences(lhs, rhs)
        self.assertEqual([(d['field'], d['index']) for d in result],
                         [('generated_geometry_samples', 0),
                          ('generated_geometry_samples', 1), ('prefix_state', 1)])
        self.assertEqual(differences({'sample': [1.0]}, {'sample': [1.0 + 5e-10]}), [])

    def test_missing_fields_and_lengths_cannot_pass(self):
        self.assertTrue(differences({'x': [1]}, {'x': [1, 2]}))
        self.assertTrue(differences({'x': []}, {}))


class EvidenceTests(unittest.TestCase):
    def test_report_requires_each_fixture_exactly_once(self):
        cases = {'a': {}, 'b': {}}
        validate_case_names([{'name': 'b'}, {'name': 'a'}], cases)
        for names in [[], ['a'], ['a', 'a', 'b'], ['a', 'b', 'unknown']]:
            with self.subTest(names=names), self.assertRaises(ValueError):
                validate_case_names([{'name': name} for name in names], cases)

    def test_analyzer_rejects_partial_reports_before_comparison(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            report = root / 'report.json'
            debug = root / 'rust.debug'
            debug.write_text('')
            cases = json.loads(FIXTURE.read_text())['cases']
            names = [case['name'] for case in cases]
            archives = {'libsmlib_api.a': 'a' * 64, 'libsmlib_core.a': 'b' * 64}
            for selected in [[], names[:-1], names + [names[0]], names + ['unknown']]:
                report.write_text(json.dumps({
                    'fork_commit': FORK, 'libraries': archives,
                    'build_stamp': {'fork_commit': FORK, 'libraries': archives},
                    'records': [{'name': name} for name in selected]}))
                with self.subTest(names=selected), self.assertRaisesRegex(ValueError, 'exactly once'):
                    analyze_report(report, debug, root / 'output.json')

    def test_archive_stamp_binds_every_archive_and_revision(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            libraries = [root / 'libsmlib_api.a', root / 'libsmlib_core.a']
            for library in libraries:
                library.write_bytes(library.name.encode())
            stamp = root / 'build.json'
            manifest = {'fork_commit': 'revision', 'libraries': {
                p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in libraries}}
            stamp.write_text(json.dumps(manifest))
            self.assertEqual(verify_archive_stamp(stamp, libraries, 'revision'), manifest)
            with self.assertRaises(ValueError):
                verify_archive_stamp(stamp, libraries, 'wrong-revision')
            del manifest['libraries'][libraries[0].name]
            stamp.write_text(json.dumps(manifest))
            with self.assertRaises(ValueError):
                verify_archive_stamp(stamp, libraries, 'revision')
            manifest['libraries'][libraries[0].name] = '0' * 64
            stamp.write_text(json.dumps(manifest))
            with self.assertRaises(ValueError):
                verify_archive_stamp(stamp, libraries, 'revision')


if __name__ == '__main__':
    unittest.main()
