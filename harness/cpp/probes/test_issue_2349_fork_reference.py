# SPDX-License-Identifier: Apache-2.0
"""Wire-format and diagnostic regressions for the independent fork probe."""
import unittest

from issue_2349_fork_reference import decode_record, differences


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


if __name__ == '__main__':
    unittest.main()
