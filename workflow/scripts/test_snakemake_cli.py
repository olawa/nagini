import os
import sys
import tempfile
import unittest

import snakemake_cli as cli


class RunTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.log = os.path.join(self.tmp.name, "job.log")

    def tearDown(self):
        self.tmp.cleanup()

    def read(self):
        with open(self.log) as handle:
            return handle.read()

    def test_arguments_are_passed_as_strings_and_output_goes_to_the_log(self):
        seen = []

        def main(argv):
            seen.extend(argv)
            print("to stdout")
            print("to stderr", file=sys.stderr)

        self.assertEqual(cli.run(main, ["--n", 4, "--x"], [self.log]), 0)
        self.assertEqual(seen, ["--n", "4", "--x"])
        self.assertEqual(self.read(), "to stdout\nto stderr\n")

    def test_exit_with_a_message_is_a_failure_that_is_logged(self):
        def main(argv):
            sys.exit("BAM file not found: x.bam")

        self.assertEqual(cli.run(main, [], [self.log]), 1)
        self.assertIn("BAM file not found: x.bam", self.read())

    def test_exit_codes_pass_through(self):
        self.assertEqual(cli.run(lambda argv: sys.exit(0), [], [self.log]), 0)
        self.assertEqual(cli.run(lambda argv: sys.exit(None), [], [self.log]), 0)
        self.assertEqual(cli.run(lambda argv: sys.exit(3), [], [self.log]), 3)

    def test_an_exception_is_a_failure_with_its_traceback_in_the_log(self):
        def main(argv):
            raise ValueError("bad table")

        self.assertEqual(cli.run(main, [], [self.log]), 1)
        self.assertIn("ValueError: bad table", self.read())
        self.assertIn("Traceback", self.read())

    def test_without_a_log_output_stays_on_the_console(self):
        self.assertEqual(cli.run(lambda argv: None, [], None), 0)


if __name__ == "__main__":
    unittest.main()
