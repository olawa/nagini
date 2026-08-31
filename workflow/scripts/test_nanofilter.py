import gzip
import tempfile
import unittest
from pathlib import Path

from nanofilter_fallback import filter_fastq, parse_len_range


class NanofilterFallbackTest(unittest.TestCase):
    def test_parse_len_range(self):
        self.assertEqual(parse_len_range("500-3000"), (500, 3000))
        self.assertEqual(parse_len_range("500..3000"), (500, 3000))
        self.assertEqual(parse_len_range("500:3000"), (500, 3000))
        self.assertEqual(parse_len_range("500-"), (500, None))
        self.assertEqual(parse_len_range("-3000"), (None, 3000))
        self.assertEqual(parse_len_range(""), (None, None))

    def test_filter_fastq_by_len_and_qv(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            in_fq = tmp / "input.fastq.gz"
            out_fq = tmp / "output.fastq.gz"

            # Write 3 reads:
            # Read 1: len 100 bp (too short for 500-3000), mean Q ~30
            # Read 2: len 1000 bp, mean Q ~10 (too low for Q20)
            # Read 3: len 1000 bp, mean Q ~30 (PASS)
            with gzip.open(in_fq, "wt") as handle:
                # Read 1: len 100, Q ~30 (ASCII '?')
                handle.write("@read1\n" + "A" * 100 + "\n+\n" + "?" * 100 + "\n")
                # Read 2: len 1000, Q ~10 (ASCII '+')
                handle.write("@read2\n" + "A" * 1000 + "\n+\n" + "+" * 1000 + "\n")
                # Read 3: len 1000, Q ~30 (ASCII '?')
                handle.write("@read3\n" + "A" * 1000 + "\n+\n" + "?" * 1000 + "\n")

            filter_fastq(in_fq, out_fq, min_qv=20.0, min_len=500, max_len=3000)

            with gzip.open(out_fq, "rt") as handle:
                lines = handle.readlines()

            # Only read3 should pass
            self.assertEqual(len(lines), 4)
            self.assertEqual(lines[0].strip(), "@read3")


if __name__ == "__main__":
    unittest.main()
