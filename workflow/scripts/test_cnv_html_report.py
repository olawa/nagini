import json
import os
import tempfile
import unittest
from pathlib import Path

from cnv_html_report import create_report


class TestCnvHtmlReport(unittest.TestCase):

    def setUp(self):
        self.template_content = (
            "<!DOCTYPE html><html><head><style>{{ css }}</style></head>"
            "<body>"
            "<p>Sample: {{ metadata.sample }} | TC: {{ metadata.tc }} ({{ metadata.tc_method }})</p>"
            "<script>"
            "const originalTc = {{ metadata.tc if metadata.tc is not none else 1 }};"
            "{{ js }}"
            "</script></body></html>"
        )
        self.template_file = tempfile.NamedTemporaryFile("w", delete=False)
        self.template_file.write(self.template_content)
        self.template_file.close()

        self.json_file = tempfile.NamedTemporaryFile("w", delete=False)
        self.json_file.write(json.dumps([{"mock": "data"}]))
        self.json_file.close()

    def tearDown(self):
        os.unlink(self.template_file.name)
        os.unlink(self.json_file.name)

    def test_create_report_numeric_tc(self):
        wildcards = type("Wildcards", (), {"sample": "test_sample"})()
        html = create_report(
            template_filename=self.template_file.name,
            json_filename=self.json_file.name,
            css_files=[],
            js_files=[],
            show_table=False,
            extra_tables=[],
            tc=1.0,
            tc_method="tumor-only",
            wildcards=wildcards,
            wide_plot_width=2000,
        )
        self.assertIn("const originalTc = 1.0;", html)
        self.assertIn("Sample: test_sample | TC: 1.0 (tumor-only)", html)
        self.assertNotIn("const originalTc = ;", html)

    def test_create_report_none_tc(self):
        wildcards = type("Wildcards", (), {"sample": "test_sample"})()
        html = create_report(
            template_filename=self.template_file.name,
            json_filename=self.json_file.name,
            css_files=[],
            js_files=[],
            show_table=False,
            extra_tables=[],
            tc=None,
            tc_method="tumor-only",
            wildcards=wildcards,
            wide_plot_width=2000,
        )
        self.assertIn("const originalTc = 1;", html)
        self.assertNotIn("const originalTc = ;", html)


if __name__ == "__main__":
    unittest.main()
