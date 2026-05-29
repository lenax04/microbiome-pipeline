#!/usr/bin/env python3
# Integration tests for MicroSnake pipeline
# Author: Dawid X (dawidx1233)

import os
import subprocess
import unittest

class TestMicroSnake(unittest.TestCase):
    def setUp(self):
        self.project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
        
    def test_dry_run(self):
        """Test if the Snakemake pipeline is syntactically correct with a dry-run."""
        cmd = ["snakemake", "--dry-run", "--use-conda"]
        result = subprocess.run(cmd, cwd=self.project_dir, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, f"Dry-run failed: {result.stderr}")
        print("Dry-run check passed successfully!")

if __name__ == "__main__":
    unittest.main()
