#!/usr/bin/env python3
"""
Main entry point for IDEAL-GENOM-QC CLI.

This module serves as the entry point for the command-line interface,
importing and delegating to the actual CLI implementation.
"""

import sys
from pathlib import Path

# Add the project root to Python path to ensure imports work
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from .core.cli import main

if __name__ == '__main__':
    sys.exit(main())