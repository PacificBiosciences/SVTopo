#!/usr/bin/env python
"""
Global pytest configuration for the project.
This file sets up the Python path and provides shared fixtures.
"""

import sys
import os
from pathlib import Path

# Add the SVTopoVz package to the Python path
project_root = Path(__file__).parent
svtopovz_path = project_root / "SVTopoVz"
sys.path.insert(0, str(svtopovz_path))

# Also add the project root to the path for any other imports
sys.path.insert(0, str(project_root)) 