#!/usr/bin/env python3
"""
Verification script: Checks all repository files for forbidden dash characters.
Rules strictly forbid the use of em dashes (U+2014) and en dashes (U+2013)
in all code, comments, documentation, and metadata files.
"""

import sys
import os

FORBIDDEN_CHARS = {
    chr(0x2014): 'Em dash (U+2014)',
    chr(0x2013): 'En dash (U+2013)',
}

EXCLUDE_DIRS = {'.git', '_build', '__pycache__', '.pytest_cache', 'bin', 'obj'}

def scan_file(filepath):
    violations = []
    try:
        with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
            for line_no, line in enumerate(f, start=1):
                for char, name in FORBIDDEN_CHARS.items():
                    if char in line:
                        violations.append((line_no, name, line.strip()))
    except Exception as e:
        print(f"Warning: could not read {filepath}: {e}", file=sys.stderr)
    return violations

def main():
    root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    has_violations = False
    checked_count = 0

    for dirpath, dirnames, filenames in os.walk(root_dir):
        dirnames[:] = [d for d in dirnames if d not in EXCLUDE_DIRS]
        for fname in filenames:
            ext = os.path.splitext(fname)[1].lower()
            if ext in {'.fut', '.c', '.h', '.py', '.md', '.txt', '.json', '.yaml', '.yml', '.sh', '.html', '.svg', '.csv'}:
                filepath = os.path.join(dirpath, fname)
                violations = scan_file(filepath)
                checked_count += 1
                if violations:
                    has_violations = True
                    rel_path = os.path.relpath(filepath, root_dir)
                    print(f"Violation in {rel_path}:")
                    for line_no, name, line in violations:
                        print(f"  Line {line_no} [{name}]: {line[:80]}")

    if has_violations:
        print("\nVerification FAILED: Forbidden dash characters detected.")
        sys.exit(1)
    else:
        print(f"Verification PASSED: {checked_count} files scanned, 0 forbidden dash characters found.")
        sys.exit(0)

if __name__ == '__main__':
    main()
