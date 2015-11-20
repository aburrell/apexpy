#!/bin/bash
set -ev

python setup.py check --strict --metadata --restructuredtext
check-manifest
flake8 src tests
