#!/bin/bash

cd ../home && ls -a
cd ../elliptic_case && python3 -m unittest code_tests/test_20_case.py
