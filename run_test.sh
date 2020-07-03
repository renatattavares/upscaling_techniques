#!/bin/bash

cd ../home && ls
cd ../elliptic_case && python3 -m unittest code_tests/test_20_case.py
