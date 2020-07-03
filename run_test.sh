#!/bin/bash

cd ../home && ./pymoab_update.sh
cd ../elliptic_case && python3 -m unittest code_tests/test_20_case.py
