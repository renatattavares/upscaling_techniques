#!/bin/bash

cd ../home/scientific && ./pymoab_update.sh
cd ../../test && python3 -m unittest code_tests/test_20_case.py
