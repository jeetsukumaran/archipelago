#! /bin/sh

rm -rfv $(find . -name "*.pyc")
rm -rfv $(find . -name "__pycache__")
rm -rfv build
rm -rfv dist
