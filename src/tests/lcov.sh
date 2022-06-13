#!/bin/bash
lcov --capture --directory . --output-file coverage.info && 
lcov --remove coverage.info '*/catch2/*' --output-file coverage.info && 
lcov --remove coverage.info '*/Modules/*' --output-file coverage.info && 
lcov --remove coverage.info '*tests*' --output-file coverage.info && 
lcov --remove coverage.info '/usr/*' --output-file coverage.info && 
lcov --list coverage.info
