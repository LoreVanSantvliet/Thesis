#!/bin/bash
# Lore Van Santvliet 31/03/2022
# This script removes '\\n' characters in a file (e.g. Phaster output file) and replaces them by '\n', a newline character. It also replaces a line with '-' characters by spaces.

sed -i '' 's/\\n/\'$'\n/g' $1
sed -i '' 's/--/  /g' $1

