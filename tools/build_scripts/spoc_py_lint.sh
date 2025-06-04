#!/bin/sh

set -eu  # strict mode: stop on error or undefined var

SRC_DIR=${1:-"./"}
PYCODESTYLE_CFG_PATH=${2:-"./"}

# 1. Validate source directory exists
if [ ! -d "$SRC_DIR" ]; then
  echo "ERROR: Source directory not found: $SRC_DIR"
  exit 1
fi

# 2. Validate config file exists
if [ ! -f "$PYCODESTYLE_CFG_PATH/.pycodestyle" ]; then
  echo "ERROR: .pycodestyle not found in: $PYCODESTYLE_CFG_PATH"
  exit 1
fi

echo "Running pycodestyle in: $SRC_DIR"
echo "Using config: $PYCODESTYLE_CFG_PATH/.pycodestyle"

# 3. Find Python files recursively under SRC_DIR
PY_FILES=$(find "$SRC_DIR" -name '*.py' -o -name '*.py.in')

if [ -z "$PY_FILES" ]; then
  echo "No Python files found in: $SRC_DIR"
  exit 0
fi

# 4. Run pycodestyle
pycodestyle -v --config="$PYCODESTYLE_CFG_PATH/.pycodestyle" $PY_FILES

