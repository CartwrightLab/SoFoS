#!/bin/bash

echo "Running Integration Test: usage ..."

SOFOS="${1:-../sofos}"

"${SOFOS}" -h | grep -q '^Usage:'

if [[ $? -ne 0 ]]; then
	echo "  ERROR: No usage information found on -h."
	exit 1
fi

"${SOFOS}" -? 2>&1 >/dev/null | grep -q '^Usage:'

if [[ $? -ne 0 ]]; then
	echo "  ERROR: No usage information found on flag error."
	exit 1
fi
