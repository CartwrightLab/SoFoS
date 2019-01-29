#!/bin/bash

SOFOS="${1:-../sofos}"

bash test-00.bash "${SOFOS}" || exit 1

bash test-01.bash "${SOFOS}" || exit 1

bash test-02.bash "${SOFOS}" || exit 1
