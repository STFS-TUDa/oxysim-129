#!/bin/bash
# Small shell helpers used by the Alltest scripts.

printNiceHeader() {
    echo ""
    echo "======================================================================"
    echo " $*"
    echo "======================================================================"
}

printWarning() {
    echo "WARNING: $*" >&2
}
