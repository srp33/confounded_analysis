#!/bin/bash

# Redirect to the BatchQC analysis in the batchqc subfolder
exec "$(dirname "$0")/batchqc/run_batchqc.sh" "$@"