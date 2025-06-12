#!/bin/bash

set -e

# bash /scripts/adjust/autoclass.sh
# bash /scripts/adjust/tampor.sh
# bash /scripts/adjust/normalize.sh
# bash /scripts/adjust/limma.sh
# bash /scripts/adjust/combat.sh
bash /scripts/adjust/icvae.sh
bash /scripts/adjust/vfae.sh

