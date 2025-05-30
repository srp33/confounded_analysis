#!/bin/bash

set -e

bash /scripts/adjust/limma.sh
#bash /scripts/adjust/tampor.sh
#bash /scripts/adjust/scale.sh
bash /scripts/adjust/combat.sh
#bash /scripts/adjust/confounded.sh
