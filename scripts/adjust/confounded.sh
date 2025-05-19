 #!/bin/bash

set -e

printf "\033[0;32mAdjusting the data with Confounded\033[0m\n"

# We also save a log for our loss chart
confounded /data/gse20194/unadjusted.csv -o /data/gse20194/confounded.csv -b batch -c 1000 -f /outputs/metrics/gse20194_confounded_log.csv -s "sigmoid" -w 1.0 -g 0.0001 -m 256 -i 10000 -a 2

confounded /data/gse24080/unadjusted.csv -o /data/gse24080/confounded.csv -b batch -c 1000
confounded /data/gse49711/unadjusted.csv -o /data/gse49711/confounded.csv -b  Class -c 1000
