# ==============================================================================
# SCRIPT 1: profile_adjust.sh
#
# This script automates the process of performance testing your adjust.R script.
# It creates subsets of your data, runs the adjustment, and logs the time.
# ==============================================================================
#!/bin/bash

# We remove 'set -e' to handle errors within the loop,
# which prevents the script from exiting prematurely if one R run fails.

printf "\033[0;32mStarting Performance Profiling...\033[0m\n"

# --- Configuration ---
ADJUST_SCRIPT="/scripts/adjust/adjust.R"
SUBSET_SCRIPT="scripts/adjust/subset_data.R" 
INPUT_FILE="/data/gold/gse49711/unadjusted.csv"
BATCH_COLUMN="meta_Sex" 
PRED_COLUMN="meta_er_status"
OUTPUT_DIR="/tmp"
PERFORMANCE_LOG="/outputs/performance_log.csv"
ROW_SIZES=(10 50 100 150 200 250 272)
COL_SIZES=(10 50 100 300 500 1000 5000 10000)

# --- Setup ---
mkdir -p $OUTPUT_DIR
echo "rows,cols,real_time,user_time,sys_time" > $PERFORMANCE_LOG
printf "Configuration complete. Output will be saved in '%s'\n" $OUTPUT_DIR

# --- Profiling by Number of Rows (Samples) ---
FIXED_COLS=50
printf "\n\033[0;34m--- Profiling by Number of Rows (Samples) ---\033[0m\n"
printf "Columns fixed at: %d\n" $FIXED_COLS

for N_ROWS in "${ROW_SIZES[@]}"; do
    TEMP_INPUT="$OUTPUT_DIR/subset_${N_ROWS}rows_${FIXED_COLS}cols.csv"
    TEMP_OUTPUT="$OUTPUT_DIR/adjusted_${N_ROWS}rows_${FIXED_COLS}cols.csv"
    
    printf "\nCreating subset: %d rows, %d columns...\n" $N_ROWS $FIXED_COLS
    Rscript $SUBSET_SCRIPT $INPUT_FILE $TEMP_INPUT -r $N_ROWS -c $FIXED_COLS

    printf "Running adjustment for subset with %d rows...\n" $N_ROWS
    
    ADJUST_CMD="Rscript $ADJUST_SCRIPT $TEMP_INPUT $TEMP_OUTPUT -b $BATCH_COLUMN -a fairadapt -c $PRED_COLUMN"
    
    # ========================= MODIFICATION START =========================
    # This portable method uses the shell's built-in 'time' and pipes its
    # output to 'tee'. 'tee' displays the output on the screen in real-time
    # AND saves it to a variable for parsing.

    TIME_OUTPUT=$((time -p $ADJUST_CMD) 2>&1 | tee /dev/tty)
    EXIT_CODE=${PIPESTATUS[0]} # Use PIPESTATUS to get the R script's exit code, not tee's

    if [ $EXIT_CODE -ne 0 ]; then
        # Handle the error. The full error output from the R script has
        # already been displayed on the screen by 'tee'.
        printf "\n\033[0;31mERROR: R script failed with exit code %d for %d rows, %d columns.\033[0m\n" $EXIT_CODE $N_ROWS $FIXED_COLS
        # Log the failure with a -1 for time
        echo "$N_ROWS,$FIXED_COLS,-1,-1,-1" >> $PERFORMANCE_LOG
    else
        # If successful, parse the time from the captured output and log it
        REAL_TIME=$(echo "$TIME_OUTPUT" | grep real | awk '{print $2}')
        USER_TIME=$(echo "$TIME_OUTPUT" | grep user | awk '{print $2}')
        SYS_TIME=$(echo "$TIME_OUTPUT" | grep sys | awk '{print $2}')

        printf "Time taken: %s seconds\n" $REAL_TIME
        echo "$N_ROWS,$FIXED_COLS,$REAL_TIME,$USER_TIME,$SYS_TIME" >> $PERFORMANCE_LOG
    fi
    # ========================== MODIFICATION END ==========================

    rm -f $TEMP_INPUT $TEMP_OUTPUT
done

# --- Profiling by Number of Columns (Features) ---
FIXED_ROWS=272
printf "\n\033[0;34m--- Profiling by Number of Columns (Features) ---\033[0m\n"
printf "Rows fixed at: %d\n" $FIXED_ROWS

for N_COLS in "${COL_SIZES[@]}"; do
    TEMP_INPUT="$OUTPUT_DIR/subset_${FIXED_ROWS}rows_${N_COLS}cols.csv"
    TEMP_OUTPUT="$OUTPUT_DIR/adjusted_${FIXED_ROWS}rows_${N_COLS}cols.csv"
    
    printf "\nCreating subset: %d rows, %d columns...\n" $FIXED_ROWS $N_COLS
    Rscript $SUBSET_SCRIPT $INPUT_FILE $TEMP_INPUT -r $FIXED_ROWS -c $N_COLS

    printf "Running adjustment for subset with %d columns...\n" $N_COLS
    
    ADJUST_CMD="Rscript $ADJUST_SCRIPT $TEMP_INPUT $TEMP_OUTPUT -b $BATCH_COLUMN -a fairadapt"
    
    # ========================= MODIFICATION START =========================
    TIME_OUTPUT=$((time -p $ADJUST_CMD) 2>&1 | tee /dev/tty)
    EXIT_CODE=${PIPESTATUS[0]}

    if [ $EXIT_CODE -ne 0 ]; then
        printf "\n\033[0;31mERROR: R script failed with exit code %d for %d rows, %d columns.\033[0m\n" $EXIT_CODE $FIXED_ROWS $N_COLS
        echo "$FIXED_ROWS,$N_COLS,-1,-1,-1" >> $PERFORMANCE_LOG
    else
        REAL_TIME=$(echo "$TIME_OUTPUT" | grep real | awk '{print $2}')
        USER_TIME=$(echo "$TIME_OUTPUT" | grep user | awk '{print $2}')
        SYS_TIME=$(echo "$TIME_OUTPUT" | grep sys | awk '{print $2}')

        printf "Time taken: %s seconds\n" $REAL_TIME
        echo "$FIXED_ROWS,$N_COLS,$REAL_TIME,$USER_TIME,$SYS_TIME" >> $PERFORMANCE_LOG
    fi
    # ========================== MODIFICATION END ==========================
    
    rm -f $TEMP_INPUT $TEMP_OUTPUT
done

printf "\n\033[0;32mProfiling complete. Results are in %s\033[0m\n" $PERFORMANCE_LOG

Rscript /scripts/adjust/plot_performance.R