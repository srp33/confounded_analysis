#!/bin/bash

echo "Analyzing failed jobs..."
echo ""

failed_logs=()
for log in grp_batch_effects/outputs/book_chapter/logs/classify_adjusters/*.log; do
    if [ -f "$log" ]; then
        base=$(basename "$log" .log)
        expected_output="grp_batch_effects/outputs/book_chapter/results/adjusters/individual/${base}.csv"
        
        if [ ! -f "$expected_output" ]; then
            failed_logs+=("$log")
        fi
    fi
done

echo "Total failed jobs: ${#failed_logs[@]}"
echo ""

if [ ${#failed_logs[@]} -gt 0 ]; then
    echo "Sample error from first failed log:"
    echo "===================================="
    tail -5 "${failed_logs[0]}"
    echo ""
    echo "Failed log: ${failed_logs[0]}"
fi
