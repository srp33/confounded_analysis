#!/bin/bash

echo "Checking for missing output files..."
echo ""

# Check classify_adjusters outputs
echo "=== Missing classify_adjusters outputs ==="
missing_count=0
for log in grp_batch_effects/outputs/book_chapter/logs/classify_adjusters/*.log; do
    if [ -f "$log" ]; then
        # Extract the base name and construct expected output path
        base=$(basename "$log" .log)
        expected_output="grp_batch_effects/outputs/book_chapter/results/adjusters/individual/${base}.csv"
        
        if [ ! -f "$expected_output" ]; then
            echo "Missing: $expected_output (log: $log)"
            ((missing_count++))
            if [ $missing_count -ge 20 ]; then
                echo "... (showing first 20)"
                break
            fi
        fi
    fi
done

if [ $missing_count -eq 0 ]; then
    echo "None found"
fi

echo ""
echo "Total missing classify_adjusters outputs: $missing_count"
