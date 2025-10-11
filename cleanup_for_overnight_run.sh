#!/bin/bash

echo "=== Cleaning up for overnight run ==="

# Remove all generated results
echo "Removing simulation results..."
rm -rf scripts/evaluations/robustifying/results/*.csv 2>/dev/null || true

echo "Removing 4-study real data results..."
rm -rf scripts/evaluations/robustifying/results_real_4studies/ 2>/dev/null || true

echo "Removing 6-study real data results..."
rm -rf scripts/evaluations/robustifying/results_real_6studies/ 2>/dev/null || true

echo "Removing simulation results directory..."
rm -rf scripts/evaluations/robustifying/results_sim/ 2>/dev/null || true

echo "Removing all generated figures..."
rm -rf scripts/evaluations/robustifying/figures/*.png 2>/dev/null || true
rm -rf scripts/evaluations/robustifying/figures/*.pdf 2>/dev/null || true

echo "Removing debug files..."
rm -f scripts/evaluations/robustifying/debug_gmm_data.R 2>/dev/null || true

echo "Cleaning SLURM output files..."
rm -f grp_batch_effects/slurm_scripts/slurm_outputs/robustifying/*.out 2>/dev/null || true

echo "=== Cleanup complete ==="
echo "Ready for fresh overnight run!"
echo ""
echo "To start the overnight run:"
echo "sbatch grp_batch_effects/slurm_scripts/robustifying_complete_workflow.sh"