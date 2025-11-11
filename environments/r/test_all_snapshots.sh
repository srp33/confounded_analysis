#!/bin/bash
# Test multiple nixpkgs snapshots to find one with cache coverage
# This tests if we can even install rix itself (prerequisite for everything)

set -euo pipefail

NIX_ROOT="/grphome/grp_batch_effects/nix"
NIX_CHROOT_CMD="$NIX_ROOT/nix-user-chroot $NIX_ROOT"

# Priority dates to test (major releases and recent snapshots)
DATES_TO_TEST=(
    "2024-10-01"  # Last major release before 2024-12-14
    # "2024-08-19"  # Previous major release
    # "2024-06-14"  # Previous major release
    # "2024-04-29"  # Previous major release
    # "2024-02-29"  # Previous major release
    # "2025-11-10"  # Today
    # "2025-11-03"  # Recent
    # "2025-11-01"  # Recent major release
    # "2025-10-27"  # Recent
    # "2025-10-20"  # Recent
    # Others
    # "2019-03-14"
    # "2019-05-05"
    # "2019-07-22"
    # "2019-12-19"
    # "2020-03-12"
    # "2020-04-27"
    # "2020-06-22"
    # "2020-08-20"
    # "2020-10-30"
    # "2021-02-26"
    # "2021-04-01"
    # "2021-05-29"
    # "2021-08-03"
    # "2021-10-28"
    # "2022-01-16"
    # "2022-04-19"
    # "2022-06-22"
    # "2022-08-22"
    # "2022-10-20"
    # "2022-12-20"
    # "2023-02-13"
    # "2023-04-01"
    # "2023-06-15"
    # "2023-08-15"
    # "2023-10-30"
    # "2023-12-30"
    # "2024-02-29"
    # "2024-04-29"
    # "2024-06-14"
    # "2024-08-19"
    # "2024-10-01"
    # "2024-12-14"
    # "2025-01-14"
    # "2025-01-24"
    # "2025-01-27"
    # "2025-02-03"
    # "2025-02-10"
    # "2025-02-17"
    # "2025-02-24"
    # "2025-02-28"
    # "2025-03-03"
    # "2025-03-10"
    # "2025-03-17"
    # "2025-03-24"
    # "2025-03-31"
    # "2025-04-07"
    # "2025-04-11"
    # "2025-04-14"
    # "2025-04-16"
    # "2025-04-29"
    # "2025-05-05"
    # "2025-05-16"
    # "2025-05-19"
    # "2025-05-26"
    # "2025-06-02"
    # "2025-06-09"
    # "2025-06-13"
    # "2025-06-24"
    # "2025-07-02"
    # "2025-07-07"
    # "2025-07-14"
    # "2025-07-21"
    # "2025-07-28"
    # "2025-08-04"
    # "2025-08-11"
    # "2025-08-18"
    # "2025-08-25"
    # "2025-09-01"
    # "2025-09-04"
    # "2025-09-09"
    # "2025-09-11"
    # "2025-09-16"
    # "2025-09-22"
    # "2025-09-29"
    # "2025-10-07"
    # "2025-10-14"
    # "2025-10-20"
    # "2025-10-27"
    # "2025-11-01"
    # "2025-11-03"
    # "2025-11-10"
)

echo "=========================================="
echo "Testing nixpkgs Snapshots for Cache Coverage"
echo "=========================================="
echo ""
echo "This will test if we can install rix (prerequisite for everything)"
echo "Testing ${#DATES_TO_TEST[@]} snapshots..."
echo ""

RESULTS_FILE="snapshot_test_results.txt"
> "$RESULTS_FILE"  # Clear file

for date in "${DATES_TO_TEST[@]}"; do
    {
    echo "----------------------------------------"
    echo "Testing: $date"
    echo "----------------------------------------"
    
    # Try to install rix with this snapshot
    echo "Attempting: nix-shell -p R rPackages.rix (snapshot: $date)"
    
    if $NIX_CHROOT_CMD bash -c "
        source ~/.nix-profile/etc/profile.d/nix.sh
        nix-shell -I nixpkgs=https://github.com/rstats-on-nix/nixpkgs/archive/refs/heads/$date.tar.gz \
                  -p R rPackages.rix \
                  --run 'R --version' 2>&1
    " 2>&1 | tee -a "$RESULTS_FILE" | grep -q "R version"; then
        echo "✅ SUCCESS: $date works!" | tee -a "$RESULTS_FILE"
        echo ""
        echo "=========================================="
        echo "FOUND WORKING SNAPSHOT: $date"
        echo "=========================================="
        echo ""
        echo "You can use this date in your generate_env.R:"
        echo "  rix(date = \"$date\", ...)"
        echo ""
        exit 0
    else
        echo "❌ FAILED: $date (other error)" | tee -a "$RESULTS_FILE"
    fi
    
    echo ""
    } &
done
wait

/grphome/grp_batch_effects/nix/nix-user-chroot /grphome/grp_batch_effects/nix bash -c "
        source ~/.nix-profile/etc/profile.d/nix.sh
        nix-shell -I nixpkgs=https://github.com/rstats-on-nix/nixpkgs/archive/refs/heads/2025-11-10.tar.gz \
                  -p R rPackages.rix \
                  --run 'R --version' 2>&1
        "

echo "=========================================="
echo "CONCLUSION: No Working Snapshot Found"
echo "=========================================="
echo ""
echo "Tested ${#DATES_TO_TEST[@]} snapshots - none have cache coverage in nix-user-chroot"
echo ""
echo "This means:"
echo "  ❌ rstats-on-nix cache doesn't have pre-built binaries for nix-user-chroot"
echo "  ❌ All snapshots try to build from source"
echo "  ❌ Source builds fail with permission errors"
echo "  ❌ R migration to Nix is not viable in this environment"
echo ""
echo "RECOMMENDATION: Hybrid Approach"
echo "  ✅ Python: Use uv (already 100% working)"
echo "  ✅ R: Keep Apptainer (proven, stable)"
echo ""
echo "See environments/r/NEXT_STEPS.md for details"
echo ""
echo "Full results saved to: $RESULTS_FILE"
