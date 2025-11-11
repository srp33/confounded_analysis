#!/bin/bash
# Monitor what Nix is doing during a build
# Run this in a separate terminal while the build is running

echo "=== Nix Build Monitor ==="
echo "Press Ctrl+C to stop monitoring"
echo ""

while true; do
  clear
  echo "=== $(date) ==="
  echo ""
  
  # Check if nix processes are running
  echo "--- Nix Processes ---"
  ps aux | grep -E "(nix-shell|nix-build|nix-daemon)" | grep -v grep || echo "No nix processes found"
  echo ""
  
  # Check network activity (downloads)
  echo "--- Network Activity ---"
  netstat -an | grep ESTABLISHED | grep -E ":(80|443)" | wc -l | xargs echo "Active connections:"
  echo ""
  
  # Check Nix store size
  echo "--- Nix Store Size ---"
  du -sh /grphome/grp_batch_effects/nix/store 2>/dev/null || echo "Store not accessible"
  echo ""
  
  # Check CPU and memory
  echo "--- System Resources ---"
  top -bn1 | head -n 5
  echo ""
  
  sleep 5
done
