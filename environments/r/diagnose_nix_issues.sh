#!/bin/bash
# Diagnose why nix-user-chroot is failing on HPC
# Tests three potential issues + actual R package build scenario

set +e  # Don't exit on errors

echo "=========================================="
echo "Nix-User-Chroot Diagnostic Tests"
echo "=========================================="
echo ""

# Test 1: User namespaces
echo "TEST 1: User Namespaces"
echo "----------------------------------------"
if unshare --user --pid --mount --fork true 2>/dev/null; then
    echo "✅ User namespaces are ENABLED"
    USERNS_OK=1
else
    echo "❌ User namespaces are DISABLED or RESTRICTED"
    echo "   Error: $(unshare --user --pid --mount --fork true 2>&1)"
    USERNS_OK=0
fi
echo ""

# Test 2: Filesystem features
echo "TEST 2: Filesystem Features"
echo "----------------------------------------"
TEST_DIR="/tmp/nix-fs-test-$$"
mkdir -p "$TEST_DIR"

# Test chmod
touch "$TEST_DIR/testfile"
if chmod 755 "$TEST_DIR/testfile" 2>/dev/null; then
    echo "✅ chmod works"
    CHMOD_OK=1
else
    echo "❌ chmod fails: $(chmod 755 "$TEST_DIR/testfile" 2>&1)"
    CHMOD_OK=0
fi

# Test chown (will likely fail without root, but let's see)
if chown $(id -u):$(id -g) "$TEST_DIR/testfile" 2>/dev/null; then
    echo "✅ chown works"
else
    echo "⚠️  chown fails (expected without root): $(chown $(id -u):$(id -g) "$TEST_DIR/testfile" 2>&1)"
fi

# Test filesystem type
FS_TYPE=$(df -T "$TEST_DIR" | tail -1 | awk '{print $2}')
echo "   Filesystem type: $FS_TYPE"
if [[ "$FS_TYPE" == "nfs"* ]] || [[ "$FS_TYPE" == "cifs" ]] || [[ "$FS_TYPE" == "lustre" ]]; then
    echo "⚠️  Network/distributed filesystem detected - may cause issues"
    FS_OK=0
else
    echo "✅ Local filesystem"
    FS_OK=1
fi

rm -rf "$TEST_DIR"
echo ""

# Test 3: Permission model inside nix-user-chroot
echo "TEST 3: Permissions Inside nix-user-chroot"
echo "----------------------------------------"
NIX_ROOT="/grphome/grp_batch_effects/nix"
NIX_CHROOT="$NIX_ROOT/nix-user-chroot"

if [ ! -x "$NIX_CHROOT" ]; then
    echo "❌ nix-user-chroot not found at $NIX_CHROOT"
    CHROOT_OK=0
else
    # Test if we can create and chmod files inside the chroot
    TEST_RESULT=$($NIX_CHROOT "$NIX_ROOT" bash -c '
        TEST_FILE="/tmp/chroot-test-$$"
        touch "$TEST_FILE" 2>/dev/null || exit 1
        chmod 755 "$TEST_FILE" 2>/dev/null || exit 2
        rm "$TEST_FILE"
        exit 0
    ' 2>&1)
    
    EXIT_CODE=$?
    
    if [ $EXIT_CODE -eq 0 ]; then
        echo "✅ Permissions work inside nix-user-chroot"
        CHROOT_OK=1
    elif [ $EXIT_CODE -eq 1 ]; then
        echo "❌ Cannot create files inside chroot"
        echo "   Error: $TEST_RESULT"
        CHROOT_OK=0
    elif [ $EXIT_CODE -eq 2 ]; then
        echo "❌ Cannot chmod files inside chroot"
        echo "   Error: $TEST_RESULT"
        CHROOT_OK=0
    else
        echo "❌ Unknown error inside chroot"
        echo "   Error: $TEST_RESULT"
        CHROOT_OK=0
    fi
fi
echo ""

# Test 4: Actual Nix build sandbox test
echo "TEST 4: Nix Build Sandbox"
echo "----------------------------------------"
if [ -x "$NIX_CHROOT" ]; then
    echo "Testing if Nix can unpack and chmod a simple tarball..."
    
    $NIX_CHROOT "$NIX_ROOT" bash -c '
        source ~/.nix-profile/etc/profile.d/nix.sh 2>/dev/null
        
        # Create a test derivation that mimics R package unpacking
        TEST_DRV=$(mktemp -d)
        cd "$TEST_DRV"
        
        # Create a simple tarball
        mkdir test-pkg
        echo "test" > test-pkg/file.txt
        tar czf test.tar.gz test-pkg
        
        # Try to unpack and chmod like Nix does
        tar xzf test.tar.gz
        chmod -R u+w test-pkg 2>/dev/null
        
        if [ $? -eq 0 ]; then
            echo "✅ Tarball unpack + chmod works"
            exit 0
        else
            echo "❌ Tarball unpack + chmod fails"
            exit 1
        fi
    ' 2>&1
    
    if [ $? -eq 0 ]; then
        NIX_UNPACK_OK=1
    else
        NIX_UNPACK_OK=0
    fi
else
    echo "⚠️  Skipping (nix-user-chroot not available)"
    NIX_UNPACK_OK=0
fi
echo ""

# Test 5: R Package Build Scenario
echo "TEST 5: R Package Build Scenario"
echo "----------------------------------------"
if [ -x "$NIX_CHROOT" ]; then
    echo "Testing actual R package build process..."
    
    $NIX_CHROOT "$NIX_ROOT" bash -c '
        source ~/.nix-profile/etc/profile.d/nix.sh 2>/dev/null
        
        cd /tmp
        
        # Check if curl tarball derivation exists
        CURL_DRV="/nix/store/mv0nd3v3maqsiqxvcbvvy0g1pv4698yn-curl_5.2.3.tar.gz.drv"
        if [ -f "$CURL_DRV" ]; then
            echo "Found curl derivation: $CURL_DRV"
            
            # Try to realize it
            echo "Attempting to fetch curl tarball..."
            nix-store -r "$CURL_DRV" 2>&1 | head -10
            
            # Check if tarball is now in store
            TARBALL=$(ls /nix/store/*-curl_5.2.3.tar.gz 2>/dev/null | head -1)
            if [ -n "$TARBALL" ]; then
                echo "✅ Tarball fetched: $TARBALL"
                
                # Try manual unpack
                UNPACK_DIR=$(mktemp -d)
                cd "$UNPACK_DIR"
                echo "Unpacking to: $UNPACK_DIR"
                
                tar xzf "$TARBALL" 2>&1
                UNPACK_EXIT=$?
                echo "Unpack exit code: $UNPACK_EXIT"
                
                if [ $UNPACK_EXIT -eq 0 ] && [ -d "curl" ]; then
                    echo "Attempting chmod on unpacked files..."
                    chmod -R u+w curl 2>&1
                    CHMOD_EXIT=$?
                    echo "chmod exit code: $CHMOD_EXIT"
                    
                    if [ $CHMOD_EXIT -eq 0 ]; then
                        echo "✅ R package unpack + chmod works"
                        exit 0
                    else
                        echo "❌ chmod failed on R package"
                        ls -la curl/ | head -10
                        exit 2
                    fi
                else
                    echo "❌ Unpack failed"
                    exit 3
                fi
            else
                echo "⚠️  Tarball not in store after fetch attempt"
                exit 4
            fi
        else
            echo "⚠️  Curl derivation not found (not built yet)"
            exit 5
        fi
    ' 2>&1
    
    R_BUILD_EXIT=$?
    
    if [ $R_BUILD_EXIT -eq 0 ]; then
        echo "✅ R package build scenario works"
        R_BUILD_OK=1
    elif [ $R_BUILD_EXIT -eq 5 ]; then
        echo "⚠️  Cannot test (derivation not built yet)"
        R_BUILD_OK=2
    else
        echo "❌ R package build scenario fails (exit: $R_BUILD_EXIT)"
        R_BUILD_OK=0
    fi
else
    echo "⚠️  Skipping (nix-user-chroot not available)"
    R_BUILD_OK=2
fi
echo ""

# Test 6: Binary Cache Connectivity
echo "TEST 6: Binary Cache Connectivity"
echo "----------------------------------------"
echo "Testing cache.nixos.org..."
if curl -s -I --max-time 5 "https://cache.nixos.org" | grep -q "200\|301\|302"; then
    echo "✅ cache.nixos.org reachable"
    CACHE_NIXOS_OK=1
else
    echo "❌ cache.nixos.org unreachable"
    CACHE_NIXOS_OK=0
fi

echo "Testing rstats-on-nix.cachix.org..."
if curl -s -I --max-time 5 "https://rstats-on-nix.cachix.org" | grep -q "200\|301\|302"; then
    echo "✅ rstats-on-nix.cachix.org reachable"
    CACHE_RSTATS_OK=1
else
    echo "❌ rstats-on-nix.cachix.org unreachable"
    CACHE_RSTATS_OK=0
fi

# Test if specific package is in cache
echo "Checking if r-curl is in binary cache..."
if curl -s -I --max-time 5 "https://cache.nixos.org/0zdkc97hljp1vn85l4z269qf8b68ica9.narinfo" | grep -q "200"; then
    echo "✅ r-curl found in cache.nixos.org"
    RCURL_CACHE_OK=1
elif curl -s -I --max-time 5 "https://rstats-on-nix.cachix.org/0zdkc97hljp1vn85l4z269qf8b68ica9.narinfo" | grep -q "200"; then
    echo "✅ r-curl found in rstats-on-nix.cachix.org"
    RCURL_CACHE_OK=1
else
    echo "❌ r-curl NOT in any configured cache"
    RCURL_CACHE_OK=0
fi
echo ""

# Test 7: Nix Configuration
echo "TEST 7: Nix Configuration"
echo "----------------------------------------"
if [ -f ~/.config/nix/nix.conf ]; then
    echo "✅ nix.conf exists"
    
    if grep -q "^sandbox = false" ~/.config/nix/nix.conf; then
        echo "✅ sandbox disabled"
        SANDBOX_OK=1
    else
        echo "⚠️  sandbox not explicitly disabled"
        SANDBOX_OK=0
    fi
    
    if grep -q "^builders =" ~/.config/nix/nix.conf; then
        echo "✅ builders set (should prevent local builds)"
        BUILDERS_OK=1
    else
        echo "⚠️  builders not set (may attempt local builds)"
        BUILDERS_OK=0
    fi
    
    if grep -q "substituters.*cache.nixos.org" ~/.config/nix/nix.conf; then
        echo "✅ cache.nixos.org configured"
        SUBST_NIXOS_OK=1
    else
        echo "⚠️  cache.nixos.org not in substituters"
        SUBST_NIXOS_OK=0
    fi
    
    if grep -q "substituters.*rstats-on-nix.cachix.org" ~/.config/nix/nix.conf; then
        echo "✅ rstats-on-nix.cachix.org configured"
        SUBST_RSTATS_OK=1
    else
        echo "⚠️  rstats-on-nix.cachix.org not in substituters"
        SUBST_RSTATS_OK=0
    fi
else
    echo "❌ nix.conf not found"
    SANDBOX_OK=0
    BUILDERS_OK=0
    SUBST_NIXOS_OK=0
    SUBST_RSTATS_OK=0
fi
echo ""

# Test 8: Nix Store Space
echo "TEST 8: Nix Store Space"
echo "----------------------------------------"
if [ -x "$NIX_CHROOT" ]; then
    STORE_INFO=$($NIX_CHROOT "$NIX_ROOT" bash -c 'df -h /nix/store 2>/dev/null' | tail -1)
    if [ -n "$STORE_INFO" ]; then
        echo "Nix store filesystem:"
        echo "$STORE_INFO"
        
        USAGE=$(echo "$STORE_INFO" | awk '{print $5}' | tr -d '%')
        if [ "$USAGE" -lt 90 ]; then
            echo "✅ Sufficient space available"
            SPACE_OK=1
        else
            echo "⚠️  Nix store is ${USAGE}% full"
            SPACE_OK=0
        fi
    else
        echo "⚠️  Could not check store space"
        SPACE_OK=2
    fi
else
    echo "⚠️  Skipping (nix-user-chroot not available)"
    SPACE_OK=2
fi
echo ""

# Test 9: Nix Daemon/Remote
echo "TEST 9: Nix Build Mode"
echo "----------------------------------------"
if [ -x "$NIX_CHROOT" ]; then
    BUILD_MODE=$($NIX_CHROOT "$NIX_ROOT" bash -c '
        source ~/.nix-profile/etc/profile.d/nix.sh 2>/dev/null
        if [ -n "$NIX_REMOTE" ]; then
            echo "NIX_REMOTE=$NIX_REMOTE"
        else
            echo "local"
        fi
    ')
    echo "Build mode: $BUILD_MODE"
    if [[ "$BUILD_MODE" == "local" ]]; then
        echo "✅ Using local builds (expected for nix-user-chroot)"
        BUILD_MODE_OK=1
    else
        echo "⚠️  Non-standard build mode"
        BUILD_MODE_OK=0
    fi
else
    echo "⚠️  Skipping (nix-user-chroot not available)"
    BUILD_MODE_OK=2
fi
echo ""
#!/bin/bash
# ==========================================================
# Nix-User-Chroot Diagnostic Tests (Improved)
# ==========================================================
# This script verifies if local Nix builds are allowed on an HPC.
# It runs inside nix-user-chroot and performs minimal derivation tests.
# ==========================================================

set -euo pipefail

log() { echo -e "\033[1m$1\033[0m"; }
status_ok() { echo -e "✅ $1"; }
status_fail() { echo -e "❌ $1"; }

# --- Detect nix executable ---
NIX_BIN=$(command -v nix || true)
if [ -z "$NIX_BIN" ]; then
    status_fail "nix command not found. Must be run inside nix-user-chroot environment."
    exit 1
fi
status_ok "Found nix executable at $NIX_BIN"

# --- Basic Environment Info ---
log "System Info"
echo "User: $(whoami)"
echo "Host: $(hostname)"
echo "Current dir: $(pwd)"
echo "Nix version: $(nix --version)"
echo

# --- Test 1: User namespaces ---
log "TEST 1: User Namespaces"
if unshare --user true 2>/dev/null; then
    status_ok "User namespaces enabled"
else
    status_fail "User namespaces disabled (nix-user-chroot will fail)"
fi

echo
# --- Test 2: Filesystem features ---
log "TEST 2: Filesystem features"
TMPFILE=$(mktemp)
chmod 600 "$TMPFILE" && status_ok "chmod works" || status_fail "chmod failed"
chown $(whoami) "$TMPFILE" && status_ok "chown works" || status_fail "chown failed"
FSTYPE=$(stat -f -c %T "$TMPFILE")
echo "Filesystem type: $FSTYPE"
rm -f "$TMPFILE"

# --- Test 3: Trivial Nix derivation ---
log "TEST 3: Trivial Nix derivation"
cat > trivial.nix <<'EOF'
with import <nixpkgs> {};
stdenv.mkDerivation {
  name = "trivial-test";
  builder = "/bin/sh";
  args = [ "-c" "echo Hello, Nix! > $out" ];
}
EOF

if nix build -f trivial.nix 2> build.log; then
    status_ok "Local derivation built successfully"
    cat ./result
else
    status_fail "Local build failed"
    echo "--- build.log ---"
    cat build.log
fi

echo
# --- Test 4: Dry-run to detect explicit refusal ---
log "TEST 4: Nix dry-run refusal check"
if nix build -f trivial.nix --dry-run 2>dryrun.log; then
    status_ok "Dry-run completed (no refusal)"
else
    if grep -q "refusing to build" dryrun.log; then
        status_fail "Nix explicitly refuses local builds"
    else
        status_fail "Dry-run failed for another reason"
    fi
    cat dryrun.log
fi

# --- Test 5: Config Summary ---
log "TEST 5: nix.conf summary"
if [ -f "$HOME/.config/nix/nix.conf" ]; then
    grep -E 'builders|sandbox|substituters' "$HOME/.config/nix/nix.conf" || true
else
    echo "No nix.conf found in user config"
fi




# Summary
echo "=========================================="
echo "SUMMARY"
echo "=========================================="
echo ""
echo "User Namespaces:     $([ $USERNS_OK -eq 1 ] && echo '✅ OK' || echo '❌ FAIL')"
echo "Filesystem:          $([ $FS_OK -eq 1 ] && echo '✅ OK' || echo '⚠️  WARNING')"
echo "chmod:               $([ $CHMOD_OK -eq 1 ] && echo '✅ OK' || echo '❌ FAIL')"
echo "Chroot Permissions:  $([ $CHROOT_OK -eq 1 ] && echo '✅ OK' || echo '❌ FAIL')"
echo "Nix Unpack Test:     $([ $NIX_UNPACK_OK -eq 1 ] && echo '✅ OK' || echo '❌ FAIL')"
echo "R Package Build:     $([ $R_BUILD_OK -eq 1 ] && echo '✅ OK' || ([ $R_BUILD_OK -eq 2 ] && echo '⚠️  SKIPPED' || echo '❌ FAIL'))"
echo "Cache Connectivity:  $([ $CACHE_NIXOS_OK -eq 1 ] && [ $CACHE_RSTATS_OK -eq 1 ] && echo '✅ OK' || echo '⚠️  PARTIAL')"
echo "r-curl in Cache:     $([ $RCURL_CACHE_OK -eq 1 ] && echo '✅ YES' || echo '❌ NO')"
echo "Nix Config:          $([ $SANDBOX_OK -eq 1 ] && [ $BUILDERS_OK -eq 1 ] && echo '✅ OK' || echo '⚠️  NEEDS REVIEW')"
echo "Store Space:         $([ $SPACE_OK -eq 1 ] && echo '✅ OK' || ([ $SPACE_OK -eq 2 ] && echo '⚠️  UNKNOWN' || echo '⚠️  LOW'))"
echo ""

if [ $USERNS_OK -eq 0 ]; then
    echo "ROOT CAUSE: User namespaces are disabled"
    echo "  → nix-user-chroot cannot work without user namespaces"
    echo "  → Contact HPC admins to enable: sysctl kernel.unprivileged_userns_clone=1"
elif [ $CHMOD_OK -eq 0 ] || [ $CHROOT_OK -eq 0 ]; then
    echo "ROOT CAUSE: Permission operations fail"
    echo "  → Filesystem or security policies prevent chmod"
    echo "  → This breaks Nix's build process"
elif [ $NIX_UNPACK_OK -eq 0 ]; then
    echo "ROOT CAUSE: Nix build operations fail in chroot"
    echo "  → Even though basic operations work, Nix's build sandbox fails"
    echo "  → Likely due to nested namespace restrictions"
elif [ $R_BUILD_OK -eq 0 ]; then
    echo "ROOT CAUSE: R package builds specifically fail"
    echo "  → Basic operations work, but R package unpack/chmod fails"
    echo "  → This is the exact error you're seeing in builds"
elif [ $RCURL_CACHE_OK -eq 0 ]; then
    echo "ROOT CAUSE: Package not in binary cache"
    echo "  → r-curl (and likely other packages) not available pre-built"
    echo "  → With 'builders =' set, Nix should refuse to build but tries anyway"
    echo "  → Solutions:"
    echo "    1. Build on a working Nix system, copy closure to HPC"
    echo "    2. Try different nixpkgs snapshot with better cache coverage"
    echo "    3. Remove 'builders =' to allow local builds (if permissions work)"
elif [ $CACHE_NIXOS_OK -eq 0 ] || [ $CACHE_RSTATS_OK -eq 0 ]; then
    echo "ROOT CAUSE: Binary cache unreachable"
    echo "  → Network connectivity issues to binary caches"
    echo "  → Check firewall/proxy settings"
elif [ $R_BUILD_OK -eq 2 ]; then
    echo "INCONCLUSIVE: R package test skipped (derivation not available)"
    echo "  → All basic tests pass"
    echo "  → Likely issue: packages not in binary cache"
    echo "  → Solution: Build on working Nix system, copy closure to HPC"
else
    echo "ALL TESTS PASS!"
    echo "  → Environment should work for R package builds"
    echo "  → If builds still fail, issue is with binary cache availability"
    echo "  → Try: nix-build with --keep-going to see which packages fail"
fi
echo ""
