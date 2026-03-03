#!/usr/bin/env bash
# test_integration.sh — End-to-end integration tests for bat_convert.
#
# Tests run:
#   1. Build the C++ project (Release).
#   2. Download 1GIA.pdb and 6CRK.pdb from RCSB (skipped if already present).
#   3. Detect chain composition of each PDB.
#   4. For each PDB / each chain:
#       a. test_pdb_roundtrip  — C++ in-memory + I/O round-trip + physical sanity.
#       b. bat_convert          — produce a .bat file.
#       c. validate_bat.py      — Python independent cross-validation of BAT values.
#   5. Chain-filtering sanity: verify that different chains produce different
#      atom counts (proves --chain flag is honoured).
#
# Prerequisites: cmake (≥3.15), curl or wget, python3 (stdlib only, no pip).
#
# Usage:
#   ./test_integration.sh [--no-build] [--no-download]
#
# Options:
#   --no-build      Skip the cmake build step (assumes binaries already built).
#   --no-download   Skip downloading PDB files (assumes they are already in test_pdbs/).
#
# Exit codes:
#   0  all tests passed
#   1  one or more tests failed

set -uo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CPP_DIR="$SCRIPT_DIR/cpp"
BUILD_DIR="$CPP_DIR/build"
PDB_DIR="$SCRIPT_DIR/test_pdbs"
VALIDATE_PY="$SCRIPT_DIR/validate_bat.py"
CMAKE_BIN="${CMAKE_BIN:-/opt/homebrew/bin/cmake}"
PYTHON="${PYTHON:-python3}"

DO_BUILD=1
DO_DOWNLOAD=1
for arg in "$@"; do
    case "$arg" in
        --no-build)    DO_BUILD=0 ;;
        --no-download) DO_DOWNLOAD=0 ;;
        *) echo "Unknown option: $arg" >&2; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Coloured output helpers
# ---------------------------------------------------------------------------

GREEN='\033[0;32m'
RED='\033[0;31m'
CYAN='\033[0;36m'
NC='\033[0m'

section()  { echo; echo -e "${CYAN}=== $* ===${NC}"; }
pass_msg() { echo -e "  ${GREEN}PASS${NC}  $*"; }
fail_msg() { echo -e "  ${RED}FAIL${NC}  $*"; }

FAILURES=()
record_pass() { pass_msg "$*"; }
record_fail() { fail_msg "$*"; FAILURES+=("$*"); }

# ---------------------------------------------------------------------------
# 0. Build
# ---------------------------------------------------------------------------

if [[ $DO_BUILD -eq 1 ]]; then
    section "Build (Release)"
    mkdir -p "$BUILD_DIR"
    "$CMAKE_BIN" -S "$CPP_DIR" -B "$BUILD_DIR" \
        -DCMAKE_BUILD_TYPE=Release \
        --log-level=WARNING \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=OFF \
        > /dev/null
    "$CMAKE_BIN" --build "$BUILD_DIR" --parallel
    record_pass "cmake build"
else
    section "Build (skipped)"
fi

BAT_CONVERT="$BUILD_DIR/bat_convert"
TEST_ROUNDTRIP="$BUILD_DIR/test_pdb_roundtrip"

for bin in "$BAT_CONVERT" "$TEST_ROUNDTRIP"; do
    if [[ ! -x "$bin" ]]; then
        echo "ERROR: binary not found: $bin" >&2
        exit 1
    fi
done

# ---------------------------------------------------------------------------
# 1. Download PDB files
# ---------------------------------------------------------------------------

mkdir -p "$PDB_DIR"

download_pdb() {
    local id="$1"
    local dest="$PDB_DIR/${id}.pdb"
    if [[ -f "$dest" ]]; then
        echo "  $id.pdb already present."
        return 0
    fi
    local url="https://files.rcsb.org/download/${id}.pdb"
    echo "  Downloading $url ..."
    if command -v curl &>/dev/null; then
        curl -fsSL "$url" -o "$dest" || { echo "  ERROR: curl failed for $id" >&2; return 1; }
    elif command -v wget &>/dev/null; then
        wget -q "$url" -O "$dest" || { echo "  ERROR: wget failed for $id" >&2; return 1; }
    else
        echo "  ERROR: neither curl nor wget is available." >&2
        return 1
    fi
}

if [[ $DO_DOWNLOAD -eq 1 ]]; then
    section "Download PDB files"
    download_pdb 1GIA || record_fail "download 1GIA.pdb"
    download_pdb 6CRK || record_fail "download 6CRK.pdb"
fi

for id in 1GIA 6CRK; do
    if [[ ! -f "$PDB_DIR/${id}.pdb" ]]; then
        echo "ERROR: $PDB_DIR/${id}.pdb not found. Run without --no-download." >&2
        exit 1
    fi
done
record_pass "PDB files present"

# ---------------------------------------------------------------------------
# 2. Detect chains present in each PDB (using Python for robustness)
# ---------------------------------------------------------------------------

section "Detecting chain composition"

get_protein_chains() {
    # Returns space-separated sorted unique chain IDs present in ATOM records.
    "$PYTHON" - "$1" <<'PYEOF'
import sys, re
chains = sorted(set(
    line[21]
    for line in open(sys.argv[1])
    if len(line) > 21 and line[:6] == 'ATOM  '
))
print(' '.join(chains))
PYEOF
}

CHAINS_1GIA=$(get_protein_chains "$PDB_DIR/1GIA.pdb")
CHAINS_6CRK=$(get_protein_chains "$PDB_DIR/6CRK.pdb")
echo "  1GIA ATOM chain IDs: '$CHAINS_1GIA'"
echo "  6CRK ATOM chain IDs: '$CHAINS_6CRK'"

# Convert to arrays
read -r -a CHAINS_1GIA_ARR <<< "$CHAINS_1GIA"
read -r -a CHAINS_6CRK_ARR <<< "$CHAINS_6CRK"

if [[ ${#CHAINS_1GIA_ARR[@]} -eq 0 ]]; then
    record_fail "No chains found in 1GIA.pdb — is the file valid?"
fi
if [[ ${#CHAINS_6CRK_ARR[@]} -eq 0 ]]; then
    record_fail "No chains found in 6CRK.pdb — is the file valid?"
fi

echo "  1GIA has ${#CHAINS_1GIA_ARR[@]} chain(s)."
echo "  6CRK has ${#CHAINS_6CRK_ARR[@]} chain(s)."

if [[ ${#CHAINS_1GIA_ARR[@]} -eq 1 ]]; then
    echo "  1GIA confirmed as single-chain protein."
else
    echo "  NOTE: 1GIA has multiple chains — tests will run per chain."
fi
if [[ ${#CHAINS_6CRK_ARR[@]} -lt 2 ]]; then
    echo "  NOTE: 6CRK has fewer than 2 chains — chain-filter test will be limited."
fi

# ---------------------------------------------------------------------------
# Helper: run one full test suite for a given PDB + chain
# ---------------------------------------------------------------------------

run_tests_for_chain() {
    local pdb_path="$1"
    local pdb_id="$2"
    local chain="$3"
    local bat_out="$PDB_DIR/${pdb_id}_chain${chain}.bat"

    echo
    echo "  --- ${pdb_id}  chain ${chain} ---"

    # (a) C++ round-trip test
    if "$TEST_ROUNDTRIP" "$pdb_path" "$chain" 2>&1 | \
        sed 's/^/    /'; then
        record_pass "${pdb_id}:${chain}  C++ round-trip"
    else
        record_fail "${pdb_id}:${chain}  C++ round-trip"
        # Don't abort; still try the other tests.
    fi

    # (b) bat_convert → produce .bat file
    if "$BAT_CONVERT" "$pdb_path" "$bat_out" --chain "$chain" 2>&1 | \
        sed 's/^/    /'; then
        if [[ -f "$bat_out" ]]; then
            record_pass "${pdb_id}:${chain}  bat_convert"
        else
            record_fail "${pdb_id}:${chain}  bat_convert (output file missing)"
            return
        fi
    else
        record_fail "${pdb_id}:${chain}  bat_convert"
        return
    fi

    # (c) Python independent cross-validation
    if "$PYTHON" "$VALIDATE_PY" "$pdb_path" "$bat_out" "$chain" 2>&1 | \
        sed 's/^/    /'; then
        record_pass "${pdb_id}:${chain}  Python validate_bat"
    else
        record_fail "${pdb_id}:${chain}  Python validate_bat"
    fi
}

# ---------------------------------------------------------------------------
# 3. Tests for 1GIA
# ---------------------------------------------------------------------------

section "1GIA tests"
for ch in "${CHAINS_1GIA_ARR[@]}"; do
    run_tests_for_chain "$PDB_DIR/1GIA.pdb" "1GIA" "$ch"
done

# ---------------------------------------------------------------------------
# 4. Tests for 6CRK (all chains)
# ---------------------------------------------------------------------------

section "6CRK tests"
for ch in "${CHAINS_6CRK_ARR[@]}"; do
    run_tests_for_chain "$PDB_DIR/6CRK.pdb" "6CRK" "$ch"
done

# ---------------------------------------------------------------------------
# 5. Chain-filtering sanity: different chains → different NATOMS counts
# ---------------------------------------------------------------------------

section "Chain-filtering sanity (6CRK)"

if [[ ${#CHAINS_6CRK_ARR[@]} -ge 2 ]]; then
    ch_a="${CHAINS_6CRK_ARR[0]}"
    ch_b="${CHAINS_6CRK_ARR[1]}"
    bat_a="$PDB_DIR/6CRK_chain${ch_a}.bat"
    bat_b="$PDB_DIR/6CRK_chain${ch_b}.bat"

    get_natoms() { grep '^NATOMS' "$1" | awk '{print $2}'; }

    if [[ -f "$bat_a" && -f "$bat_b" ]]; then
        n_a=$(get_natoms "$bat_a")
        n_b=$(get_natoms "$bat_b")
        echo "  Chain $ch_a: $n_a atoms"
        echo "  Chain $ch_b: $n_b atoms"

        if [[ "$n_a" -gt 0 && "$n_b" -gt 0 ]]; then
            record_pass "Both chains produced non-empty output"
        else
            record_fail "One or both chains produced empty output"
        fi

        if [[ "$n_a" -ne "$n_b" ]]; then
            record_pass "Chain $ch_a ($n_a atoms) ≠ chain $ch_b ($n_b atoms): filter is working"
        else
            echo "  NOTE: Both chains have the same atom count ($n_a). This could be"
            echo "        correct if both chains are structurally identical; not a failure."
        fi

        # Verify that a .bat file produced with --chain A contains no chain-B
        # atoms in its ANCHOR section (spot-check via atom index bounds).
        # We check that NATOMS from each per-chain .bat is strictly less than
        # the all-chain total (if such a .bat exists).
        bat_all="$PDB_DIR/6CRK_all.bat"
        if "$BAT_CONVERT" "$PDB_DIR/6CRK.pdb" "$bat_all" 2>&1 | sed 's/^/    /'; then
            n_all=$(get_natoms "$bat_all")
            echo "  All chains combined: $n_all atoms"
            if [[ "$n_a" -lt "$n_all" && "$n_b" -lt "$n_all" ]]; then
                record_pass "Per-chain counts ($n_a, $n_b) < all-chain count ($n_all)"
            else
                record_fail "Per-chain counts are not smaller than all-chain count"
            fi
        fi
    else
        echo "  Skipped (one or both .bat files missing — earlier test may have failed)."
    fi
else
    echo "  Skipped (6CRK has fewer than 2 chains)."
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

section "Summary"
echo "  Total failures: ${#FAILURES[@]}"
if [[ ${#FAILURES[@]} -eq 0 ]]; then
    echo -e "  ${GREEN}ALL TESTS PASSED${NC}"
    exit 0
else
    echo -e "  ${RED}FAILED TESTS:${NC}"
    for f in "${FAILURES[@]}"; do
        echo "    - $f"
    done
    exit 1
fi
