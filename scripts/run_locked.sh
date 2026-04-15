#!/usr/bin/env bash
# run_locked.sh — global compute lock for experiment scripts
#
# Wraps a command with an mkdir-based lock so only ONE experiment runs
# on the GPU/CPU at a time across all parallel research agents.
#
# mkdir is atomic on both POSIX and NTFS, so this works on Git Bash for
# Windows where `flock` is not available.
#
# Usage:
#   bash scripts/run_locked.sh <command> [args...]
#   bash scripts/run_locked.sh python experiments/iter_5.py
#
# Env vars:
#   LOCK_DIR       — lock directory path (default: /tmp/chameleon-compute.lock)
#   LOCK_POLL_SEC  — seconds between acquire attempts (default: 5)
#   LOCK_TIMEOUT   — max seconds to wait (default: 7200 = 2 hours)
#   LOCK_VERBOSE   — if set, print wait status every minute

set -uo pipefail

LOCK_DIR="${LOCK_DIR:-/tmp/chameleon-compute.lock}"
LOCK_POLL_SEC="${LOCK_POLL_SEC:-5}"
LOCK_TIMEOUT="${LOCK_TIMEOUT:-7200}"
LOCK_VERBOSE="${LOCK_VERBOSE:-}"

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <command> [args...]" >&2
    exit 2
fi

AGENT_TAG="${AGENT_ID:-unknown}"
WAITED=0
FIRST_WAIT=1

while ! mkdir "$LOCK_DIR" 2>/dev/null; do
    if [ "$FIRST_WAIT" -eq 1 ]; then
        HOLDER="unknown"
        if [ -f "$LOCK_DIR/owner" ]; then
            HOLDER=$(cat "$LOCK_DIR/owner" 2>/dev/null || echo unknown)
        fi
        echo "[$AGENT_TAG] waiting for compute lock (held by: $HOLDER)..." >&2
        FIRST_WAIT=0
    fi
    sleep "$LOCK_POLL_SEC"
    WAITED=$((WAITED + LOCK_POLL_SEC))

    if [ -n "$LOCK_VERBOSE" ] && [ $((WAITED % 60)) -lt "$LOCK_POLL_SEC" ]; then
        echo "[$AGENT_TAG] still waiting ($WAITED s)..." >&2
    fi

    if [ "$WAITED" -ge "$LOCK_TIMEOUT" ]; then
        echo "[$AGENT_TAG] LOCK TIMEOUT after ${WAITED}s. Giving up." >&2
        exit 124
    fi
done

# Lock acquired. Record owner and ensure release on exit.
echo "$AGENT_TAG (pid $$) at $(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$LOCK_DIR/owner"
trap 'rm -rf "$LOCK_DIR"' EXIT INT TERM

if [ -n "$LOCK_VERBOSE" ] || [ "$FIRST_WAIT" -eq 0 ]; then
    echo "[$AGENT_TAG] acquired compute lock." >&2
fi

# Run the command, preserving its exit code.
"$@"
RC=$?
exit $RC
