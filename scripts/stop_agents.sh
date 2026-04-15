#!/usr/bin/env bash
# stop_agents.sh — stop all parallel research agents spawned via spawn_agents.sh
#
# Reads .agents.pids and kills each PID, then cleans up the pids file.
# Also drops the compute lock if stale.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PIDS_FILE="$PROJECT_DIR/.agents.pids"
LOCK_DIR="${LOCK_DIR:-/tmp/chameleon-compute.lock}"

if [ ! -f "$PIDS_FILE" ]; then
    echo "No .agents.pids file found. Nothing to stop."
    exit 0
fi

echo "=== Stopping parallel research agents ==="
while read -r id pid; do
    [ -z "$pid" ] && continue
    if kill -0 "$pid" 2>/dev/null; then
        echo "  -> killing agent $id (pid $pid)"
        kill "$pid" 2>/dev/null || true
        # Also kill any direct children (opencode/python that the agent spawned)
        # Best-effort; on Git Bash pkill -P is flaky, so try taskkill too.
        taskkill //PID "$pid" //T //F >/dev/null 2>&1 || true
    else
        echo "  -> agent $id (pid $pid) already gone"
    fi
done < "$PIDS_FILE"

rm -f "$PIDS_FILE"

# Clean stale compute lock if nobody's holding it
if [ -d "$LOCK_DIR" ]; then
    echo "  -> removing stale compute lock $LOCK_DIR"
    rm -rf "$LOCK_DIR" 2>/dev/null || true
fi

echo "=== Done ==="
