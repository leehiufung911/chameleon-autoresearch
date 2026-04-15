#!/usr/bin/env bash
# monitor_loop.sh — Health check for the chameleon research loop.
# Called by the monitoring agent to assess whether the loop is healthy.
# Outputs structured status for the monitoring agent to interpret.

set -euo pipefail

PROJECT_DIR="C:/Users/mic23/prototype-embed/chameleon-research-loop"
LOG_DIR="$PROJECT_DIR/logs"
EXP_DIR="$PROJECT_DIR/experiments"

echo "=== LOOP HEALTH CHECK $(date) ==="
echo ""

# 1. Process check: is opencode or a python compute job running?
echo "--- PROCESSES ---"
OPENCODE_PROCS=$(tasklist 2>/dev/null | grep -i opencode | wc -l || echo 0)
PYTHON_PROCS=$(tasklist 2>/dev/null | grep -i python | wc -l || echo 0)
echo "  opencode processes: $OPENCODE_PROCS"
echo "  python processes:   $PYTHON_PROCS"

# 2. Latest log file and its age
echo ""
echo "--- LATEST LOG ---"
LATEST_LOG=$(ls -t "$LOG_DIR"/iter_*.log 2>/dev/null | head -1)
if [ -n "$LATEST_LOG" ]; then
    LOG_SIZE=$(wc -c < "$LATEST_LOG")
    LOG_MOD=$(stat -c %Y "$LATEST_LOG" 2>/dev/null || stat -f %m "$LATEST_LOG" 2>/dev/null || echo 0)
    NOW=$(date +%s)
    LOG_AGE_MIN=$(( (NOW - LOG_MOD) / 60 ))
    echo "  Latest: $(basename "$LATEST_LOG")"
    echo "  Size:   $LOG_SIZE bytes"
    echo "  Age:    ${LOG_AGE_MIN} minutes since last modification"

    # Check last few lines for errors
    LAST_LINES=$(tail -5 "$LATEST_LOG" 2>/dev/null | tr -d '\000')
    echo "  Last 3 lines:"
    tail -3 "$LATEST_LOG" 2>/dev/null | sed 's/^/    /' | head -3
else
    echo "  No log files found"
fi

# 3. Latest iteration number from git
echo ""
echo "--- GIT STATUS ---"
cd "$PROJECT_DIR"
LATEST_COMMIT=$(git log --oneline -1 2>/dev/null || echo "no commits")
echo "  Latest commit: $LATEST_COMMIT"
UNCOMMITTED=$(git status --short 2>/dev/null | wc -l)
echo "  Uncommitted changes: $UNCOMMITTED files"

# 4. Latest output files
echo ""
echo "--- LATEST OUTPUT FILES ---"
LATEST_OUTPUT=$(ls -t "$EXP_DIR"/iter_*_output.txt 2>/dev/null | head -1)
if [ -n "$LATEST_OUTPUT" ]; then
    OUT_SIZE=$(wc -c < "$LATEST_OUTPUT")
    OUT_LINES=$(wc -l < "$LATEST_OUTPUT")
    echo "  Latest: $(basename "$LATEST_OUTPUT") ($OUT_LINES lines, $OUT_SIZE bytes)"
    if [ "$OUT_SIZE" -eq 0 ]; then
        echo "  WARNING: Output file is EMPTY (stdout capture failed)"
    fi
else
    echo "  No output files found"
fi

# 5. CPU usage snapshot (2-second sample)
echo ""
echo "--- RESOURCE USAGE (2s sample) ---"
# Get CPU usage via typeperf (Windows) or top
if command -v typeperf &>/dev/null; then
    CPU_LINE=$(typeperf "\\Processor(_Total)\\% Processor Time" -sc 2 2>/dev/null | tail -2 | head -1)
    echo "  CPU: $CPU_LINE"
elif command -v wmic &>/dev/null; then
    CPU=$(wmic cpu get loadpercentage 2>/dev/null | grep -o '[0-9]*' | head -1)
    echo "  CPU: ${CPU:-unknown}%"
else
    echo "  CPU: (no monitoring tool available)"
fi

# 6. GPU check
if command -v nvidia-smi &>/dev/null; then
    GPU_UTIL=$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits 2>/dev/null || echo "N/A")
    GPU_MEM=$(nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader,nounits 2>/dev/null || echo "N/A")
    echo "  GPU util: ${GPU_UTIL}%"
    echo "  GPU mem:  ${GPU_MEM} MiB"
fi

# 7. Summary assessment
# NOTE: Be conservative with alerts. Many "idle" states are normal:
# - Zero CPU/GPU: model is thinking (cloud inference), not hung
# - Long iterations (30-60+ min): compute jobs legitimately take time
# - No opencode process: could be between iterations or cooldown
# Only flag genuinely abnormal situations.
echo ""
echo "--- ASSESSMENT ---"

# Check for any recent file changes in experiments/ (more reliable than process checks)
LATEST_ANY=$(ls -t "$EXP_DIR"/* "$LOG_DIR"/* 2>/dev/null | head -1)
if [ -n "$LATEST_ANY" ]; then
    ANY_MOD=$(stat -c %Y "$LATEST_ANY" 2>/dev/null || stat -f %m "$LATEST_ANY" 2>/dev/null || echo 0)
    NOW_TS=$(date +%s)
    ANY_AGE_MIN=$(( (NOW_TS - ANY_MOD) / 60 ))
    echo "  Most recent file change: $(basename "$LATEST_ANY") (${ANY_AGE_MIN} min ago)"
fi

if [ "$OPENCODE_PROCS" -gt 0 ]; then
    echo "  Loop status: RUNNING (opencode active)"
    if [ -n "$LATEST_LOG" ] && [ "$LOG_AGE_MIN" -gt 120 ]; then
        echo "  ALERT: Log file unchanged for ${LOG_AGE_MIN} minutes (>2 hours)."
        echo "  AND opencode still running. Likely hung — investigate."
    elif [ -n "$LATEST_LOG" ] && [ "$LOG_AGE_MIN" -gt 60 ]; then
        echo "  NOTE: Log file last modified ${LOG_AGE_MIN} min ago."
        echo "  Could be long compute, or model reasoning. Not necessarily a problem."
    else
        echo "  Normal: opencode active, log recently updated."
    fi
elif [ "$PYTHON_PROCS" -gt 0 ]; then
    echo "  Loop status: COMPUTING (python running, no opencode)"
    echo "  A compute job is running. This is normal."
else
    echo "  Loop status: IDLE (no opencode or python processes)"
    if [ -n "$LATEST_LOG" ] && [ "$LOG_AGE_MIN" -lt 10 ]; then
        echo "  Just finished or in cooldown. Normal."
    elif [ -n "$LATEST_LOG" ] && [ "$LOG_AGE_MIN" -gt 120 ]; then
        echo "  ALERT: No processes running and no activity for ${LOG_AGE_MIN} min."
        echo "  Loop likely completed or crashed. Check last log for errors."
    else
        echo "  Idle. Loop may be between batches or completed."
    fi
fi

echo ""
echo "=== END HEALTH CHECK ==="
