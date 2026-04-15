#!/usr/bin/env bash
# run_monitor.sh — Standalone monitoring agent for the research loop
#
# Runs Claude Code headless every INTERVAL seconds to check loop health.
# Each check writes a log file to monitor/logs/.
#
# Usage:
#   bash monitor/run_monitor.sh [INTERVAL_SECONDS]
#   bash monitor/run_monitor.sh 1800    # every 30 min (default)
#   bash monitor/run_monitor.sh 900     # every 15 min
#
# To run in the background:
#   nohup bash monitor/run_monitor.sh 1800 >> monitor/monitor_runner.log 2>&1 &

set -uo pipefail

INTERVAL="${1:-1800}"  # Default: 30 minutes
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
MONITOR_DIR="$SCRIPT_DIR"
LOG_DIR="$MONITOR_DIR/logs"

mkdir -p "$LOG_DIR"

echo "=== Research Loop Monitor ==="
echo "  Interval:  ${INTERVAL}s ($(( INTERVAL / 60 )) min)"
echo "  Project:   $PROJECT_DIR"
echo "  Logs:      $LOG_DIR"
echo "  Started:   $(date)"
echo ""

CHECK_NUM=0

while true; do
    CHECK_NUM=$((CHECK_NUM + 1))
    TIMESTAMP=$(date "+%Y-%m-%d_%H-%M")
    LOG_FILE="$LOG_DIR/check_${TIMESTAMP}.md"

    echo "[$(date)] Check #${CHECK_NUM} starting..."

    # Claude Code headless: cd to monitor dir so CLAUDE.md is discovered,
    # and use --add-dir to grant access to parent project for reading logs
    cd "$MONITOR_DIR"
    claude -p "Run your monitoring checks now. Write findings to logs/check_${TIMESTAMP}.md. Be concise." \
        --model claude-sonnet-4-6 \
        --allowedTools "Bash,Read,Write,Glob,Grep" \
        --add-dir "$PROJECT_DIR" \
        2>&1 | tee "$LOG_DIR/check_${TIMESTAMP}_raw.log"

    EXIT_CODE=${PIPESTATUS[0]}

    if [ -f "$LOG_FILE" ]; then
        echo "[$(date)] Check #${CHECK_NUM} complete -> $LOG_FILE"
        # Print a brief summary to stdout
        grep -E "^\*\*Status\*\*:" "$LOG_FILE" 2>/dev/null || true
    else
        echo "[$(date)] Check #${CHECK_NUM} finished (exit=$EXIT_CODE) but no log file written"
    fi

    echo "[$(date)] Next check in ${INTERVAL}s..."
    echo ""
    sleep "$INTERVAL"
done
