#!/usr/bin/env bash
# start_session.sh — one-shot launcher for a full research session
#
# Kicks off:
#   1. N=3 parallel research agents (via scripts/spawn_agents.sh)
#   2. The Sonnet 4.6 monitor (via monitor/run_monitor.sh)
#
# The agents run ITERS_PER_AGENT iterations each, then exit. The monitor
# runs every 30 min until this script (or the user) drops monitor/STOP.
#
# Backend-flip logic already lives inside research_loop.sh — each agent
# independently falls back gemini -> gemini-pro -> opencode on quota errors.
#
# Usage:
#   bash start_session.sh [N_AGENTS] [ITERS_PER_AGENT] [BACKEND] [COOLDOWN]
#
# Examples:
#   bash start_session.sh                     # 3 agents, 5 iters each, opencode
#   bash start_session.sh 3 8 gemini          # 3 agents, 8 iters, start on gemini
#   bash start_session.sh 1 5 opencode        # solo-agent mode (no swarm)

set -uo pipefail

N_AGENTS="${1:-3}"
ITERS_PER_AGENT="${2:-5}"
BACKEND="${3:-opencode}"
COOLDOWN="${4:-90}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "============================================================"
echo "  Chameleon Autoresearch Session"
echo "============================================================"
echo "  Agents:     $N_AGENTS  (iters each: $ITERS_PER_AGENT)"
echo "  Backend:    $BACKEND (auto-fallback on quota errors)"
echo "  Monitor:    Sonnet 4.6 (every 30 min, writes monitor/logs/)"
echo "  Stop:       touch monitor/STOP     (stops the monitor)"
echo "              bash scripts/stop_agents.sh   (kills all agents)"
echo "============================================================"
echo ""

# Ensure no stale STOP file
rm -f monitor/STOP

# 1. Launch parallel agents (backgrounded, detached via nohup)
bash scripts/spawn_agents.sh "$N_AGENTS" "$ITERS_PER_AGENT" "$COOLDOWN" "$BACKEND"
echo ""

# 2. Launch monitor (backgrounded, detached via nohup)
echo "=== Launching Sonnet monitor (every 30 min) ==="
mkdir -p monitor/logs
nohup bash monitor/run_monitor.sh 1800 >> monitor/monitor_runner.log 2>&1 &
MONITOR_PID=$!
echo "$MONITOR_PID" > monitor/.monitor.pid
echo "  Monitor PID: $MONITOR_PID"
echo "  Monitor runner log: monitor/monitor_runner.log"
echo ""

echo "============================================================"
echo "  Session launched. Detach this shell at will."
echo "  View progress:"
echo "    git log --oneline -10"
echo "    tail -f logs/agent_a1.log"
echo "    ls -lt monitor/logs/check_*.md | head -5"
echo "============================================================"
