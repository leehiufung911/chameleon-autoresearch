#!/usr/bin/env bash
# spawn_agents.sh — launch N parallel research agents
#
# Fires off N background instances of research_loop.sh, each tagged with a
# unique AGENT_ID (a1, a2, ...), each running ITERS_PER_AGENT iterations.
# A single compute lock (scripts/run_locked.sh) serializes GPU/CPU use.
#
# Usage:
#   bash scripts/spawn_agents.sh [N_AGENTS] [ITERS_PER_AGENT] [COOLDOWN] [BACKEND]
#
# Examples:
#   bash scripts/spawn_agents.sh                     # 3 agents, 5 iters each, opencode
#   bash scripts/spawn_agents.sh 3 5                 # 3 agents, 5 iters each, opencode
#   bash scripts/spawn_agents.sh 3 10 90 gemini      # 3 agents, 10 iters, gemini
#
# Output:
#   logs/agent_<id>.log  — per-agent stdout/stderr
#   .agents.pids         — PIDs of spawned agents (one per line: "id pid")
#
# Agents run independently. Exit when all agents finish, or kill via:
#   bash scripts/stop_agents.sh

set -uo pipefail

N_AGENTS="${1:-3}"
ITERS_PER_AGENT="${2:-5}"
COOLDOWN="${3:-90}"
BACKEND="${4:-opencode}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
LOG_DIR="$PROJECT_DIR/logs"
PIDS_FILE="$PROJECT_DIR/.agents.pids"

mkdir -p "$LOG_DIR" "$PROJECT_DIR/claims/active" "$PROJECT_DIR/claims/done"
cd "$PROJECT_DIR"

# Fresh pids file (append-mode below)
: > "$PIDS_FILE"

echo "=== Spawning $N_AGENTS parallel research agents ==="
echo "  Backend:          $BACKEND"
echo "  Iters per agent:  $ITERS_PER_AGENT  (total: $(( N_AGENTS * ITERS_PER_AGENT )))"
echo "  Cooldown:         ${COOLDOWN}s"
echo "  Per-agent logs:   $LOG_DIR/agent_<id>.log"
echo "  Compute lock:     /tmp/chameleon-compute.lock (one experiment at a time)"
echo ""

for idx in $(seq 1 "$N_AGENTS"); do
    AGENT_ID="a${idx}"
    AGENT_LOG="$LOG_DIR/agent_${AGENT_ID}.log"

    echo "  -> Launching agent $AGENT_ID  (log: $AGENT_LOG)"

    # nohup so the agent survives parent shell exit. Each agent gets a
    # small pre-launch stagger so opencode's rate limiter doesn't see 3
    # simultaneous auth/first-prompt calls.
    (
        sleep $(( (idx - 1) * 10 ))
        nohup bash research_loop.sh "$ITERS_PER_AGENT" "$COOLDOWN" "$BACKEND" "$AGENT_ID" \
            >> "$AGENT_LOG" 2>&1 &
        echo "$AGENT_ID $!" >> "$PIDS_FILE"
    ) &
done

wait  # wait only for the launcher subshells, not the agents themselves

echo ""
echo "=== All agents launched ==="
cat "$PIDS_FILE"
echo ""
echo "Tail any agent:     tail -f $LOG_DIR/agent_<id>.log"
echo "Stop all agents:    bash scripts/stop_agents.sh"
echo "Current claims:     ls claims/active/"
