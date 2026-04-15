#!/usr/bin/env bash
# research_loop.sh — autonomous chameleonicity research loop (Phase 2)
#
# Runs an AI coding agent in headless mode, one iteration per invocation.
# Supports multiple backends: gemini CLI, gemini-pro, opencode + NIM.
# Auto-falls-back on quota exhaustion.
#
# Usage:
#   bash research_loop.sh [MAX_ITERATIONS] [COOLDOWN_SECONDS] [BACKEND] [AGENT_ID]
#
# Examples:
#   bash research_loop.sh 5                   # 5 iterations, gemini-flash, solo
#   bash research_loop.sh 10 60               # 10 iterations, 60s cooldown
#   bash research_loop.sh 5 90 opencode       # use opencode + NIM directly
#   bash research_loop.sh 5 90 opencode a1    # as agent "a1" in a parallel swarm
#
# When AGENT_ID is set (non-"solo"), per-iteration artifacts are tagged with the
# ID (iter_N_<id>.log etc.) and the agent is instructed to wrap its compute via
# scripts/run_locked.sh so only one experiment runs at a time across the swarm.

set -euo pipefail

MAX_ITER="${1:-5}"
COOLDOWN="${2:-30}"
BACKEND="${3:-gemini}"  # gemini | gemini-pro | opencode
AGENT_ID="${4:-solo}"   # solo | a1 | a2 | a3 | ...
PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
LOG_DIR="$PROJECT_DIR/logs"
EXPERIMENTS_DIR="$PROJECT_DIR/experiments"
CLAIMS_DIR="$PROJECT_DIR/claims"
PYTHON_EXE="C:/Users/mic23/miniconda3/envs/chameleon/python.exe"

# File naming: solo uses bare iter_N.log; swarm agents tag their files.
if [ "$AGENT_ID" = "solo" ]; then
    TAG=""
else
    TAG="_${AGENT_ID}"
fi

mkdir -p "$LOG_DIR" "$EXPERIMENTS_DIR" "$CLAIMS_DIR/active" "$CLAIMS_DIR/done"

# Gemini/opencode must run from the project dir to find instructions files
cd "$PROJECT_DIR"

# Sanity checks
for f in GEMINI.md RESEARCH_STATE.md; do
    if [ ! -f "$f" ]; then
        echo "ERROR: $f not found in $PROJECT_DIR" >&2
        exit 1
    fi
done

# Count existing iterations to continue numbering.
# Global count across all agents so iteration numbers remain unique per script
# (even if two agents run simultaneously, each picks its own N from a snapshot).
LAST_ITER=$(ls "$LOG_DIR"/iter_*.log 2>/dev/null | \
            sed -E 's|.*/iter_([0-9]+)(_[a-zA-Z0-9]+)?\.log|\1|' | \
            sort -n | tail -1 || true)
START_ITER=$(( ${LAST_ITER:-0} + 1 ))
END_ITER=$(( START_ITER + MAX_ITER - 1 ))

# Backend configuration
get_model() {
    case "$1" in
        gemini)     echo "gemini-3-flash-preview" ;;
        gemini-pro) echo "gemini-3.1-pro-preview" ;;
        opencode)   echo "nim/moonshotai/kimi-k2.5" ;;
        *)          echo "$1" ;;  # allow arbitrary model strings
    esac
}

next_backend() {
    case "$1" in
        gemini)     echo "gemini-pro" ;;
        gemini-pro) echo "opencode" ;;
        opencode)   echo "" ;;  # no more fallbacks
        *)          echo "" ;;
    esac
}

run_agent() {
    local backend="$1"
    local prompt="$2"
    local logfile="$3"
    local model
    model=$(get_model "$backend")

    case "$backend" in
        gemini|gemini-pro)
            gemini -p "$prompt" \
                   -m "$model" \
                   --yolo \
                   -o text \
                   2>&1 | tee "$logfile"
            ;;
        opencode)
            if [ -z "${NVIDIA_NIM_API_KEY:-}" ]; then
                echo "ERROR: NVIDIA_NIM_API_KEY not set. Export it before running opencode backend." >&2
                return 1
            fi
            PYTHONUNBUFFERED=1 \
            opencode run "$prompt" \
                   -m "$model" \
                   --dir "$PROJECT_DIR" \
                   2>&1 | tee "$logfile"
            ;;
    esac
}

# Check if output indicates quota exhaustion
is_quota_error() {
    local logfile="$1"
    grep -qiE 'quota|RESOURCE_EXHAUSTED|rate.limit|429|Too Many Requests|exhausted' "$logfile" 2>/dev/null
}

CURRENT_BACKEND="$BACKEND"

echo "=== Chameleon Research Loop (Phase 2: Novel Methods) ==="
echo "  Agent ID:   $AGENT_ID"
echo "  Backend:    $CURRENT_BACKEND (model: $(get_model "$CURRENT_BACKEND"))"
echo "  Iterations: $START_ITER to $END_ITER"
echo "  Cooldown:   ${COOLDOWN}s between iterations"
echo "  Project:    $PROJECT_DIR"
echo "  Logs:       $LOG_DIR"
echo ""

for i in $(seq "$START_ITER" "$END_ITER"); do
    ITER_LOG="$LOG_DIR/iter_${i}${TAG}.log"
    echo "--- Iteration $i / $END_ITER  [agent=$AGENT_ID]  [$(get_model "$CURRENT_BACKEND")]  ($(date)) ---"

    # Parallel-swarm-only instructions. Empty string in solo mode so the prompt
    # is identical to the well-tested single-loop version.
    if [ "$AGENT_ID" = "solo" ]; then
        SWARM_BLOCK=""
    else
        SWARM_BLOCK="PARALLEL SWARM CONTEXT:
You are agent '$AGENT_ID' in a 3-agent swarm. Two other agents may be working
concurrently. Coordinate as follows:

1. BEFORE picking a hypothesis, list claims/active/ and read any files there.
   Those are hypotheses other agents are currently testing — pick something
   DIFFERENT so we don't duplicate compute. Also glance at claims/done/ for
   recently completed work.
2. ONCE you pick a hypothesis, write a short claim file:
       echo '<one-line hypothesis summary>' > claims/active/iter_${i}_${AGENT_ID}.md
   Do this BEFORE running the experiment.
3. When you run Python, WRAP it with the compute lock so only one experiment
   runs on the GPU/CPU at a time across the swarm:
       bash scripts/run_locked.sh $PYTHON_EXE experiments/<your_script>.py
   (Not doing this will cause GPU contention and slow all agents.)
4. AFTER the experiment finishes and you have appended findings, move your
   claim to done:
       mv claims/active/iter_${i}_${AGENT_ID}.md claims/done/
5. All file names below must include the _${AGENT_ID} suffix to avoid
   clobbering sibling agents' files.
"
    fi

    PROMPT="You are iteration $i of an autonomous research loop (Phase 2: Novel Methods).

${SWARM_BLOCK}Read GEMINI.md for full instructions, then read RESEARCH_STATE.md carefully.

This is Phase 2. Phase 1 (implicit-solvent descriptor engineering, iterations 1-16)
is COMPLETE and proved to be a dead end that doesn't generalize. You must explore
FUNDAMENTALLY NEW approaches — new simulation methods, techniques from other fields,
novel physical insights, anything creative that hasn't been tried.

Pick an untested idea from 'Next directions' in RESEARCH_STATE.md, OR propose a
completely new idea that nobody has listed. The best contributions are surprises.

Write a self-contained Python experiment in experiments/iter_${i}${TAG}_<short_name>.py.

CRITICAL RULES:
1. State your hypothesis and physical reasoning in RESEARCH_STATE.md BEFORE coding.
2. Do NOT compute AUC or any classification metric on the 49-molecule benchmark.
   Report raw descriptor values for protac_1/2/3 only (qualitative check).
   Write descriptors for all 49 molecules to experiments/iter_${i}${TAG}_descriptors.json.
3. Do NOT import sklearn.metrics or use evaluate.py.
4. Do NOT tweak CI weights, add vetoes, or do any descriptor engineering.

STDOUT CAPTURE (CRITICAL — without this you will see NO output):
1. Put this at the VERY TOP of your Python script (before ANY other imports):
   import sys, os
   os.environ['PYTHONUNBUFFERED'] = '1'
   sys.stdout = sys.stderr  # Forces output through stderr (unbuffered in Bun)
2. Write ALL results to experiments/iter_${i}${TAG}_output.txt as reliable backup:
   with open('experiments/iter_${i}${TAG}_output.txt', 'w') as f:
       f.write(output_text)
   print(output_text)  # Now goes through stderr, which Bun doesn't buffer

To run Python: $PYTHON_EXE experiments/<your_script>.py
  (In swarm mode — prepend 'bash scripts/run_locked.sh' as shown above.)
Do NOT use cmd.exe //c. Use forward slashes in paths.
Pipeline code is in chameleon_local/ (LOCAL COPY within this project).
User PROTACs: chameleon_local/user_protacs.tsv
Keep experiments under 200 lines and under 15 minutes runtime.
Do NOT revisit dead ends listed in RESEARCH_STATE.md or GEMINI.md.
Ignore 'AttachConsole failed' errors — harmless noise.

After running, append findings to RESEARCH_STATE.md under '## Iteration $i (agent $AGENT_ID)'.
Update 'Next directions' (remove tested, add new ideas). Commit with:
'iter $i [$AGENT_ID]: <what you tested> -- <result>'"

    # Run agent, capturing exit code even on failure
    set +e
    run_agent "$CURRENT_BACKEND" "$PROMPT" "$ITER_LOG"
    EXIT_CODE=${PIPESTATUS[0]}
    set -e

    echo ""
    echo "  -> Iteration $i [$AGENT_ID] finished (exit=$EXIT_CODE) at $(date)"

    # --- Post-iteration: check for output file ---
    OUTPUT_FILE="$EXPERIMENTS_DIR/iter_${i}${TAG}_output.txt"
    if [ -f "$OUTPUT_FILE" ]; then
        echo "  -> Output file found: $OUTPUT_FILE ($(wc -l < "$OUTPUT_FILE") lines)"
    else
        echo "  -> WARNING: No output file. Experiment may have failed silently."
    fi

    # --- Post-iteration: blind evaluation (agent never sees this) ---
    DESCRIPTOR_FILE="$EXPERIMENTS_DIR/iter_${i}${TAG}_descriptors.json"
    if [ -f "$DESCRIPTOR_FILE" ]; then
        echo "  -> Running blind evaluation on descriptors..."
        "$PYTHON_EXE" scripts/blind_evaluate.py \
            "$DESCRIPTOR_FILE" \
            chameleon_local/labelled_set.tsv \
            holdout_set.tsv \
            > "$LOG_DIR/iter_${i}${TAG}_eval.log" 2>&1 || true
        echo "  -> Evaluation saved to $LOG_DIR/iter_${i}${TAG}_eval.log"
    fi

    # --- Post-iteration: sweep any stale claim files for this iteration ---
    # (agent should have moved them itself, but be defensive if it crashed)
    if [ "$AGENT_ID" != "solo" ]; then
        STALE_CLAIM="$CLAIMS_DIR/active/iter_${i}_${AGENT_ID}.md"
        if [ -f "$STALE_CLAIM" ]; then
            mv "$STALE_CLAIM" "$CLAIMS_DIR/done/" 2>/dev/null || true
        fi
    fi

    # Check for quota exhaustion and auto-fallback
    if [ "$EXIT_CODE" -ne 0 ] && is_quota_error "$ITER_LOG"; then
        NEXT=$(next_backend "$CURRENT_BACKEND")
        if [ -n "$NEXT" ]; then
            echo "  -> QUOTA EXHAUSTED on $CURRENT_BACKEND. Switching to $NEXT ($(get_model "$NEXT"))."
            CURRENT_BACKEND="$NEXT"
            # Retry this iteration with the new backend
            echo "  -> Retrying iteration $i with $CURRENT_BACKEND..."
            set +e
            run_agent "$CURRENT_BACKEND" "$PROMPT" "$ITER_LOG"
            EXIT_CODE=${PIPESTATUS[0]}
            set -e
            echo "  -> Retry finished (exit=$EXIT_CODE) at $(date)"
        else
            echo "  -> ALL BACKENDS EXHAUSTED. Stopping loop."
            echo "  -> Completed iterations: $START_ITER to $((i - 1))"
            exit 1
        fi
    elif [ "$EXIT_CODE" -ne 0 ]; then
        echo "  -> ERROR: Agent exited with code $EXIT_CODE. Pausing 120s."
        sleep 120
    fi

    # Cooldown between iterations
    if [ "$i" -lt "$END_ITER" ]; then
        # Longer cooldown for opencode (40 RPM limit)
        if [ "$CURRENT_BACKEND" = "opencode" ]; then
            ACTUAL_COOLDOWN=$(( COOLDOWN > 90 ? COOLDOWN : 90 ))
        else
            ACTUAL_COOLDOWN="$COOLDOWN"
        fi
        echo "  -> Cooling down ${ACTUAL_COOLDOWN}s before next iteration..."
        sleep "$ACTUAL_COOLDOWN"
    fi
done

echo ""
echo "=== Research loop complete: iterations $START_ITER-$END_ITER ==="
echo "  Backend used: $CURRENT_BACKEND"
echo "  Review findings in: $PROJECT_DIR/RESEARCH_STATE.md"
echo "  Iteration logs in:  $LOG_DIR/"
echo "  Evaluation logs in: $LOG_DIR/iter_*_eval.log"
echo "  Run 'git log --oneline' to see commits."
