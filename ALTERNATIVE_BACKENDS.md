# Alternative Backends: opencode + NVIDIA NIM

If you hit Gemini quota limits, switch to **opencode** with free NVIDIA NIM models.

## Setup

opencode v1.2.27 is already installed (`C:\Users\mic23\AppData\Roaming\npm\opencode`).
NVIDIA NIM API key is in `C:\Users\mic23\free-claude-code\.env`.

### Verify it works

```bash
# Kimi K2.5
opencode run "Say hello" -m nim/moonshotai/kimi-k2.5

# GLM 5 (being deprecated April 20, 2026 — prefer GLM 5.1 when available)
opencode run "Say hello" -m nim/z-ai/glm5
```

### Available NIM models (confirmed working 2026-04-13)

| Model ID | Context | Notes |
|----------|---------|-------|
| `nim/moonshotai/kimi-k2.5` | 262K | Best reasoning (has thinking mode). Recommended. |
| `nim/z-ai/glm5` | 202K | Fast, direct. Deprecated April 20 -> use glm-5.1. |
| `nim/stepfun-ai/step-3.5-flash` | — | Current default in free-claude-code |
| `nim/qwen/qwen3.5-397b-a17b` | — | Alternative |
| `nim/nvidia/nemotron-3-super-120b-a12b` | — | Alternative |

Rate limit: **40 requests per 60-second window** (free tier).

## Switching the loop script

Replace the Gemini invocation in `research_loop.sh`. Change:

```bash
MODEL="gemini-3-flash-preview"
```

to:

```bash
MODEL="nim/moonshotai/kimi-k2.5"
BACKEND="opencode"  # instead of gemini
```

And replace the Gemini call block with:

```bash
    opencode run "$PROMPT" \
           -m "$MODEL" \
           --dir "$PROJECT_DIR" \
           2>&1 | tee "$LOG_DIR/iter_${i}.log"
    EXIT_CODE=${PIPESTATUS[0]}
```

### Full drop-in replacement script

Save as `research_loop_nim.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail

MAX_ITER="${1:-5}"
COOLDOWN="${2:-60}"  # longer cooldown for 40 RPM limit
MODEL="${3:-nim/moonshotai/kimi-k2.5}"
PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
LOG_DIR="$PROJECT_DIR/logs"
EXPERIMENTS_DIR="$PROJECT_DIR/experiments"

mkdir -p "$LOG_DIR" "$EXPERIMENTS_DIR"
cd "$PROJECT_DIR"

for f in GEMINI.md RESEARCH_STATE.md; do
    if [ ! -f "$f" ]; then
        echo "ERROR: $f not found in $PROJECT_DIR" >&2
        exit 1
    fi
done

LAST_ITER=$(ls "$LOG_DIR"/iter_*.log 2>/dev/null | \
            sed 's/.*iter_\([0-9]*\)\.log/\1/' | sort -n | tail -1 || true)
START_ITER=$(( ${LAST_ITER:-0} + 1 ))
END_ITER=$(( START_ITER + MAX_ITER - 1 ))

echo "=== Chameleon Research Loop (opencode + NIM) ==="
echo "  Model:      $MODEL"
echo "  Iterations: $START_ITER to $END_ITER"
echo "  Cooldown:   ${COOLDOWN}s between iterations"
echo ""

for i in $(seq "$START_ITER" "$END_ITER"); do
    echo "--- Iteration $i / $END_ITER  ($(date)) ---"

    PROMPT="You are iteration $i of an autonomous research loop.

Read GEMINI.md for full instructions, then read RESEARCH_STATE.md carefully.

Pick the SINGLE most promising untested hypothesis from the 'Next directions'
section. Write a self-contained Python experiment in experiments/iter_${i}_<short_name>.py
to test it. Run it. Append your findings to RESEARCH_STATE.md under
'## Iteration $i'. Update the 'Next directions' section (remove what you tested,
add new ideas). Commit everything with: 'iter $i: <what you tested> -- <result>'.

IMPORTANT:
- To run Python, use:
  cmd.exe //c \"call C:\\Users\\mic23\\miniconda3\\Scripts\\activate.bat chameleon && python experiments/<your_script>.py\"
- The pipeline code is at ../chameleon/ (relative to this directory).
- The benchmark is at ../chameleon/labelled_set.tsv and ../chameleon/benchmark.tsv
- User PROTACs are at ../chameleon/user_protacs.tsv
- Keep your experiment under 200 lines. Report numbers, not feelings.
- Do NOT revisit dead ends listed in RESEARCH_STATE.md."

    set +e
    opencode run "$PROMPT" \
           -m "$MODEL" \
           --dir "$PROJECT_DIR" \
           2>&1 | tee "$LOG_DIR/iter_${i}.log"
    EXIT_CODE=${PIPESTATUS[0]}
    set -e

    echo ""
    echo "  -> Iteration $i finished (exit=$EXIT_CODE) at $(date)"

    if [ "$EXIT_CODE" -ne 0 ]; then
        echo "  -> ERROR: opencode exited with code $EXIT_CODE. Pausing 120s."
        sleep 120
    fi

    if [ "$i" -lt "$END_ITER" ]; then
        echo "  -> Cooling down ${COOLDOWN}s before next iteration..."
        sleep "$COOLDOWN"
    fi
done

echo ""
echo "=== Research loop complete: iterations $START_ITER-$END_ITER ==="
```

## Key differences from Gemini

| | Gemini CLI | opencode + NIM |
|---|---|---|
| Cost | Free (Google AI Pro quota) | Free (NIM free tier) |
| Rate limit | Higher (Pro sub) | 40 RPM / 5000 credits |
| YOLO mode | `--yolo` | All tools auto-allowed by default |
| Project instructions | Reads `GEMINI.md` | Reads `.opencode/instructions.md` |
| Sandbox | `--sandbox` (Docker) | No built-in sandbox |
| Worktree | `-w` flag | No built-in worktree |

### opencode project instructions

opencode reads `.opencode/instructions.md` (not `GEMINI.md`). Create a symlink
or copy:

```bash
mkdir -p .opencode
cp GEMINI.md .opencode/instructions.md
```

Or set it in `.opencode/config.json`:
```json
{
  "instructions": "GEMINI.md"
}
```

## Mixing backends

You can run some iterations on Gemini and some on NIM. The research state
persists in `RESEARCH_STATE.md` and git — both backends read and write the
same files. Just switch scripts:

```bash
# First 5 on Gemini
bash research_loop.sh 5

# Hit quota? Switch to NIM for next 5
bash research_loop_nim.sh 5
```
