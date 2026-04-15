# Chameleon Autoresearch

Can Agents perform the scientific method? (Spoiler: not very well/kinda?)

TLDR: I tried to get Agents to generate hypotheses and test them, rinse and repeat, Ralph-loop style (Except unlike a ralph loop, it's not really iterative improvement towards a defined goal, it's more open-ended) 
Inspired by Terence Tao: "AI has driven the cost of idea generation to almost zero", and karpathy's autoresearch

Preliminary results are... mixed. It can generate ideas and test them, but they're not amazing? Anyway, it hasn't really found anything good. 

Feels like a chemistry undergrad kind of level of research


--yes, claude wrote everything below. I won't apologize for the slop, however, the ideas are my own, and I think autonomous computational chemistry is a worthwhile thing to play with. 




A small experiment in **autonomous scientific research**. Three coding agents
generate hypotheses, write Python experiments, run them, interpret results, and
commit findings — concurrently and without human input — while a fourth agent
monitors them and a human meta-optimizer steers the direction.

The scientific target is narrow: improve a published pipeline for predicting
**chameleonicity** (conformational shape-shifting) of PROTACs and macrocycles.
The autoresearch mechanism is the actually-interesting part; the chemistry is
the substrate it runs on.

This repository is a **snapshot**, not a finished product. It is published for
transparency about what the setup looks like at this stage of the work. No
scientific breakthroughs are claimed.

## Architecture at a glance

```
┌─────────────────────────────────────────────┐
│  Human + Opus meta-optimizer (interactive)  │
│  edits GEMINI.md, RESEARCH_STATE.md, scripts│
└────────────┬────────────────────────────────┘
             │
             ↓
┌─────────────────────────────────────────────┐
│  Sonnet 4.6 monitor (headless, 30 min cycle)│
│  writes monitor/logs/check_*.md             │
│  sanity-check only, no authority to edit    │
└────────────┬────────────────────────────────┘
             │
             ↓
┌─────────────────────────────────────────────┐
│  3 × research agents (opencode + Kimi K2.5  │
│  or Gemini 3 Flash)                         │
│                                             │
│   Agent 1     Agent 2     Agent 3           │
│      │          │            │              │
│      └──────────┴────────────┘              │
│        claims/ dir (coordinate ideas)       │
│                  │                          │
│                  ↓                          │
│       scripts/run_locked.sh (compute lock)  │
│                  │                          │
│                  ↓                          │
│              GPU / CPU                      │
│                                             │
│   each appends → RESEARCH_STATE.md          │
│   each commits → git                        │
└─────────────────────────────────────────────┘
```

**Coordination is entirely file-based.** No sockets, no daemons, no shared
memory. The filesystem is the communication primitive.

## What each layer does

### Meta-optimizer (human + Opus)
Reads commits, prunes dead ends, rewrites `GEMINI.md` (the agent's instructions),
restructures `RESEARCH_STATE.md` (findings + priorities), decides when the loop
architecture itself needs to change. Not on a schedule — session-driven.

### Monitor (Sonnet 4.6, headless every 30 min)
A separate Claude Code process invoked via `claude -p`. Runs sanity checks,
writes one-page markdown status files. Read-only: no authority to launch,
stop, or edit. Its job is to notice when the loop is genuinely broken vs.
just thinking or running a long compute job.

### Research agents (opencode or gemini, 3 in parallel)
Each agent runs the full scientific cycle:

1. Read `GEMINI.md` and `RESEARCH_STATE.md`
2. Check `claims/active/` to avoid duplicating sibling work
3. Pick a hypothesis, write a claim file, state physical reasoning first
4. Write a Python experiment in `experiments/iter_N_aX_*.py`
5. Run it via `bash scripts/run_locked.sh python ...` — a global mkdir-based
   lock serializes GPU/CPU use across the swarm, so only one experiment runs
   at a time
6. Write findings to `RESEARCH_STATE.md`
7. `git commit` with the iteration number and agent tag
8. Move the claim to `claims/done/`, cooldown, next iteration

Backend auto-fallback: if an agent hits a quota error, the loop swaps
`gemini → gemini-pro → opencode` automatically.

## Running it

Requires:
- A working `chameleon` conda env with RDKit + OpenMM (see `chameleon_local/`)
- One of: `gemini` CLI, or `opencode` + an `NVIDIA_NIM_API_KEY` env var
- `claude` CLI for the monitor
- Git Bash on Windows (or Linux/macOS — untested on macOS)

```bash
export NVIDIA_NIM_API_KEY=<your key>
bash start_session.sh 3 5 opencode    # 3 agents × 5 iters, Kimi backend
```

Watch progress:
```bash
tail -f logs/agent_a1.log
git log --oneline -10
ls claims/active/
cat monitor/logs/check_$(ls -t monitor/logs/check_*.md | head -1)
```

Stop:
```bash
bash scripts/stop_agents.sh    # kills agents, releases lock
touch monitor/STOP             # stops the monitor after its current check
```

## The scientific substrate

The loop is trying to improve predictions against a curated 49-molecule
labelled set (`chameleon_local/labelled_set.tsv`) of macrocycles, cyclic
peptides, and PROTACs, each carrying a binary chameleonicity label
compiled from multiple literature sources. The labels are informed by —
but not identical to — **ChamelogK**, a chromatography-based *experimental*
technique developed by García-Jiménez et al. (*J. Med. Chem.* **66**, 2023,
DOI [10.1021/acs.jmedchem.3c00823](https://doi.org/10.1021/acs.jmedchem.3c00823))
that yields a continuous measure of a molecule's shape-shifting behavior;
values ≥ 0.6 are conventionally called chameleonic. Our labelled set
borrows some ChamelogK assignments but also uses a broader oral-bioavailability
framing, so there are documented discrepancies (e.g. everolimus is labeled
"chameleon" here despite ChamelogK 0.45); `LITERATURE.md` in the repo lists
them.

To put it plainly: **ChamelogK is a wet-lab method, not a benchmark.** The
"benchmark" is a curated labelled set that uses ChamelogK (among other
sources) to decide labels.

A separate 29-molecule holdout set (`holdout_set.tsv`), built directly from
ChamelogK compounds *not* present in the labelled set, is used by the
overseer — never by the agents — to detect overfitting. This is the
mechanism that caught Phase 1's collapse from AUC 0.910 on the labelled set
to 0.644 on the holdout, and triggered the Phase 1 → Phase 2 transition and
the "benchmark blindness protocol" (agent never computes AUC, overseer
evaluates separately).

Phase 2 (iterations 17–28) has produced a few physically-interesting partial
successes (H-bond network cooperativity, Flory-Huggins 2D interaction density,
dihedral preference statistics from torsional propensities) and many honest
failures. No explicit-solvent MD has been attempted yet, despite being the
gold standard — agents consistently prefer cheap ideas they can ship in 15
minutes over expensive ideas they believe will be correct. Working on fixing
that bias.

See `RESEARCH_STATE.md` for the full, unedited research log. It includes the
agents' hypotheses, code, results, and self-assessments, committed atomically
iteration by iteration.

## What this is (and is not)

This **is**:
- An honest working system for autonomous idea → code → experiment → finding
  cycles on a narrow technical problem, with cheap file-based coordination
  between parallel agents and a separate monitor layer
- A demonstration that the *mechanics* of autonomous research can be stood up
  with a few hundred lines of bash + markdown instructions, using commodity
  LLM backends
- An empirical record of what the resulting output actually looks like, warts
  included

This **is not**:
- A finished product or a research breakthrough
- A system that produces PhD-level science without meta-optimization
- Free of implicit overfitting, failed experiments, duplicated work, or
  cosmetic bugs — all of those are present and visible in the logs

## Known limitations (as observed so far)

- Agents prefer cheap 2D topological descriptors over expensive MD, even when
  rules allow the latter
- No synthesis instinct — each iteration is isolated, partial wins don't
  compound automatically (meta-optimizer has to prompt that)
- When experiments fail due to implementation bugs (not science), agents
  tend to record "FAIL" and abandon the idea instead of fixing the code
- Iteration numbering occasionally collides when agents read stale state —
  mostly cosmetic, doesn't break coordination
- The meta-optimizer (a human + Opus session) does substantial lifting;
  removing it degrades the system considerably

## License

MIT. This is a research sketch. No warranty.
