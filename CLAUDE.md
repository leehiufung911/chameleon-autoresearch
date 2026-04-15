# Chameleon Research Loop

## What This Project Is

An autonomous research loop that uses **Gemini 3 Flash** (via Gemini CLI) as the
inner exploration agent and **Claude** as the meta-optimizer / human overseer.

The goal: explore approaches to improve chameleonicity prediction for PROTACs,
specifically fixing the protac_3 misclassification and improving the pipeline's
ability to distinguish chameleonic vs non-chameleonic bifunctional molecules.

**This is NOT a Ralph loop** (convergent, test-driven). It is an open-ended
scientific exploration loop where each iteration generates a hypothesis, tests
it with code, and records findings. There is no single pass/fail oracle.

## Architecture

```
Claude (this session)          Gemini CLI + gemini-3-flash-preview
 - reviews findings            - reads RESEARCH_STATE.md
 - improves GEMINI.md          - picks next hypothesis
 - improves research_loop.sh   - writes & runs experimental code
 - adjusts research direction  - appends results to RESEARCH_STATE.md
 - meta-optimizes the loop     - git commits each iteration
```

## Claude's Role (YOU)

You are the **meta-optimizer**. You have full authority to change ANY part of
this experiment setup. Nothing is sacred — the loop script, the Gemini
instructions, the research state, the hypothesis list, even the evaluation
criteria can all be modified if you think it will improve the loop.

Your core activities:

1. **Review** — Read `RESEARCH_STATE.md` and `logs/` to see what Gemini found
2. **Diagnose** — Identify if the loop is stuck, going in circles, or exploring
   unproductive directions
3. **Improve** — Edit ANY file in this project: `GEMINI.md` (agent instructions),
   `research_loop.sh` (loop mechanics), `RESEARCH_STATE.md` (seed directions,
   prune dead ends, reframe the problem), or create entirely new files
4. **Restructure** — Change the loop architecture itself if needed: add
   evaluation scripts, change the prompt template, add pre/post-iteration hooks,
   split into sub-loops, whatever serves the research
5. **Evaluate** — Run the benchmark yourself if Gemini claims an improvement:
   `cmd.exe //c "call C:\Users\mic23\miniconda3\Scripts\activate.bat chameleon && python C:\Users\mic23\prototype-embed\chameleon\evaluate.py ..."`

## Key Files

| File | Purpose | Who edits |
|------|---------|-----------|
| `GEMINI.md` | Agent instructions for Gemini CLI | Claude (meta) |
| `RESEARCH_STATE.md` | Accumulated findings, evolving context | Gemini (append), Claude (prune/steer) |
| `research_loop.sh` | The loop script | Claude (meta) |
| `logs/iter_N.log` | Raw output of each iteration | Loop script |
| `experiments/` | Experimental scripts Gemini writes | Gemini |

## The Chameleon Pipeline (Reference)

Lives at `C:\Users\mic23\prototype-embed\chameleon\`. Key files:

- `chameleon.py` — main pipeline: ETKDGv3 -> MMFF94s -> Butina -> PSA/Rg/IMHB -> CI
- `gbsa_rescore.py` — SMIRNOFF + OpenMM GBSA OBC1 rescoring
- `evaluate.py` — benchmark scorer (AUC, F1)
- `benchmark.tsv` / `labelled_set.tsv` — 49-molecule benchmark (AUC 0.910 baseline)
- `user_protacs.tsv` — 3 test PROTACs (protac_1,2 = chameleonic, protac_3 = non)
- `crest_rescore.py` — CREST wrapper (WSL2, xtb+CREST installed)

**Conda env**: `chameleon` — activate with:
```
cmd.exe //c "call C:\Users\mic23\miniconda3\Scripts\activate.bat chameleon && <command>"
```

**GPU**: RTX 3070, use OpenCL (not CUDA) for OpenMM.

## Current Benchmark State

- MMFF CI: AUC **0.910** on 49-mol labelled set (F1 0.933 at threshold 2.69)
- GBSA CI: AUC **0.867** (limited by q_size=12 from 64-conf cap)
- protac_1: CI 4.13 -> chameleonic (correct)
- protac_2: CI 4.77 -> chameleonic (correct)
- protac_3: CI 3.80 -> marginal (WRONG: should be non-chameleonic)

## What We Already Tried and Ruled Out

These are dead ends — do NOT let Gemini re-explore them:

1. **Population-weighted CI (Q20-by-MMFF-energy)**: Failed. Same-scaffold PROTACs
   have near-identical MMFF energies; reweighting randomly shuffles rank order.

2. **Chloroform dielectric in GBSA (eps 2.02 -> 4.7)**: No-op. dE-based Q20 split
   is PSA-dominated; changing dielectric rescales dE without reshuffling Q20.

3. **CREST iMTD-GC + GFN-FF + ALPB(chloroform)**: Finds only compact basin for
   protac_3 (267 confs, all Rg 4.6-5.5 A). GFN-FF prefers compact by +47.5 kcal/mol
   over extended. GFN2-xTB still +13.7 kcal/mol. This is a force-field-level bias,
   not a sampling failure.

4. **Root cause**: Every implicit-solvent, enthalpy-only, single-molecule method
   (MMFF94s, GFN-FF, GFN2-xTB) strongly prefers compact IMHB-rich states. The
   paper's MD in explicit CHCl3 finds extended Rg 7-9.5 A because of effects
   implicit solvent cannot capture (explicit solute-solvent contacts, solvent
   entropy, conformational entropy at 300K).

## Productive Directions to Explore

These are hypotheses that have NOT been tested:

- Linker-specific descriptors (rotatable bond count, topological polar surface
  area of the linker alone, gauche preference of linker type)
- Entropy-aware scoring (conformational entropy from ensemble diversity as a
  correction to enthalpy-only CI)
- Descriptor-space ML (train a lightweight classifier on the 49-mol descriptors
  instead of hand-tuned CI weights)
- Pre-MMFF descriptors (compute PSA/Rg on raw ETKDG geometries before MMFF
  collapses extended conformers — protac_3 has Rg up to 9.97 pre-MMFF)
- Hybrid CI (combine MMFF CI channels with GBSA CI channels, picking whichever
  has better per-channel AUC)
- Explicit-solvent MD via OpenMM (expensive but would match the paper's method)
- Better IMHB scoring (topology-aware, accounting for linker type)
