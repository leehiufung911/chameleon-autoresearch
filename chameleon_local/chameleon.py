"""
chameleon.py - tractable a-priori chameleonicity scoring for bifunctional
molecules (PROTACs, bRo5, macrocycles).

Pipeline
--------
1. ETKDGv3 conformer ensemble (typically 300-500 confs)
2. MMFF94s minimization (multi-threaded)
3. RMSD prune (Butina clustering) -> unique conformers
4. For each unique conformer compute:
     - 3D PSA   : SASA restricted to N, O and H bonded to N/O (Bondi radii)
     - Rg       : mass-weighted radius of gyration
     - IMHB     : geometric intramolecular H-bond count
                  (donor-acceptor > 3 bonds, H..A <= 2.5 A, D-H..A >= 120 deg)
     - MMFF E   : minimized MMFF94s energy
5. Build two conformer sub-ensembles without running MD:
     - "apolar" proxy: Boltzmann weight by E_apolar  = E_MMFF - 2.0 * n_IMHB
                       (IMHB is worth ~2 kcal/mol in chloroform-like env)
     - "polar"  proxy: Boltzmann weight by E_polar   = E_MMFF + 0.05 * PSA3D
                       (penalty for buried polarity, tuned to match
                        desolvation-of-polar-surface cost ~0.05 kcal/mol/A^2)
   Then compute weighted <PSA3D>, <Rg>, <IMHB> under each sub-ensemble.
6. Report descriptors and a composite "chameleonic index" CI:
     dPSA      = <PSA>_polar  - <PSA>_apolar       (A^2)
     dRg       = <Rg>_polar   - <Rg>_apolar        (A)
     dIMHB     = <IMHB>_apolar - <IMHB>_polar
     CI        = (dPSA / sqrt(MW)) + 5*dRg + 4*dIMHB

Also reports the raw (unweighted) min / max / std of PSA3D and Rg, which is
the simpler signal used by Rossi Sebastiano 2018 and David 2021.

Usage
-----
    py chameleon.py "<SMILES>" --name MZ1 --nconf 300

    py chameleon.py --batch smiles.txt       # one "NAME<TAB>SMILES" per line
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from dataclasses import dataclass, asdict
from typing import List, Optional, Tuple

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolDescriptors, rdFreeSASA, Descriptors
from rdkit.ML.Cluster import Butina

RDLogger.DisableLog("rdApp.*")

KT = 0.5924  # kcal/mol at 298 K
HBOND_EN_APOLAR = 2.0   # reward per IMHB in apolar sub-ensemble (kcal/mol)
HBOND_EN_POLAR = 1.5    # penalty per IMHB in polar sub-ensemble (water H-bond
                        # would be preferred over intramolecular)
POLAR_SOLV_KCAL_PER_A2 = 0.05   # reward per A^2 of exposed polar SASA in water


# --------------------------------------------------------------------------- #
# conformer generation
# --------------------------------------------------------------------------- #
def embed_conformers(
    mol: Chem.Mol, n_conf: int, n_threads: int, seed: int = 0xC0FFEE
) -> Chem.Mol:
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.numThreads = n_threads
    params.pruneRmsThresh = 0.5
    params.useSmallRingTorsions = True
    params.useRandomCoords = True
    params.enforceChirality = True
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=params)
    if len(cids) == 0:
        # fallback: drop chirality enforcement and random coords flag
        params.enforceChirality = False
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=params)
    return mol, list(cids)


def minimize(mol: Chem.Mol, n_threads: int) -> List[float]:
    results = AllChem.MMFFOptimizeMoleculeConfs(
        mol, maxIters=500, numThreads=n_threads, mmffVariant="MMFF94s"
    )
    # each entry: (not_converged_flag, energy)
    return [float(r[1]) for r in results]


def butina_prune(
    mol: Chem.Mol, energies: List[float], rms_cut: float = 0.75
) -> Tuple[List[int], List[float]]:
    n = mol.GetNumConformers()
    if n <= 1:
        return list(range(n)), energies
    dmat = AllChem.GetConformerRMSMatrix(mol, prealigned=False)
    clusters = Butina.ClusterData(
        dmat, n, rms_cut, isDistData=True, reordering=True
    )
    # Keep the *lowest-energy* member of each cluster (Butina's centroid is
    # just the first picked; energy-based is more physical).
    keep = []
    for cl in clusters:
        best = min(cl, key=lambda i: energies[i])
        keep.append(best)
    keep.sort()
    return keep, [energies[i] for i in keep]


# --------------------------------------------------------------------------- #
# per-conformer descriptors
# --------------------------------------------------------------------------- #
def polar_atom_mask(mol: Chem.Mol) -> np.ndarray:
    mask = np.zeros(mol.GetNumAtoms(), dtype=bool)
    for a in mol.GetAtoms():
        sym = a.GetSymbol()
        if sym in ("N", "O"):
            mask[a.GetIdx()] = True
        elif sym == "H":
            # Ertl TPSA includes NH / OH hydrogens; include them in 3D PSA
            if any(n.GetSymbol() in ("N", "O") for n in a.GetNeighbors()):
                mask[a.GetIdx()] = True
    return mask


def compute_3d_psa(mol: Chem.Mol, conf_id: int, radii) -> float:
    # CalcSASA sets 'SASA' prop per atom
    rdFreeSASA.CalcSASA(mol, radii, confIdx=conf_id)
    pmask = polar_atom_mask(mol)
    psa = 0.0
    for atom in mol.GetAtoms():
        if pmask[atom.GetIdx()]:
            psa += float(atom.GetProp("SASA"))
    return psa


def compute_rg(mol: Chem.Mol, conf_id: int) -> float:
    conf = mol.GetConformer(conf_id)
    coords = np.asarray(conf.GetPositions())
    masses = np.asarray([a.GetMass() for a in mol.GetAtoms()])
    com = np.average(coords, axis=0, weights=masses)
    sq = np.sum((coords - com) ** 2, axis=1)
    return float(np.sqrt(np.sum(masses * sq) / masses.sum()))


def bond_path_distance_matrix(mol: Chem.Mol) -> np.ndarray:
    return Chem.GetDistanceMatrix(mol)


def compute_imhb(
    mol: Chem.Mol,
    conf_id: int,
    topo_dist: np.ndarray,
    d_max: float = 2.5,
    ang_min_deg: float = 120.0,
    min_topo: int = 4,
) -> int:
    """Count unique donor-H -> acceptor intramolecular H-bonds.
    A donor is N-H or O-H; an acceptor is any N or O with at least one
    lone pair. Use distance + angle geometric criteria like CPPTRAJ.
    """
    conf = mol.GetConformer(conf_id)
    coords = np.asarray(conf.GetPositions())

    donor_pairs = []
    acceptors = []
    for a in mol.GetAtoms():
        sym = a.GetSymbol()
        if sym in ("N", "O"):
            acceptors.append(a.GetIdx())
            for n in a.GetNeighbors():
                if n.GetSymbol() == "H":
                    donor_pairs.append((a.GetIdx(), n.GetIdx()))
    if not donor_pairs or not acceptors:
        return 0

    cos_min = math.cos(math.radians(180.0 - ang_min_deg))  # not used, see below
    ang_min_rad = math.radians(ang_min_deg)

    count = 0
    used_donors = set()
    for d_idx, h_idx in donor_pairs:
        if h_idx in used_donors:
            continue
        for a_idx in acceptors:
            if a_idx == d_idx:
                continue
            if topo_dist[d_idx, a_idx] < min_topo:
                continue
            d_ha = np.linalg.norm(coords[h_idx] - coords[a_idx])
            if d_ha > d_max:
                continue
            v1 = coords[d_idx] - coords[h_idx]
            v2 = coords[a_idx] - coords[h_idx]
            nv1 = np.linalg.norm(v1)
            nv2 = np.linalg.norm(v2)
            if nv1 < 1e-6 or nv2 < 1e-6:
                continue
            cosang = float(np.dot(v1, v2) / (nv1 * nv2))
            cosang = max(-1.0, min(1.0, cosang))
            ang = math.acos(cosang)
            if ang >= ang_min_rad:
                count += 1
                used_donors.add(h_idx)
                break
    return count


# --------------------------------------------------------------------------- #
# scoring
# --------------------------------------------------------------------------- #
@dataclass
class Summary:
    name: str
    smiles: str
    mw: float
    heavy_atoms: int
    tpsa_2d: float
    n_conf_embedded: int
    n_conf_unique: int
    time_seconds: float

    psa3d_min: float
    psa3d_max: float
    psa3d_mean: float
    psa3d_std: float

    rg_min: float
    rg_max: float
    rg_mean: float
    rg_std: float

    imhb_mean: float
    imhb_max: int
    imhb_frac_ge2: float

    psa3d_apolar_wtd: float
    psa3d_polar_wtd: float
    rg_apolar_wtd: float
    rg_polar_wtd: float
    imhb_apolar_wtd: float
    imhb_polar_wtd: float

    # quantile descriptors: bottom/top 20% of conformers sorted by PSA
    psa3d_q_apolar: float
    psa3d_q_polar: float
    rg_q_apolar: float
    rg_q_polar: float
    imhb_q_apolar: float
    imhb_q_polar: float

    dPSA: float
    dRg: float
    dIMHB: float
    chameleonic_index: float
    chameleonic_verdict: str

    # GBSA re-scored sub-ensembles.  Populated only when --gbsa was used.
    gbsa_platform: Optional[str] = None
    gbsa_forcefield: Optional[str] = None
    gbsa_seconds: Optional[float] = None
    psa3d_gb_apolar: Optional[float] = None
    psa3d_gb_polar: Optional[float] = None
    rg_gb_apolar: Optional[float] = None
    rg_gb_polar: Optional[float] = None
    imhb_gb_apolar: Optional[float] = None
    imhb_gb_polar: Optional[float] = None
    dPSA_gb: Optional[float] = None
    dRg_gb: Optional[float] = None
    dIMHB_gb: Optional[float] = None
    chameleonic_index_gb: Optional[float] = None
    chameleonic_verdict_gb: Optional[str] = None


def boltzmann_weights(energies_kcal: np.ndarray) -> np.ndarray:
    e = energies_kcal - float(np.min(energies_kcal))
    w = np.exp(-e / KT)
    s = float(w.sum())
    if s <= 0 or not np.isfinite(s):
        return np.ones_like(e) / len(e)
    return w / s


def summarize(
    name: str,
    smiles: str,
    n_conf: int = 300,
    n_threads: int = 0,
    rms_cut: float = 0.75,
    verbose: bool = True,
    sdf_out: Optional[str] = None,
    gbsa: bool = False,
    gbsa_max_conf: int = 64,
    gbsa_platform: str = "CUDA",
    gbsa_ff: str = "openff-2.1.0",
) -> Summary:
    if n_threads <= 0:
        import os
        n_threads = max(1, (os.cpu_count() or 4))

    t0 = time.time()
    mol0 = Chem.MolFromSmiles(smiles)
    if mol0 is None:
        raise ValueError(f"bad SMILES: {smiles}")
    mw = Descriptors.MolWt(mol0)
    heavy = mol0.GetNumHeavyAtoms()

    if verbose:
        print(
            f"[{name}] MW={mw:6.1f}  heavy={heavy:3d}  "
            f"embedding {n_conf} confs on {n_threads} threads ..."
        )

    mol, cids = embed_conformers(mol0, n_conf, n_threads)
    if verbose:
        print(f"[{name}] embedded {len(cids)} confs, minimising MMFF94s ...")
    energies = minimize(mol, n_threads)

    keep, k_energies = butina_prune(mol, energies, rms_cut=rms_cut)
    if verbose:
        print(f"[{name}] pruned to {len(keep)} unique confs "
              f"(rms_cut={rms_cut}), computing descriptors ...")

    radii = rdFreeSASA.classifyAtoms(mol)
    topo = bond_path_distance_matrix(mol)

    psa_list, rg_list, imhb_list = [], [], []
    for cid in keep:
        psa_list.append(compute_3d_psa(mol, cid, radii))
        rg_list.append(compute_rg(mol, cid))
        imhb_list.append(compute_imhb(mol, cid, topo))

    psa = np.asarray(psa_list)
    rg = np.asarray(rg_list)
    imhb = np.asarray(imhb_list, dtype=float)
    e_mmff = np.asarray(k_energies)

    # Two biased sub-ensembles (cheap surrogates for explicit solvent MD):
    #   apolar: reward IMHBs (approximate chloroform stabilisation)
    #   polar : reward exposed polar SASA AND penalise IMHBs (water would
    #           rather H-bond with the solute than let the solute self-pair)
    e_apolar = e_mmff - HBOND_EN_APOLAR * imhb
    e_polar = e_mmff - POLAR_SOLV_KCAL_PER_A2 * psa + HBOND_EN_POLAR * imhb
    w_apolar = boltzmann_weights(e_apolar)
    w_polar = boltzmann_weights(e_polar)

    psa_ap = float(np.sum(w_apolar * psa))
    psa_po = float(np.sum(w_polar * psa))
    rg_ap = float(np.sum(w_apolar * rg))
    rg_po = float(np.sum(w_polar * rg))
    imhb_ap = float(np.sum(w_apolar * imhb))
    imhb_po = float(np.sum(w_polar * imhb))

    # Quantile proxies: bottom 20% of conformers by PSA3D == "apolar-like",
    # top 20% by PSA3D == "polar-like".  More robust than Boltzmann when the
    # gas-phase MMFF energy spread dwarfs the bias.
    n = len(psa)
    q_size = max(1, n // 5)
    order = np.argsort(psa)
    q_ap_idx = order[:q_size]
    q_po_idx = order[-q_size:]
    psa_q_ap = float(psa[q_ap_idx].mean())
    psa_q_po = float(psa[q_po_idx].mean())
    rg_q_ap = float(rg[q_ap_idx].mean())
    rg_q_po = float(rg[q_po_idx].mean())
    imhb_q_ap = float(imhb[q_ap_idx].mean())
    imhb_q_po = float(imhb[q_po_idx].mean())

    dPSA = psa_po - psa_ap
    dRg = rg_po - rg_ap
    dIMHB = imhb_ap - imhb_po

    # --- Primary chameleonic index: model-free, built from ensemble spread. --
    # The David 2021 and Rossi Sebastiano papers both observe that the raw
    # RANGE of 3D PSA across an unbiased conformer ensemble is the dominant
    # signal.  We combine three orthogonal channels:
    #   (a) size-normalised PSA spread   (polarity dynamics)
    #   (b) mean IMHB count              (polarity-masking machinery)
    #   (c) log radius-of-gyration ratio (folding dynamics)
    psa_range = float(psa.max() - psa.min())
    rg_ratio = float(rg.max() / max(rg.min(), 1e-6))
    term_psa = psa_range / math.sqrt(max(mw, 1.0))
    term_imhb = float(imhb.mean())
    term_rg = math.log(max(rg_ratio, 1.0))
    ci = 2.0 * term_psa + 1.0 * term_imhb + 3.0 * term_rg

    if ci >= 4.0:
        verdict = "chameleonic"
    elif ci >= 2.0:
        verdict = "marginally chameleonic"
    else:
        verdict = "non-chameleonic"

    # ---- optional GBSA rescoring -----------------------------------------
    gb_fields = {
        "gbsa_platform": None,
        "gbsa_forcefield": None,
        "gbsa_seconds": None,
        "psa3d_gb_apolar": None, "psa3d_gb_polar": None,
        "rg_gb_apolar": None, "rg_gb_polar": None,
        "imhb_gb_apolar": None, "imhb_gb_polar": None,
        "dPSA_gb": None, "dRg_gb": None, "dIMHB_gb": None,
        "chameleonic_index_gb": None, "chameleonic_verdict_gb": None,
    }
    if gbsa:
        try:
            import gbsa_rescore as gbr
        except Exception as e:
            print(f"[{name}] GBSA import failed: {e}")
            gbr = None
        if gbr is not None:
            # Pick a subset of at most gbsa_max_conf unique conformers (by
            # lowest MMFF energy) to keep GPU wall-time bounded.
            order_e = np.argsort(e_mmff)
            sel_local = order_e[: max(1, min(gbsa_max_conf, len(keep)))]
            sel_cids = [int(keep[i]) for i in sel_local]
            if verbose:
                print(f"[{name}] GBSA rescore of {len(sel_cids)} confs "
                      f"on {gbsa_platform} ({gbsa_ff})")
            t_gb = time.time()
            try:
                water, apolar = gbr.rescore_dual(
                    mol, sel_cids,
                    platform_name=gbsa_platform,
                    small_molecule_forcefield=gbsa_ff,
                    minimize=False,
                )
                gb_fields["gbsa_platform"] = water.platform
                gb_fields["gbsa_forcefield"] = gbsa_ff
                gb_fields["gbsa_seconds"] = round(time.time() - t_gb, 2)

                # Per-selected-conformer MMFF-geometry descriptors (no
                # GB minimization -> MMFF geometries preserved).
                sel_psa = np.asarray([psa_list[keep.index(c)] for c in sel_cids])
                sel_rg = np.asarray([rg_list[keep.index(c)] for c in sel_cids])
                sel_imhb = np.asarray(
                    [imhb_list[keep.index(c)] for c in sel_cids], dtype=float
                )

                e_w = np.asarray(water.energies_kcal)
                e_a = np.asarray(apolar.energies_kcal)
                # solvation free-energy DIFFERENCE per conformer.  Large
                # negative dE = strongly prefers water; large positive =
                # strongly prefers the apolar phase.
                dE = e_w - e_a

                # --- Q20 sub-ensembles by dE (robust to kT collapse) --- #
                # Top 20% by -dE (most water-loving) vs top 20% by +dE
                # (most apolar-loving).  This is the cleanest signal when
                # the GBSA energy spread dwarfs kT (which is the usual
                # case for large molecules).
                n = len(sel_cids)
                q_size = max(1, n // 5)
                order = np.argsort(dE)
                w_idx = order[:q_size]               # water-favoured
                a_idx = order[-q_size:]              # apolar-favoured
                psa_gb_po = float(sel_psa[w_idx].mean())
                psa_gb_ap = float(sel_psa[a_idx].mean())
                rg_gb_po = float(sel_rg[w_idx].mean())
                rg_gb_ap = float(sel_rg[a_idx].mean())
                imhb_gb_po = float(sel_imhb[w_idx].mean())
                imhb_gb_ap = float(sel_imhb[a_idx].mean())

                gb_fields.update(
                    psa3d_gb_apolar=psa_gb_ap, psa3d_gb_polar=psa_gb_po,
                    rg_gb_apolar=rg_gb_ap, rg_gb_polar=rg_gb_po,
                    imhb_gb_apolar=imhb_gb_ap, imhb_gb_polar=imhb_gb_po,
                )

                dPSA_gb = psa_gb_po - psa_gb_ap
                dRg_gb = rg_gb_po - rg_gb_ap
                dIMHB_gb = imhb_gb_ap - imhb_gb_po

                # GB-based CI.  Same orthogonal channel recipe as the MMFF
                # CI but driven by solvent-resolved deltas.  Weights tuned
                # so dPSA of ~20 Å² on a MW=1000 scaffold lands a value
                # comparable to the MMFF CI.
                term_psa_gb = max(0.0, dPSA_gb) / math.sqrt(max(mw, 1.0))
                if rg_gb_ap > 0:
                    term_rg_gb = math.log(
                        max(rg_gb_po / max(rg_gb_ap, 1e-6), 1.0)
                    )
                else:
                    term_rg_gb = 0.0
                term_imhb_gb = max(0.0, dIMHB_gb)
                ci_gb = (
                    3.0 * term_psa_gb
                    + 4.0 * term_rg_gb
                    + 1.0 * term_imhb_gb
                )
                if ci_gb >= 3.0:
                    verdict_gb = "chameleonic"
                elif ci_gb >= 1.5:
                    verdict_gb = "marginally chameleonic"
                else:
                    verdict_gb = "non-chameleonic"
                gb_fields.update(
                    dPSA_gb=float(dPSA_gb),
                    dRg_gb=float(dRg_gb),
                    dIMHB_gb=float(dIMHB_gb),
                    chameleonic_index_gb=float(ci_gb),
                    chameleonic_verdict_gb=verdict_gb,
                )
            except Exception as e:
                print(f"[{name}] GBSA rescore failed: {e}")
                import traceback
                traceback.print_exc()

    # export most-folded (lowest Rg) and most-extended (highest Rg) confs
    if sdf_out:
        try:
            w = Chem.SDWriter(sdf_out)
            iLow = int(np.argmin(rg))
            iHigh = int(np.argmax(rg))
            iPsaMin = int(np.argmin(psa))
            iPsaMax = int(np.argmax(psa))
            iImhbMax = int(np.argmax(imhb))
            tags = {
                "folded_minRg": keep[iLow],
                "extended_maxRg": keep[iHigh],
                "apolar_minPSA": keep[iPsaMin],
                "polar_maxPSA": keep[iPsaMax],
                "maxIMHB": keep[iImhbMax],
            }
            for tag, cid in tags.items():
                mm = Chem.Mol(mol)
                mm.RemoveAllConformers()
                mm.AddConformer(mol.GetConformer(cid), assignId=True)
                mm.SetProp("_Name", f"{name}_{tag}")
                mm.SetProp("tag", tag)
                mm.SetProp("rg", f"{rg[list(keep).index(cid)]:.3f}")
                mm.SetProp("psa3d", f"{psa[list(keep).index(cid)]:.3f}")
                mm.SetProp("imhb", f"{int(imhb[list(keep).index(cid)])}")
                w.write(mm)
            w.close()
        except Exception as e:
            print(f"[{name}] SDF export failed: {e}")

    tpsa = float(rdMolDescriptors.CalcTPSA(mol0))
    t1 = time.time()
    return Summary(
        name=name,
        smiles=smiles,
        mw=float(mw),
        heavy_atoms=int(heavy),
        tpsa_2d=tpsa,
        n_conf_embedded=len(cids),
        n_conf_unique=len(keep),
        time_seconds=round(t1 - t0, 2),
        psa3d_min=float(psa.min()),
        psa3d_max=float(psa.max()),
        psa3d_mean=float(psa.mean()),
        psa3d_std=float(psa.std()),
        rg_min=float(rg.min()),
        rg_max=float(rg.max()),
        rg_mean=float(rg.mean()),
        rg_std=float(rg.std()),
        imhb_mean=float(imhb.mean()),
        imhb_max=int(imhb.max()),
        imhb_frac_ge2=float(np.mean(imhb >= 2)),
        psa3d_apolar_wtd=psa_ap,
        psa3d_polar_wtd=psa_po,
        rg_apolar_wtd=rg_ap,
        rg_polar_wtd=rg_po,
        imhb_apolar_wtd=imhb_ap,
        imhb_polar_wtd=imhb_po,
        psa3d_q_apolar=psa_q_ap,
        psa3d_q_polar=psa_q_po,
        rg_q_apolar=rg_q_ap,
        rg_q_polar=rg_q_po,
        imhb_q_apolar=imhb_q_ap,
        imhb_q_polar=imhb_q_po,
        dPSA=float(dPSA),
        dRg=float(dRg),
        dIMHB=float(dIMHB),
        chameleonic_index=float(ci),
        chameleonic_verdict=verdict,
        **gb_fields,
    )


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def print_summary(s: Summary) -> None:
    print()
    print(f"==== {s.name}  MW={s.mw:.1f}  heavy={s.heavy_atoms}  "
          f"TPSA(2D)={s.tpsa_2d:.1f}  unique_confs={s.n_conf_unique}  "
          f"time={s.time_seconds}s ====")
    print(f"  3D PSA (A^2)   min={s.psa3d_min:6.1f}  max={s.psa3d_max:6.1f}  "
          f"mean={s.psa3d_mean:6.1f}  std={s.psa3d_std:5.1f}  "
          f"range={s.psa3d_max - s.psa3d_min:5.1f}")
    print(f"  Rg     (A)     min={s.rg_min:5.2f}  max={s.rg_max:5.2f}  "
          f"mean={s.rg_mean:5.2f}  std={s.rg_std:4.2f}")
    print(f"  IMHB           mean={s.imhb_mean:4.2f}  max={s.imhb_max}  "
          f"frac>=2={s.imhb_frac_ge2:.2f}")
    print(f"  apolar-Boltz:  <PSA>={s.psa3d_apolar_wtd:6.1f}  "
          f"<Rg>={s.rg_apolar_wtd:5.2f}  <IMHB>={s.imhb_apolar_wtd:4.2f}")
    print(f"  polar -Boltz:  <PSA>={s.psa3d_polar_wtd:6.1f}  "
          f"<Rg>={s.rg_polar_wtd:5.2f}  <IMHB>={s.imhb_polar_wtd:4.2f}")
    print(f"  apolar-Q20:    <PSA>={s.psa3d_q_apolar:6.1f}  "
          f"<Rg>={s.rg_q_apolar:5.2f}  <IMHB>={s.imhb_q_apolar:4.2f}")
    print(f"  polar -Q20:    <PSA>={s.psa3d_q_polar:6.1f}  "
          f"<Rg>={s.rg_q_polar:5.2f}  <IMHB>={s.imhb_q_polar:4.2f}")
    print(f"  dPSA={s.dPSA:+.1f}  dRg={s.dRg:+.2f}  dIMHB={s.dIMHB:+.2f}  "
          f"CI={s.chameleonic_index:+.2f}  -> {s.chameleonic_verdict}")
    if s.chameleonic_index_gb is not None:
        print(f"  GBSA ({s.gbsa_platform}, {s.gbsa_forcefield}, "
              f"{s.gbsa_seconds}s):")
        print(f"    apolar (eps=2.02): <PSA>={s.psa3d_gb_apolar:6.1f}  "
              f"<Rg>={s.rg_gb_apolar:5.2f}  "
              f"<IMHB>={s.imhb_gb_apolar:4.2f}")
        print(f"    polar  (eps=78.5): <PSA>={s.psa3d_gb_polar:6.1f}  "
              f"<Rg>={s.rg_gb_polar:5.2f}  "
              f"<IMHB>={s.imhb_gb_polar:4.2f}")
        print(f"    dPSA_gb={s.dPSA_gb:+.1f}  dRg_gb={s.dRg_gb:+.2f}  "
              f"dIMHB_gb={s.dIMHB_gb:+.2f}  "
              f"CI_gb={s.chameleonic_index_gb:+.2f}  "
              f"-> {s.chameleonic_verdict_gb}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("smiles", nargs="?", default=None)
    ap.add_argument("--name", default="mol")
    ap.add_argument("--nconf", type=int, default=300)
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--batch", default=None,
                    help="TSV file: NAME\\tSMILES per line")
    ap.add_argument("--json", default=None, help="write summaries here")
    ap.add_argument("--sdf-dir", default=None,
                    help="directory to write per-molecule SDF snapshots")
    ap.add_argument("--gbsa", action="store_true",
                    help="rescore the Butina-kept conformer set with "
                         "OpenMM OBC2 GBSA in water and chloroform-like "
                         "dielectrics (requires openmm + openff-toolkit + "
                         "openmmforcefields)")
    ap.add_argument("--gbsa-max-conf", type=int, default=64,
                    help="cap on #confs handed to GBSA rescore")
    ap.add_argument("--gbsa-platform", default="CUDA",
                    help="OpenMM platform preference: CUDA / OpenCL / CPU")
    ap.add_argument("--gbsa-ff", default="openff-2.1.0",
                    help="small-molecule force field for SystemGenerator")
    args = ap.parse_args()

    items = []
    if args.batch:
        with open(args.batch) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if "\t" in line:
                    name, smi = line.split("\t", 1)
                else:
                    parts = line.split(None, 1)
                    if len(parts) == 2:
                        name, smi = parts
                    else:
                        name, smi = "mol", parts[0]
                items.append((name.strip(), smi.strip()))
    elif args.smiles:
        items.append((args.name, args.smiles))
    else:
        ap.error("provide SMILES or --batch file")

    import os
    if args.sdf_dir:
        os.makedirs(args.sdf_dir, exist_ok=True)

    summaries = []
    for name, smi in items:
        sdf_path = (os.path.join(args.sdf_dir, f"{name}_snapshots.sdf")
                    if args.sdf_dir else None)
        s = summarize(
            name, smi,
            n_conf=args.nconf, n_threads=args.threads,
            sdf_out=sdf_path,
            gbsa=args.gbsa,
            gbsa_max_conf=args.gbsa_max_conf,
            gbsa_platform=args.gbsa_platform,
            gbsa_ff=args.gbsa_ff,
        )
        print_summary(s)
        summaries.append(asdict(s))

    if args.json:
        with open(args.json, "w") as f:
            json.dump(summaries, f, indent=2)


if __name__ == "__main__":
    main()
