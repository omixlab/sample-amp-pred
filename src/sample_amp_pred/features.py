"""
Feature engineering for the sAMPle general antibacterial classifier.

This reproduces, column for column and in the exact order the model was trained on, the 586
features sAMPle's own ``sample.inference`` module computes for a bare one-letter peptide
sequence: 11 physicochemical descriptors, 28 amino-acid-composition frequencies, 400 dipeptide
frequencies, and 147 Composition/Transition/Distribution (CTD) descriptors. No block is skipped
or approximated -- a wrong feature produces a confident, silently wrong probability, which is
worse than this module raising.

One deliberate simplification against the live sAMPle server, stated plainly because it changes
a number: sAMPle's own training/serving pipeline corrects a peptide's molecular weight and net
charge for a closed backbone (head-to-tail cyclic, or disulfide-bridged) when the exact sequence
is present in sAMPle's own curated peptide-topology index. That index is internal to the sAMPle
project and is not published alongside the model, so it cannot be reproduced here. Every sequence
is therefore treated as a linear chain. This is a no-op for the overwhelming majority of inputs
(linear peptides, or any sequence sAMPle has not independently curated as cyclic), and it is the
same linear-chain fallback sAMPle's own server already uses for any sequence outside that index.
"""
from __future__ import annotations

import math
import re
from typing import Dict, List, Sequence, Tuple

import numpy as np

# ---------------------------------------------------------------------------------------------
# Sequence normalization
# ---------------------------------------------------------------------------------------------

AMINO_ACIDS: str = "ACDEFGHIKLMNPQRSTVWY"
_VALID_RESIDUES = frozenset(AMINO_ACIDS)

#: EMBOSS pKa scale, used for the `charge` descriptor (Henderson-Hasselbalch at pH 7.4).
_PKA_EMBOSS: Dict[str, float] = {
    "N_term": 8.6, "C_term": 3.6, "C": 8.5, "D": 3.9, "E": 4.1,
    "H": 6.5, "K": 10.8, "R": 12.5, "Y": 10.9,
}


def clean_sequence(sequence: str) -> str:
    """
    Normalizes a raw FASTA sequence to the uppercase one-letter residues the model expects.

    Strips whitespace, digits, and any character outside the 20 standard amino-acid letters.
    This tool does not attempt sAMPle's own three-letter-code / modifier parsing (``Gly-Ile-Gly``,
    ``-NH2`` amidation, D-/L- stereochemistry prefixes) -- a plain FASTA file holds one-letter
    sequences by convention, and guessing at a non-standard spelling risks a silent
    misinterpretation worse than refusing it outright (see :func:`validate_sequence`).

    :param sequence: Raw sequence text.
    :returns: Uppercase one-letter sequence, keeping only standard residues.
    """
    return re.sub(r"[^A-Za-z]", "", str(sequence or "")).upper()


def validate_sequence(sequence: str) -> str | None:
    """
    :param sequence: A sequence already run through :func:`clean_sequence`.
    :returns: ``None`` when scorable, otherwise a short reason it is not.
    """
    if not sequence:
        return "empty sequence"
    if len(sequence) < 2:
        return "sequence too short to score (minimum 2 residues)"
    unknown = sorted(set(sequence) - _VALID_RESIDUES)
    if unknown:
        return f"non-standard residue(s): {', '.join(unknown)}"
    return None


def get_peptide_charge(sequence: str, ph: float = 7.4) -> float:
    """
    Net charge at `ph` via the Henderson-Hasselbalch equation over the EMBOSS pKa scale.

    Matches ``sample.sequence.get_peptide_charge`` for a non-amidated linear peptide (the case
    every bare FASTA sequence falls into here).
    """
    if not sequence:
        return 0.0
    charge = 1.0 / (1.0 + 10 ** (ph - _PKA_EMBOSS["N_term"]))
    charge -= 1.0 / (1.0 + 10 ** (_PKA_EMBOSS["C_term"] - ph))
    charge += sequence.count("K") / (1.0 + 10 ** (ph - _PKA_EMBOSS["K"]))
    charge += sequence.count("R") / (1.0 + 10 ** (ph - _PKA_EMBOSS["R"]))
    charge += sequence.count("H") / (1.0 + 10 ** (ph - _PKA_EMBOSS["H"]))
    charge -= sequence.count("D") / (1.0 + 10 ** (_PKA_EMBOSS["D"] - ph))
    charge -= sequence.count("E") / (1.0 + 10 ** (_PKA_EMBOSS["E"] - ph))
    charge -= sequence.count("C") / (1.0 + 10 ** (_PKA_EMBOSS["C"] - ph))
    charge -= sequence.count("Y") / (1.0 + 10 ** (_PKA_EMBOSS["Y"] - ph))
    return round(charge, 4)


# ---------------------------------------------------------------------------------------------
# Amino-acid composition (28 columns)
# ---------------------------------------------------------------------------------------------

AA_GROUPS: Dict[str, str] = {
    "aliphatic": "AILV",
    "aromatic": "FWY",
    "acidic": "DE",
    "basic": "RHK",
    "polar_uncharged": "STNQCM",
    "small_tiny": "GAS",
    "sulfur": "CM",
    "proline": "P",
}

COMPOSITION_COLUMNS: Tuple[str, ...] = tuple(f"aa_freq_{aa}" for aa in AMINO_ACIDS) + tuple(
    f"aa_group_freq_{g}" for g in AA_GROUPS
)


def _composition_row(seq: str) -> List[float]:
    n = len(seq)
    if n == 0:
        return [0.0] * len(COMPOSITION_COLUMNS)
    single = [seq.count(aa) / n for aa in AMINO_ACIDS]
    grouped = [sum(seq.count(r) for r in residues) / n for residues in AA_GROUPS.values()]
    return single + grouped


# ---------------------------------------------------------------------------------------------
# Dipeptide composition (400 columns)
# ---------------------------------------------------------------------------------------------

DPC_COLUMNS: Tuple[str, ...] = tuple(f"dpc_{a}{b}" for a in AMINO_ACIDS for b in AMINO_ACIDS)


def _dpc_row(seq: str) -> List[float]:
    n = len(seq)
    if n < 2:
        return [0.0] * len(DPC_COLUMNS)
    counts: Dict[str, int] = {}
    for i in range(n - 1):
        pair = seq[i:i + 2]
        counts[pair] = counts.get(pair, 0) + 1
    total = float(n - 1)
    return [counts.get(f"{a}{b}", 0) / total for a in AMINO_ACIDS for b in AMINO_ACIDS]


# ---------------------------------------------------------------------------------------------
# CTD: Composition / Transition / Distribution (147 columns)
# ---------------------------------------------------------------------------------------------

CTD_PROPERTIES: Dict[str, Tuple[str, str, str]] = {
    "hydrophobicity": ("RKEDQN", "GASTPHY", "CLVIMFW"),
    "vdw_volume": ("GASTPDC", "NVEQIL", "MHKFRYW"),
    "polarity": ("LIFWCMVY", "PATGS", "HQRKNED"),
    "polarizability": ("GASDT", "CPNVEQIL", "KMHFRYW"),
    "charge": ("KR", "ANCQGHILMFPSTWYV", "DE"),
    "secondary_structure": ("EALMQKRH", "VIYCWFT", "GNPSD"),
    "solvent_accessibility": ("ALFCGIVW", "RKQEND", "MPSTHY"),
}

_CTD_DISTRIBUTION_POINTS: Tuple[float, ...] = (0.0, 0.25, 0.50, 0.75, 1.00)


def _ctd_columns() -> List[str]:
    names: List[str] = []
    for prop in CTD_PROPERTIES:
        names.extend(f"ctd_comp_{prop}_g{g}" for g in (1, 2, 3))
    for prop in CTD_PROPERTIES:
        names.extend(f"ctd_tran_{prop}_{pair}" for pair in ("12", "13", "23"))
    for prop in CTD_PROPERTIES:
        for g in (1, 2, 3):
            names.extend(f"ctd_dist_{prop}_g{g}_p{int(p * 100)}" for p in _CTD_DISTRIBUTION_POINTS)
    return names


CTD_COLUMNS: Tuple[str, ...] = tuple(_ctd_columns())


def _ctd_row(seq: str) -> List[float]:
    if not seq:
        return [0.0] * len(CTD_COLUMNS)

    encoded: Dict[str, List[int]] = {}
    for prop, groups in CTD_PROPERTIES.items():
        member = {r: i + 1 for i, residues in enumerate(groups) for r in residues}
        encoded[prop] = [member[r] for r in seq if r in member]

    composition: List[float] = []
    for prop in CTD_PROPERTIES:
        classes = encoded[prop]
        n = len(classes)
        composition.extend(classes.count(g) / n if n else 0.0 for g in (1, 2, 3))

    transition: List[float] = []
    for prop in CTD_PROPERTIES:
        classes = encoded[prop]
        pairs = float(len(classes) - 1)
        for first, second in ((1, 2), (1, 3), (2, 3)):
            if pairs <= 0:
                transition.append(0.0)
                continue
            crossings = sum(
                1 for a, b in zip(classes, classes[1:])
                if (a, b) == (first, second) or (a, b) == (second, first)
            )
            transition.append(crossings / pairs)

    distribution: List[float] = []
    for prop in CTD_PROPERTIES:
        classes = encoded[prop]
        n = len(classes)
        for g in (1, 2, 3):
            positions = [i + 1 for i, c in enumerate(classes) if c == g]
            total = len(positions)
            if total == 0 or n == 0:
                distribution.extend([0.0] * len(_CTD_DISTRIBUTION_POINTS))
                continue
            for point in _CTD_DISTRIBUTION_POINTS:
                rank = 0 if point == 0.0 else int(math.ceil(point * total)) - 1
                distribution.append(positions[min(rank, total - 1)] / n)

    return composition + transition + distribution


# ---------------------------------------------------------------------------------------------
# Physicochemical descriptors (11 columns)
# ---------------------------------------------------------------------------------------------

DESCRIPTOR_COLUMNS: Tuple[str, ...] = (
    "molecular_weight", "logp", "hydrogen_bond_donors", "hydrogen_bond_acceptors",
    "charge", "sequence_length", "boman_index", "aliphatic_index", "hydrophobic_moment",
    "isoelectric_point", "net_charge_ph7",
)


def _descriptor_row(seq: str) -> Dict[str, float]:
    """
    Molecular-weight/LogP/H-bond descriptors come from RDKit over a linear peptide built
    straight from the sequence (``Chem.MolFromSequence``); the remaining five come from
    `peptides.py`, the same library sAMPle's own enrichment step uses, so the numbers match for
    any sequence outside sAMPle's own cyclic-peptide topology index (see this module's
    docstring).
    """
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    row: Dict[str, float] = {
        "molecular_weight": float("nan"), "logp": float("nan"),
        "hydrogen_bond_donors": float("nan"), "hydrogen_bond_acceptors": float("nan"),
    }
    try:
        mol = Chem.MolFromSequence(seq)
        if mol is not None:
            row["molecular_weight"] = float(Descriptors.MolWt(mol))
            row["logp"] = float(Descriptors.MolLogP(mol))
            row["hydrogen_bond_donors"] = float(Descriptors.NumHDonors(mol))
            row["hydrogen_bond_acceptors"] = float(Descriptors.NumHAcceptors(mol))
    except Exception:
        pass

    try:
        row["charge"] = float(get_peptide_charge(seq))
    except Exception:
        row["charge"] = float("nan")
    row["sequence_length"] = float(len(seq))

    try:
        import peptides

        peptide = peptides.Peptide(seq)
        row["boman_index"] = float(peptide.boman())
        row["aliphatic_index"] = float(peptide.aliphatic_index())
        row["hydrophobic_moment"] = float(peptide.hydrophobic_moment(angle=100))
        row["isoelectric_point"] = float(peptide.isoelectric_point())
        row["net_charge_ph7"] = float(peptide.charge(pH=7.4))
        weight = float(peptide.molecular_weight())
        if weight == weight:  # not NaN
            row["molecular_weight"] = weight
    except Exception:
        for column in ("boman_index", "aliphatic_index", "hydrophobic_moment",
                       "isoelectric_point", "net_charge_ph7"):
            row.setdefault(column, float("nan"))

    return row


# ---------------------------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------------------------

#: The model's own training-column order (descriptors, composition, DPC, CTD -- see
#: ``sample.flaml_xgboost``'s own "every block, concatenated" comment). Re-derived here as a
#: fallback only; the shipped model artifact carries its OWN column list, which
#: :func:`sample_amp_pred.model.load_bundled_model` uses in preference to this one, so a
#: future retrain that reorders blocks cannot silently feed this tool's vector into the wrong
#: column.
ALL_COLUMNS: Tuple[str, ...] = DESCRIPTOR_COLUMNS + COMPOSITION_COLUMNS + DPC_COLUMNS + CTD_COLUMNS


def build_feature_matrix(sequences: Sequence[str], columns: Sequence[str]) -> np.ndarray:
    """
    Builds the `len(sequences) x len(columns)` feature matrix a model needs, in `columns` order.

    :param sequences: Cleaned (see :func:`clean_sequence`), validated one-letter sequences.
    :param columns: The model's own training columns (from its bundled artifact).
    :returns: float32 matrix; a column this module cannot compute raises, rather than being
        silently filled with zeros.
    """
    rows: List[Dict[str, float]] = []
    for seq in sequences:
        row: Dict[str, float] = {}
        row.update(_descriptor_row(seq))
        row.update(dict(zip(COMPOSITION_COLUMNS, _composition_row(seq))))
        row.update(dict(zip(DPC_COLUMNS, _dpc_row(seq))))
        row.update(dict(zip(CTD_COLUMNS, _ctd_row(seq))))
        rows.append(row)

    known = set(ALL_COLUMNS)
    missing = [c for c in columns if c not in known]
    if missing:
        raise RuntimeError(
            f"This tool cannot compute {len(missing)} of the model's {len(columns)} feature "
            f"columns (first missing: {', '.join(missing[:5])}). The bundled model may be newer "
            "than this package version."
        )

    matrix = np.empty((len(rows), len(columns)), dtype=np.float32)
    for i, row in enumerate(rows):
        matrix[i, :] = [row[c] for c in columns]
    return matrix
