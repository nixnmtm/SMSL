# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)


def _alpha_code(i):
    """0 -> A, 1 -> B, ..., 25 -> Z, 26 -> AA, ..."""
    letters = []
    i = int(i)
    while True:
        i, rem = divmod(i, 26)
        letters.append(chr(ord("A") + rem))
        if i == 0:
            break
        i -= 1
    return "".join(reversed(letters))


def _is_psf_safe_segid(segid, maxlen=8):
    if segid is None:
        return False
    segid = str(segid).strip()
    if not segid:
        return False
    if " " in segid:
        return False
    if len(segid) > maxlen:
        return False
    return True


def _first_molnum(bead):
    try:
        molnums = bead.atoms.molnums
        if molnums is not None and len(molnums) > 0:
            return int(molnums[0])
    except Exception:
        pass
    return None


def _first_chainid(bead):
    try:
        chainids = bead.atoms.chainIDs
        if chainids is not None and len(chainids) > 0:
            chainid = str(chainids[0]).strip()
            if chainid:
                return chainid
    except Exception:
        pass
    return None


def _first_fragment(bead):
    try:
        fragindices = bead.atoms.fragindices
        if fragindices is not None and len(fragindices) > 0:
            return int(fragindices[0])
    except Exception:
        pass

    try:
        fragindex = bead.atoms[0].fragment.fragindex
        return int(fragindex)
    except Exception:
        pass

    return None


def _first_existing_segid(bead):
    try:
        segids = bead.atoms.segids
        if segids is not None and len(segids) > 0:
            segid = str(segids[0]).strip()
            if segid:
                return segid
    except Exception:
        pass
    return None


def _first_resname(bead):
    try:
        resname = str(bead.resnames[0]).strip().upper()
        if resname:
            return resname
    except Exception:
        pass
    return None


def _protein_group_key(bead):
    """
    For protein beads, prefer fragment first because GRO/TPR may not have
    useful segids, and fragment separates protein bodies in this workflow.

    Fallback order:
        fragment -> molnum -> chainID -> existing segid -> generic protein
    """
    frag = _first_fragment(bead)
    if frag is not None:
        return ("frag", frag)

    molnum = _first_molnum(bead)
    if molnum is not None:
        return ("molnum", molnum)

    chainid = _first_chainid(bead)
    if chainid and len(chainid) <= 4 and " " not in chainid:
        return ("chain", chainid)

    segid = _first_existing_segid(bead)
    if _is_psf_safe_segid(segid) and segid.upper() not in {"SYS", "SYSTEM", "PROT"}:
        return ("segid", segid)

    return ("fallback", "protein")


def _safe_segid(name, bead, user_segid_map=None, protein_group_map=None, default="SYS"):
    """
    Stable segid assignment.

    Priority:
      1. explicit user override by residue name
      2. preserve existing valid segid
      3. protein bead -> fragment/molnum/chain grouped PROA/PROB/...
      4. fallback to residue name
      5. fallback default
    """
    if user_segid_map is None:
        user_segid_map = {}
    if protein_group_map is None:
        protein_group_map = {}

    resname = _first_resname(bead)

    # 1. explicit user override wins
    if resname in user_segid_map:
        return str(user_segid_map[resname])[:8]

    # 2. preserve existing safe segid if present
    raw = _first_existing_segid(bead)
    if _is_psf_safe_segid(raw) and raw.upper() not in {
        "SYSTEM",
        "PROTEININWATER",
        "PROTEIN IN WATER",
    }:
        return raw[:8]

    # 3. protein grouping
    if name in {"CA", "CB", "N", "O"}:
        key = _protein_group_key(bead)
        if key not in protein_group_map:
            protein_group_map[key] = "PRO{}".format(_alpha_code(len(protein_group_map)))
        return protein_group_map[key][:8]

    # 4. residue-name fallback for ions/ligands/other non-protein beads
    if resname:
        return resname[:8]

    # 5. final fallback
    return default[:8]