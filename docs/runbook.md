# Runbook for the public analysis subset

Keep all scientific inputs and generated outputs **outside this checkout**. Commands below are public analysis entrypoints, not the private end-to-end launchers. Private design inputs and historical orchestration scripts are intentionally absent.

## Analysis environment

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-analysis.txt
```

Use separate licensed/model-specific environments for PyRosetta and inference software. Do not upgrade a shared working installation just to run this subset. Save version inventories for reproducibility.

## Input contracts

- PDB audit: standard amino acids, complete heavy atoms, explicit chain IDs and residue numbering.
- Retained antibody geometry: H/L chains with the campaign hu4D5 canonical cysteine numbering; default CDR bounds are campaign-specific, not universal IMGT coordinates. Re-annotate before another framework.
- TNF packing: antigen chains A/B/C and antibody H/L; compact-binder support uses chain D. Set `TNF_H3_LENGTH` from your own numbered antibody. Supply complete contiguous residues expected by the CDR definitions.
- pKa wrapper: explicitly specify binder and target chains; all present chains must be assigned. Heteroatom/glycan handling requires separate review.
- Sequence flags: a private CSV with columns `name,VH,VL`; no such file is supplied.

```bash
python tools/sequence_flags.py --input /path/outside/checkout/input.csv --output /path/outside/checkout/flags.csv
python tools/audit_structure.py --input /path/outside/checkout/complex.pdb --output /path/outside/checkout/audit.json
python tools/ph_linkage.py --input /path/outside/checkout/complex.pdb --binder-chains H L --target-chains A B C --output /path/outside/checkout/pka
```

In a separately authorized PyRosetta environment:

```bash
export TNF_H3_LENGTH=17  # replace with your own annotated loop length
python tools/pack_interface.py --input /path/outside/checkout/complex.pdb --output /path/outside/checkout/packing
```

The value above is an illustrative configuration, not a published design. `--sequences` accepts a private H/L mapping for local mutation/repacking; the mapping is not included. After packing, independently rerun geometry, intended-site, glycan and context checks. The packing score gate alone does not certify a binder.

## GPU and inference workflow

```bash
nvidia-smi --query-gpu=index,name,memory.used,memory.total,utilization.gpu --format=csv
nvidia-smi --query-compute-apps=gpu_uuid,pid,used_memory --format=csv
```

Check availability immediately before each launch and use scheduler/reservation mechanisms on shared systems. The supplied GPU helper is a snapshot check, not an atomic reservation. Never terminate another user's job.

For generation, sequence design and independent prediction use the installed upstream CLI and its version-specific input schema. Record chain formatting, design masks, target context, hotspots, model versions and seeds. This public repository deliberately omits private launch manifests; no generic one-line command can reproduce the withheld campaign. Follow the stage gates in workflow.md rather than assuming all generated outputs deserve expensive prediction.

## Publication checks

```bash
python tools/publication_guard.py
python -m unittest discover -s tests -v
# Install the local guard for future commits and pushes:
git config core.hooksPath hooks
```

The guard checks staged files and names; it supplements human review and cannot prove the absence of every short or encoded sequence. Never force-add design files. The initial publication was also checked locally against withheld candidate sequences without publishing that comparison database.
