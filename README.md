# De-Novo-Conditional-EGFR-Binder

Methods and selected analysis tools from computational conditional-antibody design campaigns against **EGFR** and **TNF-alpha**. The repository name reflects the first antigen; the TNF-alpha workflow is documented separately.

**This is a methods-only release. No de novo antibody designs are included:** no designed sequences, CDR sequences, mutation panels, candidate identifiers, structures, prediction outputs, ranking tables, submission packages, or private colleague data. Inputs and generated outputs belong outside this checkout.

## Contents

- [Shared workflow](docs/workflow.md): target preparation, generation, sequence design and staged evaluation.
- [EGFR methods](docs/egfr.md): domain/epitope mapping, glycan-aware interface analysis and cross-species evaluation.
- [TNF-alpha methods](docs/tnf_alpha.md): trimer geometry, C3 symmetry and neutral-ON/acid-OFF hypotheses.
- [Runbook](docs/runbook.md): setup, private inputs and executable analysis commands.
- [Evidence interpretation](docs/evidence.md): confidence, physical screening, novelty and limitations.
- [Software provenance](docs/software.md): tool roles and reproducibility boundaries.
- [Publication policy](DATA_POLICY.md): safeguards against publishing designs.
- `tools/`: selected campaign analysis utilities and a standalone pKa wrapper.

## Scope

The campaigns explored RFantibody/RFdiffusion and RFdiffusion3 generation, antibody-aware MPNN sequence optimization, AF3 and ESMFold2 evaluation, and selected additional structural checks. RFD3 is a generator; RF3 is a predictor. Tool use varied by branch; this does not imply that every candidate passed every method.

Target-relative geometry, all-atom physical compatibility and sequence provenance were evaluated alongside confidence scores. ipTM alone was never sufficient to establish intended-site binding. Known-complex controls validate pipeline operation but may benefit from training familiarity. Conditioned local models, independent predictions and experimental evidence are distinct.

This is **not** the complete private campaign repository or a turnkey reproduction of the final designs. Model software, licensed dependencies, weights, reference databases and private inputs are not redistributed. No experimentally validated affinity, cross-reactivity or pH selectivity is claimed.

## Quick start

Use a separate analysis environment and external input/output directories. See the runbook before applying the tools: the retained geometry and packing utilities have campaign-specific chain/numbering assumptions.

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-analysis.txt
python tools/sequence_flags.py --help
python tools/audit_structure.py --help
python tools/ph_linkage.py --help
```

The optional packing script requires a separately licensed PyRosetta installation. Its full runtime was not reinstalled or exercised for this publication.
