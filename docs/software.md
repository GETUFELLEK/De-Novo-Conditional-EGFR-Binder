# Software roles and provenance boundaries

| Component | Campaign role | Redistribution |
|---|---|---|
| RFantibody / RFdiffusion | Antibody backbone generation | Upstream source and weights required separately |
| RFdiffusion3 (RFD3) | Alternative generation with all-atom context | Upstream installation required separately |
| AbMPNN / ProteinMPNN | Antibody-aware / general sequence design | Separate tools and checkpoints |
| AlphaFold3 | Independent complex prediction | Authorized software, data and parameters required |
| ESMFold2 | Orthogonal complex prediction | Separate model installation and weights |
| RoseTTAFold-family predictors | Selected calibration/orthogonal checks | Separate installation |
| PyRosetta | Local packing and interface metrics | Licensed dependency, not redistributed |
| PROPKA | Structure-conditioned pKa estimates | Install separately |
| ANARCI / IMGT | Antibody numbering and CDR assignment | Numbering dependencies and references required |
| PLAbDab / TheraSAbDab | Reference-limited novelty screening | Databases are not redistributed |
| Biopython / NumPy / SciPy | Coordinates, alignment and numerical analysis | Minimal dependencies listed |

The release includes selected project-authored utilities, not upstream model implementations. `tools/geometry.py`, `tools/pack_interface.py` and `tools/gpu_resources.py` are retained from the campaign analysis code. `tools/ph_linkage.py` extracts the campaign's rigid linkage calculation into a configurable CLI. `tools/audit_structure.py` and `tools/sequence_flags.py` expose reusable checks without private project inputs.

Record source commits, model-weight hashes, Python/CUDA versions, container digests, environment inventories, seeds and full command lines for each new run. This repository does not claim a tested clean-machine reinstall of the historical inference stack. The analysis dependency list is not a complete historical lock. No blanket license for third-party software or data is granted here.
