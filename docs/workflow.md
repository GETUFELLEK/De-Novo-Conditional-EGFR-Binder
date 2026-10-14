# Shared computational workflow

## 1. Target and framework preparation

Record the target accession, experimental/structural source, biological assembly, species and construct boundaries. Map every prepared residue back to source and UniProt numbering. Preserve resolved glycans and mark unresolved glycosylation anchors. A truncated local target is useful for generation only if the complete biological context is restored before advancement.

Number antibody variable domains with IMGT/ANARCI. Keep framework/CDR boundaries explicit. Audit all consecutive C–N distances, bond angles, chain continuity, canonical disulfides, framework packing and self-clashes. Check compressed peptide bonds as well as long breaks. A known antibody structure may serve as a calibration reference; framework reuse must be distinguished from copying its parental CDRs or docking pose.

## 2. Backbone generation

RFantibody/RFdiffusion and RFD3 were explored in bounded branches. Vary approach direction, rotation, docking depth, CDR participation and hotspot chemistry rather than merely increasing seeds for a failed pose family. Use target sidechain and glycan context during generation where supported; otherwise restore them immediately before screening. RFD3 all-atom capability does not itself establish a usable interface.

## 3. Precursor gates and sequence design

Reject broken covalent geometry, major framework deformation, severe clashes, inappropriate target surfaces and incompatible biological context. Evaluate hotspot realization and CDR participation. Distinguish a repairable sidechain/packing problem from an invalid backbone, and document every repair. A backbone that appears clean with antigen sidechains omitted is not a clean precursor.

Apply AbMPNN to eligible antibody precursors with explicit designable CDR positions and fixed framework/target positions. Verify masks and exact output sequences. Evaluate multiple sequence alternatives, introduced cysteines, glycosylation sequons and other liabilities before expensive prediction. Generic ProteinMPNN and antibody-aware MPNN are not interchangeable labels.

## 4. Local packing and independent prediction

Local repacking tests compatibility with an assumed pose. It does not demonstrate that the sequence adopts that pose. Use AF3/ESMFold2 sequence-only predictions to test recovery; use other predictors selectively to diagnose disagreement. First sample broadly across qualified sequences, then deepen sampling of justified survivors. Track unique sequences separately from predicted models.

Align targets before evaluating antibody pose, contact footprints and residue pairs. Quantify per-CDR contacts, precursor retention, confidence, protein clashes, glycan compatibility and recurrence. Retain failed predictions in the evidence record. Confidence scores support geometric evidence; they do not replace it.

## 5. Cross-species and pH hypotheses

Test the exact same H/L pair against the homologous target. Map contact partners and compare equivalent geometry, rather than merely counting contacts in each species. After a credible parent interaction is available, prioritize titratable interface sites with a physical mechanism. Use bounded mutations and matched neutral controls. Histidine count alone is not evidence of a pH switch.

## 6. Ranking and reporting

Maintain sequence provenance, physical status, predictor evidence, novelty and uncertainty as separate fields. Annotate whether evidence is precursor-based, locally conditioned, independently predicted or experimental. Competition design intent is not a substitute for validated binding. Publication and submission rules must be checked at the time of use. No private candidate evidence is distributed here.
