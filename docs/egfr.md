# EGFR workflow and lessons

The EGFR work focused on antibody interactions with extracellular EGFR, including domain III, and subsequent human/mouse and conditional-binding evaluation.

1. Map prepared target numbering to the full-length human and mouse sequences. Retain each branch's original intended epitope; scoring every design against one common patch can misclassify success.
2. Track a per-sequence survival funnel from RF/RFD3 precursor geometry through sequence design and independent predictions. Count both epitope contacts and target-aligned orientation.
3. Restore complete antigen sidechains, resolved glycans and the full ectodomain context before accepting locally generated interfaces. Unresolved glycan anchors warrant explicit uncertainty.
4. Use native complexes to calibrate chain mapping, contact analysis and predictor execution. Native recovery may reflect model familiarity and is not evidence of prospective design accuracy.
5. Audit frameworks independently of their original antigen specificity. Framework histidines and parental CDR histidines have different provenance; parental CDR identities and poses must not silently become design inputs.
6. Compare sequence-only predictors and diagnose disagreements. A high confidence score on an alternative EGFR surface does not establish recovery of the intended epitope.
7. Require exact-sequence human/mouse comparison and an explicit interfacial mechanism for conditional-binding variants. Local energy and pKa calculations remain model-dependent hypotheses.

The private campaign included framework calibration, original-site diagnostics, orthogonal prediction, local repair and sequence-optimization branches. The public release explains the decision process without exposing the selected designs, exact mutation choices, poses or candidate-level results.
