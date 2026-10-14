# Evidence interpretation

| Evidence | Useful interpretation | Does not establish |
|---|---|---|
| Clean generated interface | Plausible geometry under generator constraints | Sequence encodability |
| Exact-sequence local repacking | Compatibility with the assumed pose | Independent binding recovery |
| Sequence-only prediction | Model support for a sequence-dependent pose | Experimental binding |
| Multiple samples | Recurrence under a specified model/configuration | Independent experimental evidence |
| Predictor agreement | Reduced model disagreement | Certainty or measured affinity |
| High ipTM / low interface PAE | Confidence supporting other evidence | Correct epitope or physical validity alone |
| Rosetta interface energy | Local model scoring/packing comparison | Measured free energy or affinity |
| Bound/free pKa linkage | A conditional-binding hypothesis | Quantitative pH selectivity |
| Sequence novelty screen | Distance from specified references | Official eligibility or biological performance |

For antibody novelty, record IMGT-numbered VH/VL, CDRH3, nearest H3 match and metric, global chain similarities, reference source/version and uncertainty. Do not substitute generic fold novelty for antibody sequence/CDR novelty. Novelty classification and competition eligibility remain separate. Public rule frameworks may change; this repository does not encode an authoritative eligibility decision.

Known native controls may be represented in training data. Keep failures, unmatched chains, incomplete jobs and missing confidence metrics explicit. Do not report absent metrics as zero. Unresolved glycans, protonation coupling, uncertain rotamers, reference-database coverage and predictor training exposure limit interpretation.
