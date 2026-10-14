# Publication boundary

This repository contains methods, documentation and reviewed reusable source code only. De novo design information remains private.

Do not commit or push:

- Designed VH/VL/CDR sequences, mutation panels or candidate identifiers.
- Generated, optimized, docked or predicted binder coordinates, including coordinates embedded in scripts or notebooks.
- Candidate registries, ranking/evidence tables, input manifests, logs, submission CSVs or archives.
- Colleague sequences, structures, screenshots or private communications.
- Model parameters, databases, environments, credentials or private infrastructure configuration.

A `.py` or `.md` extension does not make embedded design information safe to publish. Review string literals and documentation as carefully as data files. Use an independent repository with a fresh history; do not add the full working campaign tree or its existing history.

`.gitignore` uses a default-deny layout. The local pre-commit/pre-push guard enforces permitted paths and flags obvious sequence literals, candidate identifiers, coordinate records and credentials. It is defense in depth, not a guarantee. Adding a workflow CI scan would happen after content is pushed and cannot undo disclosure.

Any future release of designs requires separate explicit authorization. No experiment is launched and no prediction data are uploaded by the public utilities unless explicitly invoked with local inputs.
