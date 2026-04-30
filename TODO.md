Let me check a few specific files before giving you an accurate list.

Good news — `HI_inference.py` is fully implemented, it's just not wired into `pipeline.py`. Here's a clean todo list:

---

**Core Functionality**
- [x] Wire `HI_inference()` into `pipeline.py` step 9 (replace `pass`)
- [x] Add `HI_inference` section to `example_config.yaml`
- [x] Implement `config.py` — YAML schema validation with clear error messages for missing/invalid fields (marked `# TODO:` in pipeline.py)

**Testing**
- [ ] Write unit tests for `count_reads`, `normalize_mosdepth`, `find_neighbors`, `compute_dipcn`, `HI_inference`
- [ ] Write an end-to-end integration test with a small example CRAM dataset
- [ ] Replace `test/test_add.py` placeholder

**Documentation**
- [ ] Document all config YAML fields (types, defaults, required vs optional) — `example_config.yaml` is a start but needs prose explanations
- [ ] Document output file formats (columns, units, what each means biologically)
- [ ] Add a quickstart / worked example in the README with real or synthetic data
- [ ] Document `hardcoded_positions.txt` format and how it was generated
- [ ] Add troubleshooting section (common errors, HPC/SLURM tips)
- [ ] Publish docs to ReadTheDocs (links in README are dead)

**Code Cleanup**
- [x] Update `examples/SOL.py` — references old commented-out CLI commands, will fail if run
- [x] Remove or formally deprecate `google_cloud_copy.py`
- [x] Uncomment or clean up commented-out log statements in `pipeline.py`

**Infrastructure**
- [x] Set up CI/CD (GitHub Actions) — pyproject.toml references it but `.github/workflows/` is missing/empty
- [x] Add version bump process before tagging v1.0 release

---

## JOSS Paper

**Paper (`paper.md` / `paper.bib`)**
- [ ] Verify `mukamel2021` citation — confirm title is *Protein-Coding Repeat Polymorphisms Strongly Shape Diverse Human Phenotypes*, Science 2021
- [ ] Fill in `hujoel2026` — add actual title and replace `doi: VERIFY`
- [ ] Fill in `didericksen2024` or replace with the actual long-read citation you intend to use (or remove it)
- [ ] Fill in grant/funding info in Acknowledgements
- [ ] Fix typo in old title (`Utalizes`) — already updated in the new title but double-check rendered output
- [ ] Confirm `kamstrup2010` — the bib entry year says 2009 (JAMA), key says 2010; pick one and be consistent

**Software requirements (JOSS editorial checklist)**
- [x] Confirm the software is installable from scratch — run `pip install -e .` in a clean conda env and verify
- [x] Add a `CHANGELOG.md` or `CHANGES.md` — JOSS reviewers often look for this
- [x] Ensure `CITATION.cff` is complete and valid — it exists but verify all fields match the paper authors/version
- [x] Tag a release on GitHub matching the version in `__init__.py` (currently `0.1.0`) — JOSS requires a citable release with an archive (e.g., Zenodo DOI)
- [ ] Get a Zenodo DOI for the release — link it in the README and `CITATION.cff`

**Documentation (JOSS specifically checks this)**
- [ ] README must include: installation, a minimal working example, and link to full docs — currently missing a runnable example
- [ ] API/usage documentation must cover all config YAML fields — partial right now
- [ ] Publish docs to ReadTheDocs (README links are dead)

**Tests (JOSS requires evidence the software works)**
- [ ] At minimum, add one integration test that runs the pipeline on a small synthetic dataset
- [ ] JOSS won't require full coverage but reviewers will flag zero real tests