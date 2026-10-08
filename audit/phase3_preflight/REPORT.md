# PlaViSca final preprocessing/data acceptance closure

> **Superseded (2026-10-08):** the current release is
> `plavisca-v1-gametocytes-2026-10-08` (gametocyte stage precedence, rRNA
> removal, refreshed integration features); see
> `audit/gametocyte_rebuild/REPORT.md` and
> `acceptance_matrix_AT01_AT31.tsv`. The text below describes the
> 2026-09-21 release.

> **Authoritative DEC01/DEC02 update (2026-09-21): FINAL RELEASE DATA CANDIDATE — DEC01/DEC02 RESOLVED.**
> This section supersedes the historical release-closure narrative below.
> DEC01 and DEC02 are implemented, AT14 passes, and focused scientific
> resolution of AT25 establishes 30 PASS / 0 FAIL / 0 PENDING.

## DEC01/DEC02 final rebuild verdict

The six study objects and every membership-dependent downstream analysis were
rebuilt from the validated count artifacts without FASTQ/alignment/count
reprocessing. Mancio-Silva2022 now contains 1,480 cells; the other study counts
remain 9,766 (Sa2020), 3,294 (Hazzard2022), 80,024 (Hazzard2024), 1,438
(Ruberto2022_1), and 9,947 (Ruberto2022_2), for 105,949 total cells. The common
integration universe remains 5,203 features.

DEC01 removes exactly the 14 proven `.1` representations listed in
`audit/mancio_silva2022/dec01_removed_duplicate_cells.tsv`; every unsuffixed
source counterpart remains. DEC02 retains `D5Seq1_CCCCGATTGACG` and
`D5Seq2_CCCCGATTGACG`, whose distinct source samples (`SEQ5r1`, `SEQ5r2`) and
run sets remain recorded, without a D5 duplicate/suspicion flag. AT14 passes.

Annotation validation passes: 90,219 cells are IDC-eligible/labeled, 22,795
cells have gametocyte calls, no non-blood IDC label or impossible lifecycle
state exists, and the only seven paired missing-stage rows are the documented
low-UMI Hazzard2024 cells. Integration mixing was evaluated within comparable
broad-stage populations and improved in all three multi-study strata.

AT25's original failure is classified **VALIDATION_SEMANTICS_DEFECT**. The old
test averaged all atlas cells per cluster, whereas production calculates marker
means only among `gametocyte_eligible <- idc_eligible` cells. It also required
every called all-cell cluster to clear thresholds while separately allowing up
to 10% boundary mixing, making that allowance ineffective. The corrected test
keeps all thresholds unchanged and separates production consistency from
integrated-population coherence.

Female: n=21,801, AUC=0.8647923, median shift=0.4847102, eligible-cluster
score/margin minima=0.3735319/0.3324021, all-cell coherence=1.0 — PASS. Male:
n=994, AUC=0.9748556, median shift=1.8501185, eligible-cluster score/margin
minima=0.1284835/0.1578435, all-cell coherence=976/994=0.9818913 — PASS. The
18/994 non-clearing calls remain visible in cluster 4, which contains 1,687
total cells, 18 eligible Male calls, and 1,669 ineligible cells. No threshold
was lowered, and no atlas, Harmony, or SingleR rebuild occurred. The final
acceptance result is **30 PASS / 0 FAIL / 0 PENDING_DECISION / 0 PENDING_APP**.
All six frozen scientific/application artifacts retain their DEC01/DEC02
SHA256 values; exact comparisons are recorded in
`AT25_scientific_artifact_hash_verification.tsv`.

Scientific rebuild commit: `173431cfca8e19c716a96164ee2ea6053441855b`.
Authoritative manifest: `audit/phase3_preflight/final_release_manifest.tsv`.

## Historical release-closure report (superseded)

## Verdict

The preprocessing/data candidate is now **FINAL RELEASE DATA CANDIDATE**. All 29 technically actionable AT01–AT30 criteria pass. AT14 alone remains `PENDING_DECISION` because DEC01/DEC02 still lack explicit project-lead approval; no cell membership was changed to force closure. There are zero FAIL and zero PENDING_APP results.

No STARsolo, kallisto, FASTQ processing, Harmony, PCA, UMAP, t-SNE, SingleR, or gametocyte classification was rerun. The saved Phase-2 atlas remained the scientific input and its SHA256 is unchanged.

## Identity and scope

- Phase-2 scientific commit: `54c6b9581ff5268422ddffda42ce457bf09c5b55`
- Independent QA base: `a57be450a84d42e47a67d1120320e8a5536797dd`
- Frozen closure code snapshot: `8ce8828985aa75ea625bb4e762465f6ab8689905`
- Closure branch: `plavisca-release-closure`
- App commit inspected read-only: `c44cde046acf9791f5dd8a83f884d2f83ec08486`
- Authoritative manifest: `audit/phase3_preflight/final_release_manifest.tsv`

## Closed failures

### AT02 — PASS

`audit/mancio_silva2022/sa2020_authoritative_cell_lineage.tsv` contains exactly 9,766 retained Sa2020 cells, 9,766 unique source coordinates, 9,766 unique study IDs, and 9,766 unique final IDs, with zero missing or extra final cells. The source coordinate is the exact audited STARsolo `run_id:raw_barcode`. Author-workbook IDs are retained for observed rows; the 334 retained cells absent from the lower-threshold workbook are represented by their exact run/barcode source identity without inventing author IDs.

### AT17 — PASS

Canonical order is lexicographic sort of the final cell key. Raw, normalized, scaled, and cleaned metadata each have 105,963 unique keys with identical sets and identical order. Production `flatten_data.R` now applies and asserts this order before serialization.

Only `cleaned_dataset.rds` needed export regeneration. Its old SHA256 was `30b61ccdad19ef52eb835bf4e644800b2151497d14066438bd926a717d3f37a0`; final SHA256 is `0be6e6dcde0721d11193d89f1a122ae87b24aac5e8ffc776fcd9262ddbfba7a8`. The script preserved `pca_df`, `hiv_data`, and `top_genes_exp` identically and changed only `mr_data` ordering. Raw, normalized, scaled, study objects, and atlas hashes are unchanged.

### AT23 — PASS

`final_release_manifest.tsv` is the unambiguous authoritative selector. The historical root `build_manifest.tsv` remains unchanged. The final manifest records code/QA/scientific lineage, release branch/date, R and package versions, pixi manifest/lock checksums, seeds/parameters, six study counts, Ruberto count provenance, key scientific/reference inputs, and exact path/size/SHA256/count for every final output. `validate_final_manifest.R` independently recomputed all listed hashes and sizes.

### AT25 — PASS

The test recreates the unchanged production AddModuleScore definitions on an in-memory copy of the saved normalized atlas and relates them to integrated clusters; it does not rerun Harmony or modify the atlas.

Predeclared criteria per sex were: own-score discrimination from blood-stage uncalled cells with AUC >= 0.70; positive median score shift; every called integrated cluster meeting the production mean-score threshold >= 0.10 and own-minus-opposite margin >= 0.05; and >= 90% called-cell integrated-cluster coherence.

- Female: AUC 0.8584; median shift 0.4663; minimum called-cluster score 0.3512; minimum margin 0.3077; coherence 1.000 — PASS.
- Male: AUC 0.9772; median shift 1.8624; minimum called-cluster score 1.5281; minimum margin 1.5627; coherence 1.000 — PASS.

## Other closure findings

AT19 and AT20 independently pass at immutable app commit `c44cde0`. Its committed evidence reports all 40 configured fields resolving, zero live obsolete-schema violations, and all dedicated AT19/AT20 tests passing. Because two data files changed afterward, the app team must rerun the same suite against the final candidate hashes.

`data_source.csv` now has one verified `study_pmid` per study and exact agreement with atlas metadata. Old SHA256 was `cb1dca3ff9247634eda474dde693ce16382343cf5bd29e70a59c765c27bdaf98`; final SHA256 is `cd17242801bb863251c7782cfd240f55af26f2d54fa4ebe7366b03cb6a8e892c`.

The 21 paired stage-missing cells are not a new regression: seven Hazzard2024 cells have 5–9 UMI and are below the documented 10-UMI inference gate; 14 Mancio-Silva2022 cells are exactly the retained/flagged `.1` duplicate representations lacking source/inferred stage. Both stage fields are missing together. Schema nullability was corrected to document these cases; no stage was fabricated.

## Final hashes

- Atlas: `15fe6371810f53e3dcb0210a2fe70a19da1c7c607ec3adb8cec0c64a6a519c8b`
- Normalized: `8a87841094655148422d157635e6ed50d4bb9dba7735798cc9e415ae069a8740`
- Raw: `15df4435d666cef1e3ec4d2d977f507a4415ecf3a8bd5d952199dda42e009cf8`
- Scaled: `86d08ff3eea411c01437bc0f6ca461bd65184ebba032cab381adb700e7ea1cbf`
- Cleaned: `0be6e6dcde0721d11193d89f1a122ae87b24aac5e8ffc776fcd9262ddbfba7a8`
- Bibliography: `cd17242801bb863251c7782cfd240f55af26f2d54fa4ebe7366b03cb6a8e892c`

## Application retest handoff

The app team must replace/retest these changed files:

1. `cleaned_dataset.rds` — changed row ordering only.
2. `data_source.csv` — added `study_pmid`.

The app contract test should still load raw and normalized exports alongside cleaned metadata to verify identical order; their bytes are unchanged. Scaled export bytes are also unchanged. No deployment was performed.
