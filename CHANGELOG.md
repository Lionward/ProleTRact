# Changelog

All notable changes to ProleTRact are documented in this file.

## [1.1.0] - 2026-02-04

### Added
- **File browser**: Browse for VCF files and folders via graphical file picker instead of typing paths. Browse button next to VCF path, population folder, and cohort folder inputs.
- **Advanced filtering**: Filter panel with criteria for motif size, copy number, chromosomes, genotypes, and pathogenic-only regions.
- **Filter presets**: Save and load filter configurations (stored in `localStorage`).
- **Mode cache**: Switching between Individual and Cohort modes preserves VCF paths, filters, and navigation state.
- **Region jump by number**: Jump directly to region N from the floating navigation bar.
- **Resizable sidebar**: Adjustable sidebar width with saved preference.
- **Pathogenicity gene search**: Search pathogenic catalog by gene name via `/api/pathogenic/search`.
- **Pathogenicity panel**: Improved layout for gene, disease, inheritance, and threshold display.
- **Documentation**: bioRxiv citation (ProleTRact/TandemTwister preprint) and updated installation instructions with `proletract --install-deps`.

---

## [1.0.1] - 2026-01-01

### Added
- **Pathogenic TR catalog**: Built-in `pathogenic_TRs.bed` with context for known disease-associated loci (gene, disease, inheritance, thresholds). Pathogenicity panel overlays this information when viewing overlapping regions.

---

## [1.0.0] - 2025-12-01

### Added
- **React frontend**: Complete rewrite from Streamlit to React and TypeScript.
- **Individual mode**: Load and analyze single VCF files.
- **Cohort mode**: Reads-based and Assembly-based VCF views for cohort analysis.
- **Region visualization**: Color-coded motifs, interruption highlighting, allele comparison.
- **Pathogenicity reference**: Overlay for known pathogenic TR loci.
- **Fast navigation**: Previous/Next controls or jump to region via coordinates.
- **PyPI package**: `pip install proletract` for easy installation.

### Changed
- **Breaking**: Streamlit-based interface replaced with React. Backend API remains FastAPI-based and compatible.

---
