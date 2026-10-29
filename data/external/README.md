# Independent external observations

These files are primary-source supplementary data for cross-environment tests.
The publisher metadata and retrieval URLs are in `pmc-supplement-manifest.json`;
each file is bound to a SHA-256 and publisher MD5. The three studies are CC BY.
The associated article XML records the methods and attribution. Sources were
retrieved on 29 October 2026 through the public PMC cloud interface.

| Study | Data used | Interpretation |
|---|---|---|
| Buzzard et al. 2022, [doi:10.7717/peerj.12966](https://doi.org/10.7717/peerj.12966) | Supplements 7–8 and analysis scripts 10–11 | 48 Sonoran incubation bottles and 64 field collars; use corrected `Hydrogen.Uptake`, reverse sign for positive uptake |
| Tribbia et al. 2026, [doi:10.1128/aem.00265-26](https://doi.org/10.1128/aem.00265-26) | Tables S1, S7 and S8 | Ten Namib sites, 50 wet rates and ten quantified dry rates from two sites; no zero imputation at other dry sites |
| Jordaan et al. 2020, [doi:10.1128/mSystems.01131-20](https://doi.org/10.1128/mSystems.01131-20) | Table S4 | Descriptive three-family hydrogenase comparison in Australian, Namib, Gobi and Mojave soils; study-specific assays and normalizers |

Reji et al.'s forest, meadow and sand bottle series remains in
`../sources/reji-bottles.txt`, with its commit-pinned URL and checksum in the
source manifest. That study informed the activation threshold and is labelled
a development source, including when testing transfer between its soils.

The retrospective protocol is `reviews/v05/external-evaluation-protocol.md`.
Run `uv run python reproducibility/v05/reproduce.py` from the repository root.
It needs no network after dependency/reasoner installation. The optional
`fetch_pmc_cloud.py` acquisition script uses the publisher's declared cloud
metadata and verifies original checksums; it is not run during reproduction.

Measurement boundaries matter: historical Namib field water potentials are not
bottle measurements; gravimetric water content requires density for conversion
to volumetric water content; gene-copy measures and biomass-normalized bottle
rates are not interchangeable with areal ecosystem flux.
