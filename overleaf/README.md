# Overleaf project: manuscript + supplementary information

Import this GitHub repository into Overleaf (New Project → Import from GitHub).
Everything in the repo appears in the Overleaf file tree; only this folder is LaTeX.

| File | What it is |
|---|---|
| `main.tex` | Main manuscript — Springer Nature / *Scientific Reports* template (`sn-jnl.cls`, `sn-nature` bibliography style). **Set as Main document in Overleaf.** |
| `references.bib` | Bibliography for the manuscript (`\cite{}` keys = first author + year, e.g. `larsen2021`) |
| `figures/media/` | Figures 1–3 (image2 = Fig. 1, image1 = Fig. 2, image3 = Fig. 3) |
| `supplementary/appendix_S1.tex` | Appendix S1 (Supplementary Methods, Table S1, Figures S1–S4) — a standalone document. To compile it in Overleaf: Menu → *Main document* → choose this file; or compile locally. Its PDF (`appendix_S1.pdf`) is uploaded to the journal as a separate file. |
| `supplementary/figures/media/` | Supplementary figures |
| `sn-jnl.cls`, `sn-nature.bst`, `sn-basic.bst`, `sn-template-user-manual.pdf` | Springer Nature template files and manual |

Notes
- Text was converted from Word with pandoc; yellow `\hl{}` highlights from the draft were kept (package `soul`) — remove them (and the `soul` line in the preamble) before submission.
- *Scientific Reports* asks for a single `.tex` file for the manuscript (no `\input`), SI as a separate file.
