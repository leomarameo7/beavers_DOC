# Overleaf project: manuscript + supplementary information

Import this GitHub repository into Overleaf (New Project → Import from GitHub).
Everything in the repo appears in the Overleaf file tree; only this folder is LaTeX.

| File | What it is |
|---|---|
| `main.tex` | Main manuscript for *Environmental Science & Technology* (ES&T): Introduction, Methods, Results and Discussion (no separate Conclusion). Drafted in the Springer Nature template (`sn-jnl.cls`); reformat to the ACS template at submission. **Set as Main document in Overleaf.** |
| `cover_letter.tex` | Cover letter to the ES&T editors (standalone; items marked `[TO ADD]` need author input). |
| `references.bib` | Bibliography for the manuscript (`\cite{}` keys = first author + year, e.g. `larsen2021`) |
| `figures/media/` | Main-text figures: `figure1_v2.png` (Fig. 1), `Figure_2a.png` (Fig. 2), `figure4_v1.png` (Fig. 3), `figure3_v2.png` (Fig. 4); how each is produced is listed in the top-level `README.md` |
| `supplementary/supporting_information.tex` | Supporting Information (Supplementary Methods, Tables S1–S2, Figures S1–S14) — a standalone document. To compile it in Overleaf: Menu → *Main document* → choose this file; or compile locally. Its PDF (`supporting_information.pdf`) is uploaded to the journal as a separate file. |
| `supplementary/figures/media/` | Supplementary figures |
| `sn-jnl.cls`, `sn-nature.bst`, `sn-basic.bst`, `sn-template-user-manual.pdf` | Springer Nature template files and manual |

Notes
- Text was converted from Word with pandoc; yellow `\hl{}` highlights from the draft were kept (package `soul`) — remove them (and the `soul` line in the preamble) before submission.
- *Scientific Reports* asks for a single `.tex` file for the manuscript (no `\input`), SI as a separate file.
