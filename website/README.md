# VTA-Agent — showcase website

A single self-contained `index.html` (no build step, no dependencies) for the hackathon
demo. It replays a **real committed benchmark run** (SARS-CoV-2 Mpro) and tells the honest
"agent that catches its own inflation" story. All numbers are from committed artifacts.

## Run locally
```bash
# any static server works; e.g.
python3 -m http.server -d website 8080
# then open http://localhost:8080
```
Or just double-click `website/index.html`.

## Deploy (pick one — all free)
- **GitHub Pages:** push the repo, Settings → Pages → deploy from `/website` (or move
  `index.html` to repo root / `docs/`). Live at `https://<user>.github.io/<repo>/`.
- **Netlify:** drag the `website/` folder onto app.netlify.com/drop → instant URL.
- **Vercel:** `vercel deploy website` (framework preset: Other).
- **Cloudflare Pages:** connect repo, build command none, output dir `website`.

## What's on the page
1. Hero + tagline
2. **Interactive pipeline demo** — "Run the agent on SARS-CoV-2 Mpro" animates classify →
   structure → pockets → dock → rank → validate with the real numbers.
3. **The validation ladder** — the 3/4 → 1/4 → CIs → powered → beaten-by-2D-sim → redock arc.
4. **Headline finding** — animated BEDROC bars (random / 2D-sim / Vina).
5. Three-target results, the nucleotide meta-finding, how-it-works, honest-use footer.

## Editing
- Everything is in `index.html` (inline CSS + vanilla JS). No framework.
- Numbers to update if the benchmark changes live near the `#results` section and in the
  `steps` array of the `<script>` (the demo replay).
- Slide 11 of `docs/VTA_AGENT_DECK.md` tracks the two extended docks (Mpro 31:1, HCV-NI
  property-unmatched); drop those results into a new card here once they finish.

## Honesty note
The page states plainly that outputs are ranked, uncertainty-bearing **hypotheses**, not
clinical claims — keep that framing in any edits.
