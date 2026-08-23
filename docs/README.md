# Documentation index

Read in this order to follow how the programme actually developed; the numbering
is chronological, not thematic.

## Where it ended up

| doc | |
|---|---|
| [19-results.md](19-results.md) | Consolidated tidal/seismology findings — the physics programme |
| [23-ml-results.md](23-ml-results.md) | **Chart features vs earthquake timing: the global result** |
| [25-stratified-results.md](25-stratified-results.md) | **The same, split by focal mechanism and depth** |

Short version: sub-monthly astronomical configurations are bounded at ~4% per
standard deviation globally and ~6–9% within any mechanism class. Outer-planet
configurations are *untestable* against earthquake catalogues, which is a weaker
statement than unsupported — see 23 §6. The sealed 2017–2024 test period was
never opened.

## Pre-registrations

Both written and committed before the analyses they govern.

| doc | |
|---|---|
| [21-preregistration.md](21-preregistration.md) | Global analysis: metric, threshold, splits, failure conditions. Two amendments, both recorded before any positive result. |
| [22-model-grid.md](22-model-grid.md) | The 26 hyperparameter configurations, fixed in advance |
| [24-stratified-preregistration.md](24-stratified-preregistration.md) | Stratified analysis, with each stratum's detectable effect stated before it was tested |

## How the work was done

| doc | |
|---|---|
| [07-research-log.md](07-research-log.md) | Chronological log. **The traps are here** — each produced a confident wrong answer, the worst reporting p = 10⁻⁸⁹ where the truth was 0.70. |
| [20-data-acquisition.md](20-data-acquisition.md) | Catalogues and external datasets, and how to refetch them |
| [04-ml-architecture.md](04-ml-architecture.md) | Design of the matched case-control machinery |

## Earlier planning

[00-framing](00-framing.md) · [01-literature](01-literature.md) ·
[05-research-frontier](05-research-frontier.md) · [08-hypotheses](08-hypotheses.md) ·
[09-deep-dive-agenda](09-deep-dive-agenda.md) · [12-build-plan](12-build-plan.md) ·
[16-plan](16-plan.md) · [17-next-steps](17-next-steps.md) ·
[18-options-roster](18-options-roster.md)

## A note on reading the nulls

Every null here is reported with the effect size it could have detected, measured
by planting a known signal into the real feature matrix and re-running the whole
pipeline. A null without that number is not a bound, and several of the results
below were *uninformative* rather than negative — the conditional logistic model
could not recover a planted 50% effect, so its null said nothing at all. The
distinction is kept explicit throughout.
