# DEFECT — channel F/G/H: pooled z and MH strata use different position windows

Status: **filed, not fixed.** Found while attributing the `count_table_golden` red on
`ds_artifact_cap_masked` (FLB57, a ds blunt-end fill-in library).

## Problem

When the pos-0 fill-in artifact flag fires, the terminal window is *supposed* to start at p=1,
because a ds blunt-end fill-in prepends a non-template G at position 0. Only half the emitted
count table honours that.

- `damage_estimation_finalize_oxidation.cpp:286` sets `p0_tc_5 = position_0_artifact_5prime ? 1 : 0`.
  The primary aggregates (`:301`) and the Mantel-Haenszel strata (`:329`, `:413`) both loop
  `p in [p0_tc_5, 5)` — they drop p=0 as intended.
- `damage_estimation_update.cpp:985` accumulates the prefix-binned terminal arrays
  (`ca_pre_terminal_by_pfx`, `ca_stop_terminal_by_pfx`, `ca_deam_shadow_terminal_by_pfx`) under a
  **hardcoded `p < 5`** at *stream* time — before the flag exists, since the flag is a composition
  statistic computed in finalize.
- `damage_estimation_channels.cpp:57-70` (`recompute_fgh_excluding_adapter_prefixes`) rebuilds the
  emitted aggregates *and* `profile.channel_f_z` from exactly those prefix-binned arrays, and
  preserves the strata from the primary row (`:49`, and the comment at `:43-45`).

So on any library where the pos-0 flag fires **and** adapter-prefix exclusion is active, the
emitted `channel_f_z` is computed over p in [0,5) while `channel_f_mh_z` is computed over p in
[1,5). Two statistics whose entire job is to be compared against each other — pooled z vs
stratified MH, the Simpson-reversal check — are computed on different position windows.

The two preconditions intersect precisely in the `ds_artifact_cap_masked` fixture, which is the
fixture that exists to police pooled-z-vs-MH reversal.

## Evidence

Job 20134269: one binary, one 4.1G input (FLB57), the pos-0 threshold the only variable.

    detector ON  (thr 0.07):  g0=0.4110 gb=0.2690 diff=0.1420 flag=1
    detector OFF (thr 9.0):   g0=0.4110 gb=0.2690 diff=0.1420 flag=0

Channel F, terminal strata `a_term_stop` / `b_term_other`:

    flag=0   a= 3023111  2145707  5896122   b= 15122309  14189942  13331369
    flag=1   a= 2585427  1721629  4909845   b= 13033701  12604793  10297384

Aggregates, **byte-identical under both**: `term_pre` 25187155, `term_stop` 11060731,
`z_den_term` 53688263, `term_shadow` 17440377, `pre_clamp_z` -39.763765176654886.

Strata move, aggregates do not. That is the signature of the split above, and it is only possible
because the aggregates are being sourced from the prefix-binned (flag-blind) path.

## Fix sketch

Make the prefix-binned terminal arrays carry the position axis, so the p=0 slice can be dropped in
finalize once the flag is known. Cheapest form: split `ca_*_terminal_by_pfx` into a p=0 bin and a
p in [1,5) bin (2 x 4096 doubles x 3 arrays), and have `recompute_fgh_excluding_adapter_prefixes`
resum from `p0_tc_5` onward, matching `finalize_oxidation.cpp:301`. The p<5 constant at
`update.cpp:985` then becomes the only window definition in the code, instead of one of two.

## Why it is not fixed in the same commit as the re-capture

Fixing it **changes `channel_f_z`** on every fill-in library — a real estimator change with its own
validation burden (it must be shown not to disturb the four passing count-table fixtures, and the
verdict routes downstream of channel F/G/H must be re-measured). Bundling that into a commit whose
purpose is to make the suite green would hide an estimator change inside a test-hygiene change.

The re-captured `ds_artifact_cap_masked` fixture therefore records **current behaviour, asymmetry
included**. Fixing this defect will legitimately and loudly re-capture it again.

## Adjacent, trivial, also unfixed

`damage_estimation_finalize_init.cpp:76` says "threshold 0.10 separates them with margin"; `:85`
tests `> 0.07`. FLB57 sits at +0.142 and clears both, so the attribution above does not depend on
which number is right — but the comment and the code must be reconciled.
