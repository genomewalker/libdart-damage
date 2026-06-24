# De-circularizing deam_strata (cross-fit + null + composition gate)

## Problem
`deam_strata` bins reads by per-read `max(ga3_exc, ct3_exc)` (3' positions 1..4,
vs the read's interior) and then reports the 3' GA terminal_excess back per
stratum. Key and readout share the same terminus → the gradient is partly
tautological (selection on a noisy per-read estimator yields a monotone gradient
under H0). For ds samples whose 5' end is artifact-dead (FLB57md), the natural
cross-end witness CT5 is unavailable, so the gradient cannot be corroborated.

Panel: conductor-07f9ed8e (Opus x3) + dual-ecace1b3 (Opus + gpt-5.5, real codex).
Verified data fact: interior backgrounds CT5_bg and GA3_bg fall in lockstep
across strata (0.5531→0.4805, 0.5548→0.4815) — a GC/composition sort, a *second*
circularity that position-thinning alone does not remove.

## Design (3 layers)

### Layer 1 — cross-fit position thinning (within the 3' terminus)
Replaces split-terminus (5'-vs-3'); needs no 5' end, so FLB57md's dead 5' is moot.

Per read, deterministic fold from `fold(off) = (read_hash ^ off) & 1` where
`read_hash` = 64-bit hash of the sequence. Interleaving (not absolute even/odd)
balances the Briggs decay between folds in aggregate.

- KEY from fold-A 3' offsets (1..4): `max(ga3_exc_A, ct3_exc_A)` vs interior-A.
- READOUT raw {k,n} for GA3 (and CT3) from fold-B 3' offsets, indexed by stratum_A.
- Reciprocal: KEY fold-B → stratum_B, READOUT fold-A; **pool both directions**.

Shared-interior leakage fix: compute interior-A from fold-A interior positions
(10..14 with matching fold mask) and interior-B from fold-B. Decouples the two
estimates' null reference — this is what stops the interior-GC lockstep from
re-coupling key and readout.

Accumulators (per [length_bin][stratum], mergeable, fixed size):
  cf_k_ga3[L][S], cf_n_ga3[L][S]   // held-out GA3 successes/trials (pooled folds)
  cf_k_ct3[L][S], cf_n_ct3[L][S]   // held-out CT3 (ss axis)
  cf_reads[L][S]                    // read count
  cf_len_sum[L][S], cf_len_sumsq[L][S]
  cf_igc_sum[L][S]                  // Σ interior GC fraction (composition gate)
L = N_LEN_FINE reuse or a coarse 4-bin; S = N_OX_DEAM_STRATA (5).

### Layer 2 — read-level resampled null
Analytic binomial z is anticonservative (per-read overdispersion + cross-fit reuse).
Single-pass faithful null = bounded reservoir of per-read tuples
  {stratum, k_ga3_heldout, n_ga3_heldout, len, igc}  (reservoir ~200k).
Finalize: (a) key-label permutation within (length × igc) cells → null gradient
distribution → effect-size floor; (b) read-level bootstrap → CI on the gradient.
Reservoir reuses the existing LSD reservoir pattern (lsd_cnt) in WorkerState.

### Layer 3 — composition residualization gate
Emit `interior_gc[S]`. Gradient must survive regressing out the interior-GC trend:
`heldout_te_resid` = held-out GA3 gradient after removing the component linear in
interior_gc[S]. Null bins in Layer 2 are composition-matched on igc.

## Pass/fail (emit strata only if all hold)
1. pooled held-out cross-fit GA3 gradient monotone, same sign s0→s4;
2. exceeds read-level resampled null (>99th pct) AND bootstrap CI lower > 0;
3. survives interior-composition adjustment (heldout_te_resid keeps sign+significance).
Axis: GA3 for ds, CT3 for ss. Legacy full-data gradient demoted to
`circular_diagnostic` (descriptive only, never an authenticity claim).

## Threshold calibration (not magic constants)
Power-set from reservoir count n: permutation floor at 99th pct (data-driven, no
fixed rho); CI from bootstrap. Report `decirc_n_reservoir` so thresholds are auditable.

## Code map
- struct: libtaph/include/taph/sample_damage_profile.hpp — add cf_* arrays + reservoir hooks.
- accumulate: libtaph/src/damage_estimation_update.cpp ~190 (keying) / ~390 (readout).
- merge: same file ~1373 (add cf_* to merge loop).
- reset: same file ~2080.
- finalize: damage_estimation_finalize*.cpp — null + residual + gates.
- emit: libtaph/src/profile_json.cpp deam_strata block — held-out stats, gates,
  circular_diagnostic, interior_gc, decirc_n_reservoir.

## Validation checkpoint (BEFORE 109-sample rerun)
Build, run FLB57md merged. Require: held-out GA3 gradient still positive+monotone
and clears gate 3 (GA3 move +0.30 >> bg move −0.07). If held-out gradient
collapses to ~null, the original strata WERE largely a sort artifact — report that
honestly rather than shipping.

## Open (from panel, track but don't block)
- length-as-validator: usable only variance-matched on interior base count.
- overdispersion magnitude: confirm read-level resampling captures it.
- CpG reference-free: deferred (dinucleotide-context accumulators feasible? later).
