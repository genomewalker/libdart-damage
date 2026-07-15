# Handoff spec — unify the single-end damage profiling loop into `libtaph::profile_reads()`

**Audience:** the engineer developing **libtaph**. You have **libtaph**, **fqdup**, and **dart** checked out.
**Goal:** one code path — `taph::profile_reads()` — produces a finalized `SampleDamageProfile` from a read
source. `fqdup profile` (single-end / merged) and `dart predict --allow-self-profile` both call it, so a
byte-identical `SampleDamageProfile` is guaranteed for the same reads + same config. The two engines have
been drifting (dart self-profile lands `LOW_ABUNDANCE` where fqdup `DETECT`s on identical reads —
see `dart/src/main.cpp` comment at the `library_summary` block); collapsing to one loop removes the drift
surface.

---

## 0. Pre-flight checks — VERIFIED against the tree 2026-07-15

1. **`clip_adapters` current body — CONFIRMED, line ref corrected.** It lives at
   `libtaph/src/damage_profiler.cpp:52` (not `:64-70`). The real body matches Section 3's reconstructed
   "before" semantically: single 5' prefix strip, **single** 3' strip. Two mechanical differences from the
   §3 "before" snippet, harmless to the diff: the real code uses `s.substr(...)` (not `s.erase(...)`) and has
   **no `L >= 12` guard** on the 3' branch. The §3 "after" *adds* the `L < 12` break, so the iterative
   promotion is not a pure superset for 6–11 bp remnants — call that out in the commit message.
2. **`run_prescan` parity — CONFIRMED.** `taph::run_prescan` (`src/damage_profiler.cpp:8`) builds the 3'
   terminal hexamer histogram internally (`:31` `++hex3[code]; ++n_hex3;`) and passes it to
   `detect_adapter_stubs` (`:38`). It replaces fqdup's hand-rolled parallel pre-scan
   (`fqdup/src/profile.cpp:266-332`) with identical stub inputs.
3. **`finalize_sample_profile(dp, n_threads)` thread-invariance — NOT verified.** Run **both engines with
   `-p 1`** for the byte-identical verification (Section 6) so this is moot. Do not assume invariance.
4. **Field names — CORRECTED (Section 6 harness would not compile as written).** On `SampleDamageProfile`
   (`include/taph/sample_damage_profile.hpp`): `pi` is a `DamageEstimate` (`damage_estimate.hpp:97`, has
   `point/lo/hi/state`), and `d_max_5prime` / `d_max_3prime` / `d_max_combined` / `n_reads` exist. But:
   - `d_max_se` is **not** a profile member — it lives on `BulkDamageResult` (`bulk_damage_model.hpp:185`),
     reachable as **`p.bulk_damage.d_max_se`** (`bulk_damage` is a `SampleDamageProfile` member, `:1598`).
   - `n_elig_total` **does not exist** on the profile — the only `n_elig` is the per-position
     `PiPosCount::n_elig`; the aggregate `n_elig_total` is a *local* in `read_ancient_llr.cpp`, not a
     profile field. Section 6 Harness A has been updated to drop it and use `p.n_reads` as the abundance
     proxy. If the read-LLR eligibility total is genuinely needed for the drift comparison, source it from
     the `read_ancient_llr` path, not from `SampleDamageProfile`.

---

## 1. NEW — additions to `libtaph/damage_profiler.hpp`

Append below the existing `clip_adapters` declaration, inside `namespace taph`.

```cpp
// ── Whole-source profiling ─────────────────────────────────────────────────────
//
// profile_reads() runs the complete single-end / merged damage-profiling loop:
//   1. adapter pre-scan (delegates to run_prescan) → AdapterStubs
//   2. multi-threaded full pass: per-read adapter clip → update_sample_profile,
//      with an optional caller hook per full-length read for side channels
//   3. merge per-thread profiles → finalize_sample_profile
// The returned profile is FINALIZED and ready to read.
//
// One engine, one result: fqdup `profile` and dart `predict` self-profiling both
// call this, so the finalized SampleDamageProfile is byte-identical for identical
// reads + identical ProfileConfig (see verification plan).

// Factory returning a BatchFn that reads sequences from the START of the source.
// profile_reads() invokes it up to twice: once for the pre-scan (unless stubs are
// supplied in cfg), once for the full pass. Each returned BatchFn MUST start at the
// head of the source and append up to `max` sequences to `out`, returning true, or
// return false (leaving `out` empty) when exhausted.
using ReaderFactory = std::function<BatchFn()>;

// Called once per FULL-LENGTH read (clipped length L >= cfg.min_read_len) AFTER
// adapter clipping, from worker thread `tid` (0 <= tid < cfg.threads). Use it to
// accumulate side channels (OxoG sig histograms, length-stratified counts, etc.)
// into caller-owned, per-thread state indexed by `tid`. Do NOT touch a
// SampleDamageProfile here — profile_reads owns that. Thread-safe requirement:
// only ever touch the slot for `tid` (no cross-tid sharing → no locks needed).
using PerReadHook = std::function<void(int tid, const std::string& seq, int L)>;

struct ProfileConfig {
    int    threads              = 1;
    SampleDamageProfile::LibraryType forced_lib
        = SampleDamageProfile::LibraryType::UNKNOWN;
    int    min_read_len         = 30;         // full-length threshold (fqdup LSD_L_MIN)
    int    short_fragment_floor = 30;         // set on EVERY per-thread profile before
                                              // the first update (load-bearing, see .cpp)
    int    max_read_len         = 0;          // 0 = off; drop reads longer than this
    size_t prescan_reads        = 1'000'000;  // 0 = scan all; qualifying-read cap for
                                              // the internal run_prescan
    // Optional: skip the internal pre-scan and use these stubs verbatim (dart supplies
    // the stubs it already computed for its Pass-2 per-read clip). nullptr = pre-scan.
    const AdapterStubs* precomputed_stubs = nullptr;
};

struct ProfileResult {
    SampleDamageProfile profile;      // merged + FINALIZED
    AdapterStubs        stubs;        // detected (or the supplied precomputed_stubs)
    int64_t reads_scanned = 0;        // clipped L >= min_read_len  (full-length)
    int64_t reads_skipped = 0;        // clipped L <  min_read_len  (short/empty)
    int     len_min = INT_MAX;        // over full-length reads only
    int     len_max = 0;
    int64_t len_sum = 0;
    int64_t n_stub5_hits        = 0;  // UNCLIPPED 5' hexamer == a stub5 code
    int64_t n_stub3_hits        = 0;  // UNCLIPPED 3' hexamer == a stub3 code
    int64_t n_stub_reads_checked = 0; // reads where the stub tally ran
};

ProfileResult profile_reads(const ReaderFactory& make_reader,
                            const ProfileConfig& cfg,
                            const PerReadHook&   hook = {});
```

`<functional>`, `<climits>` (for `INT_MAX`), `<cstdint>` are already pulled in by the existing header.

---

## 2. NEW — `libtaph/src/damage_profiler.cpp` : `profile_reads()` definition

Add to the existing translation unit (the one that already defines `run_prescan`, `encode_stub_codes`,
`clip_adapters`). It reuses `run_prescan` for the pre-scan and mirrors fqdup's `worker_fn` /
`clip_worker_fn` read-handling exactly, so integer counts merge bit-exactly regardless of how reads
partition across threads.

```cpp
#include "taph/damage_profiler.hpp"
#include "taph/frame_selector_decl.hpp"     // FrameSelector, encode_hex_at
#include "taph/sample_damage_profile.hpp"

#include <algorithm>
#include <climits>
#include <condition_variable>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

namespace taph {

namespace {

// Bounded batch queue — 1:1 with fqdup WorkQueue (damage_worker.hpp).
struct SeqQueue {
    std::mutex              mtx;
    std::condition_variable cv_ne, cv_nf;
    std::vector<std::vector<std::string>> batches;
    bool done = false;
    int  max_depth;
    explicit SeqQueue(int d) : max_depth(d) {}
    void push(std::vector<std::string>&& b) {
        std::unique_lock<std::mutex> lk(mtx);
        cv_nf.wait(lk, [&]{ return (int)batches.size() < max_depth; });
        batches.push_back(std::move(b));
        cv_ne.notify_one();
    }
    bool pop(std::vector<std::string>& b) {
        std::unique_lock<std::mutex> lk(mtx);
        cv_ne.wait(lk, [&]{ return !batches.empty() || done; });
        if (batches.empty()) return false;
        b = std::move(batches.back());
        batches.pop_back();
        cv_nf.notify_one();
        return true;
    }
    void set_done() {
        std::unique_lock<std::mutex> lk(mtx);
        done = true;
        cv_ne.notify_all();
    }
};

struct ThreadAccum {
    int64_t reads_scanned = 0, reads_skipped = 0;
    int     len_min = INT_MAX, len_max = 0;
    int64_t len_sum = 0;
    int64_t h5 = 0, h3 = 0, hc = 0;   // unclipped stub-hit tallies
};

constexpr int BATCH_SZ = 8192;   // matches fqdup profile.cpp BATCH_SZ

} // namespace

ProfileResult profile_reads(const ReaderFactory& make_reader,
                            const ProfileConfig& cfg,
                            const PerReadHook&   hook) {
    const int nt = std::max(1, cfg.threads);
    ProfileResult result;

    // ---- pre-scan (or accept caller-supplied stubs) -----------------------
    if (cfg.precomputed_stubs) {
        result.stubs = *cfg.precomputed_stubs;
    } else {
        BatchFn scan_fn = make_reader();                       // fresh reader @ head
        PrescanResult pr = run_prescan(scan_fn, cfg.forced_lib,
                                       cfg.min_read_len, cfg.prescan_reads);
        result.stubs = std::move(pr.stubs);
    }
    const bool do_clip = result.stubs.adapter_clipped || result.stubs.adapter3_clipped;
    const bool have_stubs = !result.stubs.stubs5.empty() || !result.stubs.stubs3.empty();
    const std::vector<uint32_t> stub5 = encode_stub_codes(result.stubs.stubs5);
    const std::vector<uint32_t> stub3 = encode_stub_codes(result.stubs.stubs3);

    // ---- per-thread profiles ---------------------------------------------
    // LOAD-BEARING: set short_fragment_floor (and forced_library_type) on EACH
    // per-thread profile BEFORE its first update_sample_profile. update reads the
    // floor to decide, on the very first read, whether to run the terminal-only
    // short path or abstain. Setting it after the first update diverges silently.
    std::vector<std::unique_ptr<SampleDamageProfile>> tp(nt);
    std::vector<ThreadAccum> acc(nt);
    for (int t = 0; t < nt; ++t) {
        tp[t] = std::make_unique<SampleDamageProfile>();
        tp[t]->forced_library_type  = cfg.forced_lib;
        tp[t]->short_fragment_floor = cfg.short_fragment_floor;
    }

    // ---- full pass --------------------------------------------------------
    SeqQueue queue(2 * nt);
    std::vector<std::thread> workers;
    workers.reserve(nt);
    for (int t = 0; t < nt; ++t) {
        workers.emplace_back([&, t]() {
            SampleDamageProfile& prof = *tp[t];
            ThreadAccum&         a    = acc[t];
            std::vector<std::string> batch;
            while (queue.pop(batch)) {
                for (std::string& raw : batch) {
                    const int L0 = (int)raw.size();
                    if (cfg.max_read_len > 0 && L0 > cfg.max_read_len) continue;

                    // Unclipped stub-hit tally (drives JSON read fractions). Runs
                    // whenever stubs exist, independent of do_clip. Code-compare is
                    // equivalent to fqdup's seq.compare(0,6,stub) — an N in the 6-mer
                    // gives encode_hex_at == -1 (no match), same as a string mismatch.
                    if (have_stubs && L0 >= 6) {
                        const int c5 = encode_hex_at(raw, 0);
                        const int c3 = encode_hex_at(raw, L0 - 6);
                        bool m5 = false, m3 = false;
                        if (c5 >= 0) for (uint32_t c : stub5) if ((uint32_t)c5 == c) { m5 = true; break; }
                        if (c3 >= 0) for (uint32_t c : stub3) if ((uint32_t)c3 == c) { m3 = true; break; }
                        a.h5 += m5; a.h3 += m3; ++a.hc;
                    }

                    // Clip via the unified clip_adapters. min_len = 0 → "pure clip"
                    // mode: return the clipped bytes even below the floor, so the
                    // short-fragment terminal path still sees them (see Section 3).
                    std::string seq = do_clip ? clip_adapters(raw, stub5, stub3, 0)
                                              : std::move(raw);
                    const int L = (int)seq.size();

                    if (L < cfg.min_read_len) {
                        // Short/empty after clip. Route 16..(min-1) bp reads to the
                        // terminal-only short path ONLY when the floor is lowered;
                        // update_sample_profile enforces the floor + abstain itself.
                        if (prof.short_fragment_floor < cfg.min_read_len)
                            FrameSelector::update_sample_profile(prof, seq);
                        ++a.reads_skipped;
                        continue;
                    }

                    FrameSelector::update_sample_profile(prof, seq);
                    if (L < a.len_min) a.len_min = L;
                    if (L > a.len_max) a.len_max = L;
                    a.len_sum += L;
                    ++a.reads_scanned;
                    if (hook) hook(t, seq, L);     // side channels: OxoG / LSD / etc.
                }
            }
        });
    }

    // ---- producer: one forward pass, batched -----------------------------
    try {
        BatchFn feed = make_reader();              // fresh reader @ head
        std::vector<std::string> src, batch;
        batch.reserve(BATCH_SZ);
        while (feed(src, BATCH_SZ)) {
            for (auto& s : src) {
                batch.push_back(std::move(s));
                if ((int)batch.size() == BATCH_SZ) {
                    queue.push(std::move(batch));
                    batch.clear();
                    batch.reserve(BATCH_SZ);
                }
            }
            src.clear();
        }
        if (!batch.empty()) queue.push(std::move(batch));
    } catch (...) {
        queue.set_done();
        for (auto& w : workers) w.join();
        throw;
    }
    queue.set_done();
    for (auto& w : workers) w.join();

    // ---- merge (thread-index order) + finalize ---------------------------
    result.profile.forced_library_type  = cfg.forced_lib;
    result.profile.short_fragment_floor = cfg.short_fragment_floor;
    for (int t = 0; t < nt; ++t) {
        FrameSelector::merge_sample_profiles(result.profile, *tp[t]);
        const ThreadAccum& a = acc[t];
        result.reads_scanned        += a.reads_scanned;
        result.reads_skipped        += a.reads_skipped;
        result.len_sum              += a.len_sum;
        result.n_stub5_hits         += a.h5;
        result.n_stub3_hits         += a.h3;
        result.n_stub_reads_checked += a.hc;
        if (a.len_min < result.len_min) result.len_min = a.len_min;
        if (a.len_max > result.len_max) result.len_max = a.len_max;
    }
    result.profile.forced_library_type = cfg.forced_lib;   // survive the merge

    FrameSelector::finalize_sample_profile(result.profile, nt);
    return result;
}

} // namespace taph
```

Notes:
- Merge order is fixed (thread index 0..nt-1); integer count fields add associatively so the merged profile
  is identical no matter how the queue happened to hand batches to threads. This is the same guarantee
  fqdup relies on (`fqdup/src/profile.cpp:294-305` comment: "merge associatively and bit-exactly regardless
  of read partition").
- `reads_scanned` = full-length; `reads_skipped` = short/empty. This matches fqdup's `worker_fn`
  (`damage_worker.hpp` `worker_fn` / `clip_worker_fn`): short reads bump `reads_skipped` even when the
  short-fragment path profiles them.

---

## 3. Clip unification — promote fqdup's iterative 3' stripping into `clip_adapters`

fqdup's clip worker strips **multiple** stacked 3' stubs in a loop
(`fqdup/damage_worker.hpp`, `clip_worker_fn`, the `do { … } while (trimmed)` block — your line refs
`314-328`). libtaph's `clip_adapters` currently strips a **single** 3' stub. Promote the loop into
`clip_adapters` so every caller (fqdup full pass, dart Pass-2 per-read clip, `profile_reads`) shares one
behaviour. Also add the `min_len <= 0` "pure clip" mode that `profile_reads` uses.

### Before — `libtaph/src/damage_profiler.cpp:64-70` (reconstructed from the header contract + call sites)

```cpp
std::string clip_adapters(const std::string&           seq,
                          const std::vector<uint32_t>& stub5_codes,
                          const std::vector<uint32_t>& stub3_codes,
                          int                          min_len) {
    std::string s = seq;
    int L = (int)s.size();
    // 5': single hex-code prefix match
    if (L >= 6 && !stub5_codes.empty()) {
        const int c = encode_hex_at(s, 0);
        if (c >= 0)
            for (uint32_t code : stub5_codes)
                if ((uint32_t)c == code) { s.erase(0, 6); break; }
    }
    // 3': SINGLE terminal hex-code match
    L = (int)s.size();
    if (L >= 6 && !stub3_codes.empty()) {
        const int c = encode_hex_at(s, L - 6);
        if (c >= 0)
            for (uint32_t code : stub3_codes)
                if ((uint32_t)c == code) { s.erase(L - 6, 6); break; }
    }
    return ((int)s.size() < min_len) ? std::string() : s;
}
```

### After — iterative 3' stripping + `min_len <= 0` pure-clip

```cpp
std::string clip_adapters(const std::string&           seq,
                          const std::vector<uint32_t>& stub5_codes,
                          const std::vector<uint32_t>& stub3_codes,
                          int                          min_len) {
    std::string s = seq;
    int L = (int)s.size();
    // 5': single hex-code prefix match (unchanged)
    if (L >= 6 && !stub5_codes.empty()) {
        const int c = encode_hex_at(s, 0);
        if (c >= 0)
            for (uint32_t code : stub5_codes)
                if ((uint32_t)c == code) { s.erase(0, 6); break; }
    }
    // 3': ITERATIVE terminal hex-code stripping — promoted from fqdup
    // clip_worker_fn (damage_worker.hpp). Guard mirrors fqdup's `L < 12` break so a
    // 6 bp remnant is never chased below a full hexamer.
    if (!stub3_codes.empty()) {
        bool trimmed;
        do {
            trimmed = false;
            L = (int)s.size();
            if (L < 12) break;
            const int c = encode_hex_at(s, L - 6);
            if (c >= 0)
                for (uint32_t code : stub3_codes)
                    if ((uint32_t)c == code) { s.erase(L - 6, 6); trimmed = true; break; }
        } while (trimmed);
    }
    // min_len <= 0 → "pure clip": return the clipped bytes even when shorter than a
    // floor (profile_reads / the short-fragment path need them). min_len > 0 keeps
    // the legacy drop-to-empty contract for dart's Pass-2 length gate.
    return (min_len > 0 && (int)s.size() < min_len) ? std::string() : s;
}
```

Once this lands, delete the inline `do { … } while (trimmed)` block from
`fqdup/damage_worker.hpp:clip_worker_fn` — but that whole function is deleted in Section 4 anyway (fqdup's
SE full pass moves into `profile_reads`), so this is really only about libtaph owning the canonical loop
and dart's Pass-2 clip picking it up automatically.

---

## 4. REFACTOR `fqdup/src/profile.cpp` — collapse the SE path into `profile_reads()`

> ⚠️ **CORRECTION (2026-07-15) — the "delete the SE pre-scan" step below is NOT byte-identical.**
> The SE pre-scan buffer (`se_scan_buf`) has a SECOND consumer the table omits: it builds
> `prescan_lsd_hist` (`profile.cpp:399-401`) → `prov_edges` (`compute_lsd_edges`, `:408`) → each
> worker's `lbs.configure(prov_edges)` + `lsd_prov_edges` (`:414-415`). And the LSD **fusion** uses
> those provisional edges DIRECTLY as the final fusion edges: `lsd_fuse_edges = states[0].lsd_prov_edges`
> (`profile.cpp:705`). They are never recomputed from the full data — so deleting the SE pre-scan
> changes the LSD bin edges and the fused LSD result is NOT byte-identical.
>
> **Correct scope for §4:** keep fqdup's SE **pre-scan** (it does triple duty: adapter stubs +
> `prov_edges` + the buffered reads), and collapse only the **full pass** onto `profile_reads`, passing
> `precomputed_stubs` so no second stub-scan runs, and a `PerReadHook` that replicates `worker_fn`'s
> ENTIRE side-channel set (`lbs.update`, `lsd_cnt`, `lsd_hist`, `ox_fwd`/`ox_rev`, the read/len
> counters) — NOT just `lsd_hist`. Then validate byte-identity on the 435M-read gold library
> (`LV7008960924-18S-C-eDNA-ds_S41`, gold md5 `ca77a961…`) before committing. This is a separate,
> gold-gated change; it does not block dart (dart pins libtaph, not fqdup).

**Keep the paired / `states_pe` path completely untouched.** Only the single-end / merged (`have_merged`)
pieces collapse. `states` still exists — its `SampleDamageProfile` member is now unused on the SE path
(the hook fills only the side channels), but `lbs`, `lsd_cnt`, `lsd_prov_edges`, `lsd_hist`, `ox_fwd`,
`ox_rev` are still needed for the OxoG / LSD passes downstream.

### Spans that collapse

| fqdup span (verified) | What it did | After |
|---|---|---|
| `266-332` | SE pre-scan (`if (have_merged) { … }`, hand-rolled parallel accumulation + `detect_adapter_stubs`) | **Deleted.** `profile_reads` runs it via `run_prescan`. The `else` (paired) pre-scan **stays**. |
| `385-467` | SE full pass (`if (have_merged) { WorkQueue … worker_fn/clip_worker_fn … }`) | **Deleted.** Replaced by the `profile_reads` call. The `if (have_paired) { … }` block **stays**. |
| `417` (`s.profile.short_fragment_floor = short_floor` in `configure_states`) | floor on per-thread `states` profiles | **Keep** for `states_pe` (paired). On SE `states` it's now dead but harmless; leave it. |
| `527` (`dp.short_fragment_floor = short_floor`) | floor on the merged profile pre-finalize | **Delete** — `profile_reads` sets it on `result.profile` before finalize. |
| `535` (`for (auto& s : states) merge_sample_profiles(dp, s.profile); …`) | merge per-thread profiles + `lsd_hist` + read/len counters | **Replace**: `dp` comes from `profile_reads`; only `merged_lsd_hist` (hook-filled) is merged from `states` here; read/len counters come from `ProfileResult`. |
| `623` (`finalize_sample_profile(dp, n_threads)`) | finalize the SE/merged profile | **Delete** — `profile_reads` finalizes. The **unmerged** `finalize_sample_profile(*dp_unmerged_owner, …)` at the next line **stays**. |

`se_scan_buf` and `reader_se` (the single-pass buffer-reuse optimization, `profile.cpp` around
`247-332` / `385-406`) are **removed**. `profile_reads` re-reads the file for its pre-scan via the reader
factory. This costs one extra sequential decode of the first `prescan_reads` records (fqdup's own comment
notes decode is I/O-bound, ~8% of the pre-scan window); it does not change any result. If you want the
optimization back later, add a `precomputed_stubs`-style hook for reads too — out of scope here.

### Resulting thin SE caller

Replace spans `266-332` (SE branch only), `385-467`, `527`, `535` (SE portion), and `623` (SE finalize)
with:

```cpp
// states still hold the OxoG / LSD side channels; configure_states() (unchanged)
// already set lbs / lsd_cnt / lsd_prov_edges on each. Their .profile member is now
// unused on the SE path (profile_reads owns profiling).
std::vector<uint64_t> merged_lsd_hist(LSD_HIST_BINS, 0);

if (have_merged) {
    taph::ProfileConfig cfg;
    cfg.threads              = n_threads;
    cfg.forced_lib           = forced_lib;
    cfg.min_read_len         = LSD_L_MIN;
    cfg.short_fragment_floor = short_floor;
    cfg.max_read_len         = max_length;         // 0 = off
    cfg.prescan_reads        = static_cast<size_t>(adapter_scan_reads);  // 0 = all

    // Reader factory: each call opens a fresh reader at the file head.
    auto make_reader = [&]() -> taph::BatchFn {
        auto rdr = std::make_shared<std::unique_ptr<FastqReaderBase>>(
            make_fastq_reader(in_path, static_cast<size_t>(n_threads)));
        return [rdr](std::vector<std::string>& out, size_t max) -> bool {
            FastqRecord rec;
            size_t added = 0;
            while (added < max && (*rdr)->read(rec)) {
                out.push_back(std::move(rec.seq));
                ++added;
            }
            return added > 0;
        };
    };

    // Hook: accumulate ONLY the side channels into the caller's per-thread states.
    // profile_reads owns the SampleDamageProfile; we never touch states[tid].profile.
    taph::PerReadHook hook = [&](int tid, const std::string& seq, int L) {
        WorkerState& st = states[tid];
        ++st.lsd_hist[lsd_hist_bin(L)];
        if (!st.lsd_cnt.empty())
            st.lbs.update(seq, (L < LSD_L_MAX) ? L : LSD_L_MAX);
        worker_ox_accumulate(st, seq, L);
    };

    taph::ProfileResult pres = taph::profile_reads(make_reader, cfg, hook);

    // Adopt the finalized profile + counters.
    dp    = std::move(pres.profile);      // dp is the heap-allocated SampleDamageProfile
    stubs = std::move(pres.stubs);
    reads_scanned         = pres.reads_scanned;
    reads_skipped         = pres.reads_skipped;
    len_min               = pres.len_min;
    len_max               = pres.len_max;
    len_sum               = pres.len_sum;
    n_stub5_hits          = pres.n_stub5_hits;
    n_stub3_hits          = pres.n_stub3_hits;
    n_stub_reads_checked  = pres.n_stub_reads_checked;

    // lsd_hist was filled into states via the hook; fold it up here (order-invariant).
    for (auto& s : states)
        for (int hb = 0; hb < LSD_HIST_BINS; ++hb)
            merged_lsd_hist[hb] += s.lsd_hist[hb];
}
```

Everything downstream is unchanged and consumes `dp` / `stubs` / `merged_lsd_hist` / the counters exactly
as before:

- the adapter-prefix F/G/H recompute (`recompute_fgh_excluding_adapter_prefixes`),
- `compute_hex_enriched_5prime`, the OxoG second pass (merges `states[*].ox_fwd/ox_rev`),
- the LSD fusion (merges `states[*].lbs` / `states[*].lsd_cnt`, plus `states_pe` in combined mode),
- the `dp_unmerged_owner` PE finalize and the `length_contrast` JSON block.

Because `dp` is a `std::unique_ptr`-owned `SampleDamageProfile` (`dp_owner`), assign into `*dp_owner`
(`*dp_owner = std::move(pres.profile);`) rather than rebinding the reference — keep the existing ownership
shape.

The pre-scan `else` branch (paired, `profile.cpp:307-332`) and the entire `if (have_paired) { … }` full
pass (`469-517`) are **verbatim unchanged**.

---

## 5. REFACTOR `dart/src/main.cpp` — self-profile via `profile_reads()`

> ✅ **DONE (2026-07-15), branch `feat/libtaph-perread-damage`.** Implemented with two deliberate
> deviations from the plan below, both behavior-preserving:
> 1. **EM-training kept, not dropped.** §5d proposed dropping `em_training_seqs` /
>    `train_periodic_model`. That feeds Pass-2 gene prediction (not the damage gate), so dropping it
>    would silently change ORF output. Instead it is **preserved via a `PerReadHook`** that captures
>    60–200 bp clipped reads into per-thread buffers, flattened in tid order after the pass. This is
>    exactly the side-channel use the hook exists for.
> 2. **Pre-scan reads `prof_in`, not `-i`.** The standalone adapter pre-scan (§5c) now scans the same
>    `prof_in` reads that are profiled (fqdup-faithful: detect stubs on the reads you profile), then
>    applies them to `-i` in Pass 2. Stubs are library-level, so this is consistent with §8.2.
>
> Also fixed a latent contract bug found while wiring `--profile-sample 0`: `run_prescan` treated
> `prescan_reads == 0` as "scan nothing" (`0 < 0`), contradicting the documented "0 = all". Fixed in
> `libtaph/src/damage_profiler.cpp` (map 0 → unbounded). No golden used 0, so no golden moved.

### 5a. Delete Pass-1.5 metamatch (`main.cpp:437-525`)

The whole `if (!library_summary && sample_profile.damage_validated && … d_max_combined > 0.15f) { … }`
block (the second FASTQ pass computing `d_metamatch` and overwriting `sample_profile.d_max_combined`) is
**deleted**. It is a dart-only re-weighting that has no fqdup counterpart and is exactly the kind of
divergence this unification removes. The canonical damage gate is either the fqdup-fit `--library-profile`
sidecar (preferred) or the `profile_reads` self-profile (opt-in).

### 5b. New flags `--profile-input <file>` + `--profile-sample N` (`dart/src/cli/args.hpp` + `args.cpp`)

Per §8.1 the profile INPUT (a file) and the pre-scan COUNT are two distinct flags: dart must profile a
DIFFERENT file (merged pre-dedup reads) than it predicts on (`-i` = dedup reads), or the fix self-defeats.

`args.hpp` — add to `Options`:

```cpp
std::string profile_input;    // FASTQ to profile damage on; empty = -i (input_file)
size_t      profile_sample = 0; // reads sampled for self-profile adapter pre-scan; 0 = all
```

`args.cpp` — in the option parser, and in `print_usage`:

```cpp
else if (arg == "--profile-input" && i + 1 < argc) {
    opts.profile_input = argv[++i];
}
else if (arg == "--profile-sample" && i + 1 < argc) {
    long long v = std::stoll(argv[++i]);
    if (v < 0) throw ParseArgsExit(1, "--profile-sample must be >= 0 (0 = all)");
    opts.profile_sample = static_cast<size_t>(v);
}
```

```
  --profile-input <path>     FASTQ to profile damage on (default: -i). ORF prediction
                             always runs on -i; only the damage profile uses this file.
  --profile-sample N         Reads sampled for self-profile adapter pre-scan
                             (default: 0 = scan all; replaces the old 100k cap)
```

### 5c. 100k cap removal

The standalone Pass-2 adapter pre-scan (`main.cpp`, the `taph::run_prescan(scan_fn, …, 30, 100000)` call)
becomes:

```cpp
auto prescan = taph::run_prescan(
    scan_fn, to_sample_library_type(opts.forced_library_type), 30, opts.profile_sample);
```

i.e. `100000` → `opts.profile_sample` (§8.1 rename; default `0` = all). This pre-scan still feeds `stub5_codes` /
`stub3_codes` for Pass-2 per-read clipping, which stays. Pass the same `adapter_stubs` into the
self-profile config (`precomputed_stubs`) so the self-profile and Pass-2 clip use identical stubs and no
second pre-scan runs.

### 5d. Replacement self-profile branch (replaces `main.cpp:379-436`, the `if (opts.aggregate_damage && opts.use_damage)` Pass-1 body up to `calibrator.initialize`)

```cpp
// Pass 1: damage self-profile (single-engine, via libtaph::profile_reads).
// Only runs when a sidecar was NOT supplied and self-profiling is opted in.
if (opts.aggregate_damage && opts.use_damage) {
    auto pass1_start = std::chrono::high_resolution_clock::now();
    std::cerr << "[Pass 1] Self-profiling (libtaph::profile_reads)...\n";

    taph::ProfileConfig cfg;
    cfg.threads              = num_threads;
    cfg.forced_lib           = to_sample_library_type(opts.forced_library_type);
    cfg.min_read_len         = 30;
    cfg.short_fragment_floor = 30;                 // dart does not lower the floor
    cfg.prescan_reads        = opts.profile_sample; // 0 = all  (§8.1: renamed from profile_reads)
    cfg.precomputed_stubs    = clip_reads ? &adapter_stubs : nullptr;

    // §8.1 BLOCKER: profile the merged pre-dedup reads, NOT -i (the dedup set). Reading
    // opts.input_file here self-defeats the fix (dedup set lands LOW_ABUNDANCE). Use the
    // dedicated profile-input file, falling back to -i only when it is unset.
    const std::string& prof_in = opts.profile_input.empty() ? opts.input_file
                                                            : opts.profile_input;
    auto make_reader = [&]() -> taph::BatchFn {
        auto rdr = std::make_shared<dart::SequenceReader>(prof_in);
        return [rdr](std::vector<std::string>& out, size_t max) -> bool {
            dart::SequenceRecord rec;
            size_t added = 0;
            while (added < max && rdr->read_next(rec)) {
                std::string s = dart::SequenceUtils::clean(rec.sequence);
                out.push_back(std::move(s));
                ++added;
            }
            return added > 0;
        };
    };

    // No hook: dart's self-profile needs only the finalized SampleDamageProfile.
    taph::ProfileResult pres = taph::profile_reads(make_reader, cfg);
    sample_profile = std::move(pres.profile);

    auto pass1_end = std::chrono::high_resolution_clock::now();
    auto pass1_dur = std::chrono::duration_cast<std::chrono::milliseconds>(pass1_end - pass1_start);
    std::cerr << "  Reads: " << sample_profile.n_reads
              << " | 5' damage: " << std::fixed << std::setprecision(1)
              << sample_profile.max_damage_5prime * 100.0f << "%"
              << " | 3' damage: " << sample_profile.max_damage_3prime * 100.0f << "%"
              << " | " << dart::log_utils::format_duration_ms(pass1_dur.count()) << "\n\n";

    damage_model.update_from_sample_profile(sample_profile);
    calibrator.initialize(sample_profile);
}
```

Deleted along with the old body: the manual batch/`omp parallel for` loop, `thread_profiles`, the EM
training-sequence capture (`em_training_seqs`, `em_training_count`) and the
`dart::codon::train_periodic_model(...)` call, the manual `merge_sample_profiles` loop and the manual
`finalize_sample_profile` — all now inside `profile_reads`. If the periodic-model training is still wanted,
re-capture it via a `PerReadHook` in a follow-up; it is **not** part of this unification and is dropped for
now (it did not feed the damage gate).

### 5e. What stays

- `--library-profile` `apply_to` **stays exactly where it is** and remains authoritative. After the new
  Pass-1 branch, keep:
  ```cpp
  if (library_summary) {
      library_summary->apply_to(sample_profile);
      std::cerr << "[library-profile] applied sidecar as authoritative damage gate "
                   "(self-profile superseded)\n";
  }
  ```
  and the later `if (library_summary && !(opts.aggregate_damage && opts.use_damage)) { library_summary->apply_to(...); }`
  fallback for the `--no-aggregate` / `--no-damage` path. The sidecar still wins over the self-profile.
- The standalone Pass-2 adapter pre-scan + `clip_reads` + `stub5_codes`/`stub3_codes`, and every
  `seq = taph::clip_adapters(seq, stub5_codes, stub3_codes, 30)` call inside Pass 2 (gene prediction). Those
  now transparently pick up the iterative 3' stripping from Section 3.
- Pass 2 gene prediction is unchanged.

---

## 6. Verification plan — prove `fqdup profile(reads) == dart self-profile(SAME reads)` byte-identical

**Claim:** with identical `ProfileConfig` and identical reads, the finalized `SampleDamageProfile` is
byte-identical, because both engines execute the same `profile_reads` code over the same integer-count
accumulators, and `finalize_sample_profile` is deterministic.

**Control the inputs that can legitimately differ:**
- Same reads (feed the exact same FASTA/FASTQ; dart applies `SequenceUtils::clean`, so pre-clean the file
  or confirm `clean` is a no-op on it).
- `forced_lib` equal, `short_fragment_floor = 30` on both, `max_read_len = 0` on both.
- **Same `prescan_reads`** on both (else stubs → clip → profile diverge). Use `0` (all) on both:
  fqdup `--adapter-scan-reads 0`, dart `--profile-sample 0`.
- **`-p 1` on both** (removes any `finalize(dp, n_threads)` thread-count question; also makes the run
  trivially reproducible).

**Harness A — libtaph micro-test (fastest signal).** A standalone program linked against libtaph that reads
a FASTQ into a `std::vector<std::string>` once, then profiles it twice through an in-memory
`ReaderFactory`, dumping the fields:

```cpp
// tools/profile_diff.cpp  (link libtaph)
#include "taph/damage_profiler.hpp"
#include <cstdio>
#include <string>
#include <vector>

static std::vector<std::string> load(const char* path); // read seqs, no clean

static taph::BatchFn from_vec(std::shared_ptr<const std::vector<std::string>> v,
                              std::shared_ptr<size_t> pos) {
    *pos = 0;
    return [v, pos](std::vector<std::string>& out, size_t max) -> bool {
        size_t added = 0;
        while (added < max && *pos < v->size()) { out.push_back((*v)[(*pos)++]); ++added; }
        return added > 0;
    };
}

int main(int argc, char** argv) {
    auto seqs = std::make_shared<std::vector<std::string>>(load(argv[1]));
    taph::ProfileConfig cfg;
    cfg.threads = 1; cfg.min_read_len = 30; cfg.short_fragment_floor = 30; cfg.prescan_reads = 0;
    auto mk = [&]() { auto p = std::make_shared<size_t>(0); return from_vec(seqs, p); };
    auto r = taph::profile_reads(mk, cfg);
    const auto& p = r.profile;
    std::printf(
      "{\"d_max_5prime\":%.10g,\"d_max_3prime\":%.10g,\"d_max_combined\":%.10g,"
      "\"lambda_5prime\":%.10g,\"lambda_3prime\":%.10g,"
      "\"fit_baseline_5prime\":%.10g,\"fit_baseline_3prime\":%.10g,"
      "\"pi_point\":%.10g,\"pi_lo\":%.10g,\"pi_hi\":%.10g,\"pi_state\":%d,"
      "\"d_max_se\":%.10g,\"n_reads\":%lld}\n",
      p.d_max_5prime, p.d_max_3prime, p.d_max_combined,
      p.lambda_5prime, p.lambda_3prime, p.fit_baseline_5prime, p.fit_baseline_3prime,
      p.pi.point, p.pi.lo, p.pi.hi, (int)p.pi.state,
      p.bulk_damage.d_max_se, (long long)p.n_reads);   // d_max_se lives on bulk_damage; n_elig_total is not a profile field
    return 0;
}
```

Two runs on the same file must produce byte-identical stdout (self-consistency / determinism gate). Then
build the same binary in the dart tree (dart's `SequenceReader` + `clean`, same input) — its output must
match fqdup's, proving the engine boundary is transparent.

**Harness B — end-to-end, using the shipping JSON.** fqdup already writes the finalized profile:

```bash
fqdup profile -i reads.fq -p 1 --adapter-scan-reads 0 --json fqdup.json
dart predict -i reads.fq -p 1 --profile-sample 0 --allow-self-profile \
    --aggregate-damage --damage-only --summary dart.json    # dart dumps the same fields
# (to profile a DIFFERENT file than -i, add: --profile-input merged.fq.gz  — see §8.1)
```

Compare the load-bearing fields with a tiny diff (exact equality, not tolerance):

```python
# diff_profiles.py
import json, sys
FIELDS = ["d_max_5prime","d_max_3prime","d_max_combined",
          "lambda_5prime","lambda_3prime","fit_baseline_5prime","fit_baseline_3prime",
          "pi_point","pi_lo","pi_hi","pi_state","d_max_se"]   # n_elig_total dropped: not a profile field
a = json.load(open(sys.argv[1])); b = json.load(open(sys.argv[2]))
def get(d,k):           # tolerate nesting differences between the two dumps
    return d.get(k, d.get("damage",{}).get(k))
bad = [(k, get(a,k), get(b,k)) for k in FIELDS if get(a,k) != get(b,k)]
if bad:
    for k,x,y in bad: print(f"DIFF {k}: fqdup={x!r} dart={y!r}")
    sys.exit(1)
print("IDENTICAL")
```

Extend dart's `--damage-only` summary writer (`main.cpp`, the `opts.damage_only` block) to emit the same
field set as Harness A so the comparison is direct. The decisive fields for the original drift are
`pi_state` (the `LOW_ABUNDANCE` vs `DETECTED` split), `pi_point`, `d_max_combined`, and
`d_max_se` (= `p.bulk_damage.d_max_se`) — if those match, the state-gate divergence is gone. The abundance
input behind the `LOW_ABUNDANCE` gate is `p.n_reads` here (there is no `n_elig_total` on the profile); if
the read-LLR eligibility count is the real gate driver, compare it out-of-band from the `read_ancient_llr`
path, not from the profile dump.

If Harness A (pure libtaph, same in-memory reads) is identical but Harness B differs, the remaining
difference is **input** (dart's `SequenceUtils::clean` altering reads, or a different read set reaching the
profiler), not the profiler — bisect there.

---

## 7. Explicitly OUT OF SCOPE

Do **not** touch these in this change:

1. **`finalize_pi` → `combine()` delegation.** The pi finalization internals stay as-is; `profile_reads`
   calls `finalize_sample_profile` and consumes whatever pi it produces. No refactor of the pi combiner.
2. **The paired-end path.** fqdup's `-1/-2` ingestion, `paired_worker_fn`, `PairedWorkQueue`,
   `update_sample_profile_paired`, `states_pe`, the combined-mode `dp_unmerged_owner` finalize, and the
   `length_contrast` JSON block are all unchanged. `profile_reads` is single-end / merged only.
3. **The D2 9-state gate.** No changes to the damage-confidence state machine or its thresholds. The
   unification is purely structural — one loop feeding the existing gate — so the gate must see identical
   inputs and therefore reach identical states across the two engines (that is the whole point of
   Section 6).

---

## 8. Consumer (dart / convergence) side — what the downstream session needs

This section is added by the dart/convergence session. It lists requirements the §5 dart refactor must
satisfy and the cross-session coordination points, so the two sessions don't diverge.

### 8.1 BLOCKER — `--profile-reads` as written cannot do what the pipeline requires

§5b defines `--profile-reads N` as a **count** (prescan sample size) and the self-profile in §5c profiles
`-i`. But the pipeline requirement is: **dart must profile a DIFFERENT read file (the merged pre-dedup
reads) than the one it predicts ORFs on (`-i` = the deduplicated reads).** That is the entire reason the
self-profile matches fqdup — the merged set has the read COUNT that lands `DETECTED`; the dedup set
(17× fewer reads) lands `LOW_ABUNDANCE`. Profiling `-i` self-defeats the fix.

**Needed:** a distinct flag for the profile INPUT FILE, separate from the count. Proposed:

```
  --profile-input <path>   FASTQ to profile damage on (default: -i). ORF prediction
                           always runs on -i; only the damage profile uses this file.
  --profile-sample N       Reads sampled for the self-profile pre-scan (0 = all).
                           (this is the current --profile-reads count; RENAME to avoid
                            the file-vs-count name collision)
```

In §5c, dart builds its `ReaderFactory` over `opts.profile_input` when set, else `-i`. Everything else in
§5c is unchanged. `ProfileConfig::prescan_reads = opts.profile_sample`.

### 8.2 `ProfileResult::stubs` must survive the refactor (dart depends on it)

dart detects adapter stubs on the **profile-input** reads but applies them to the **`-i`** reads in Pass 2
(`dart/src/main.cpp:605`). So `profile_reads()` MUST return the `AdapterStubs` (it already does, via
`ProfileResult::stubs` in §1). dart takes `pr.stubs` and uses it for the `-i` clip. Do not drop `stubs`
from the return in the name of minimalism.

### 8.3 Pin handoff — this is what unblocks the dart rerun

dart pins libtaph: `dart/CMakeLists.txt:42` `LIBTAPH_PINNED_SHA = 14b95f1…`, and the guard is
**FATAL** on any HEAD/tree mismatch.

**DONE (2026-07-15):** `profile_reads()` + the `run_prescan` 0=all fix are committed on libtaph `main`.

- Commit 1 `3054faf` (parent `14b95f1`): "profile_reads: unified SE/merged engine".
- Commit 2 (this session, **supersedes `3054faf` as the pin**): `run_prescan` 0=all fix + docs
  (api.md, changelog.md) + this handoff reconciliation. **Pin to libtaph HEAD** —
  `git -C "$REPOS/libtaph" rev-parse HEAD` on `main`.
- Build + ctest on a compute node: **6/6 pass** (incl. `count_table_golden`,
  `golden_verdict_regression`) at both commits — no golden moved.
- ⚠️ Both commits are **unpushed and untagged** (WIRE-AND-HOLD). dart reaches HEAD only
  via the local libtaph checkout at `…/apps/repos/libtaph` on `main`, which must be **clean**
  (dart's CMake guard FATALs on a dirty tree).

**dart session action:** bump `LIBTAPH_PINNED_SHA` 14b95f1 → `<libtaph HEAD>`, rebuild
`dart/build/`, run the chain. dart will not configure until that bump. (The dart-side §5 code
refactor itself is already applied on `feat/libtaph-perread-damage`.)

### 8.4 Pin the canonical reference read set (byte-identity depends on it)

§6 proves `fqdup profile(reads) == dart self-profile(SAME reads)`. "SAME reads" must be an exact file.
Current fqdup outputs for LV7008956630 are `profile/LV7008956630/{dedup,sorted}.json` — **neither is
named "merged"**, so confirm which read set the canonical fqdup profile used before asserting identity.
Concrete inputs on disk (LV7008956630):

| role | path | size |
|---|---|---|
| damage profile input (merged, high-N) | `…/kapk_ahp56/merged/LV7008956630_merged.fq.gz` | 1.26 GB |
| ORF-prediction target (`-i`, dedup) | `…/kapk_ahp56/derep/LV7008956630/dedup.fq.gz` | 73.7 MB |

Run the §6 harness on whichever file is the pinned reference, both engines at `-p 1`, and diff the
finalized `SampleDamageProfile` fields (`d_max_*`, `pi.point`, `pi.state`, `bulk_damage.d_max_se`).
