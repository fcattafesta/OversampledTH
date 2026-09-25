# Assumptions and practical tips

## Event-level uncertainty

For bin `b`, let `a_g,b` be the sum of row weights from generated event `g` in that bin and let `F` be the oversampling factor. The histogram content is `sum_g a_g,b / F`. Its `Sumw2` variance is `sum_g (a_g,b / F)^2`, because each event contributes one output fill per bin. A `TH1` does not store covariance between different bins.

The factor is constant for the whole action. If oversampling varies by event, precompute an appropriate row weight or use a different action. Negative weights are supported; cancellation within one event is handled before the output fill.

## Which action to choose

- Use `SequentialOversampledHistogram` for one thread and contiguous generated-event IDs. It has the smallest state.
- Use `SlotLocalOversampledHistogram` only when an event stays in one slot **and** appears as one contiguous segment within that slot. Inspect both its cross-slot split and same-slot revisit diagnostics on the actual workload. It keeps seen IDs for those diagnostics, so memory rises with the number of distinct IDs.
- Use `RangeAwareOversampledHistogram` when MT scheduling can split an event at an RDF entry-range boundary. It needs correct entry-range begin/end columns and globally contiguous IDs. Boundary histograms consume memory in proportion to the number of range fragments and histogram bins. Inspect ranges with `python studies/range_layout.py -j 4`.

If a slot-local event is split into contributions `a_1, ..., a_k` for a bin, its reported variance includes `sum_j a_j^2 / F^2`; the correlated event variance is `(sum_j a_j)^2 / F^2`. For same-sign pieces, the slot-local error is too small. Mixed-sign pieces can move the error in either direction. The diagnostics report affected event IDs, not the size of the error difference. A zero cross-slot count alone does not establish exact errors: a same-slot revisit can also split an event's contribution.

## Ordering and filtering

All three actions require the same generated-event ID to identify one logical event. The sequential and slot-local actions flush when the ID changes; a later reappearance of the same ID is treated as another contribution. The range-aware action merges only its deferred boundary fragments, so it also requires globally contiguous IDs. A filter that removes rows preserves order among surviving rows; empty vector values still advance the action's event and range state.

The range action depends on `rdfsampleinfo_.EntryRange()` describing the processing ranges. If that metadata does not match the RDF scheduling, it cannot know which fragments to defer. Use `studies/range_layout.py` to inspect the ranges for a sample before relying on the range-aware result.
