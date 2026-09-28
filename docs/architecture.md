# Architecture

`include/OversampledHistogram.h` defines three ROOT `RDataFrame::Book` actions and shared histogram utilities. Each action accepts an oversampling factor, histogram name/title, and fixed-width 1D binning. Aliases ending in `1D` and `1F` instantiate `TH1D` and `TH1F`. Values may be scalar or `ROOT::VecOps::RVec`; an optional final column supplies a row weight.

For each generated event, a temporary histogram accumulates row weights per bin. The action then fills the output histogram once per event and bin with the accumulated content divided by the oversampling factor. `Sumw2` is enabled on the output histogram so its bin error is computed from these event-level contributions.

The sequential action keeps one current event. The slot-local action keeps an independent current-event histogram, output histogram, and seen-ID set for each RDF slot. At finalization it merges slot histograms and reports IDs that appeared in multiple slots or were revisited after another ID in one slot. The range-aware action additionally tracks each slot's current RDF entry range. Interior events are flushed immediately; the first and last observed event of each range are deferred, merged by event ID, then filled once. It reports how many IDs had more than one boundary piece.

All output histograms are detached from ROOT directories. RDF copies start with empty result and slot state. Per-slot state is independently allocated; final merging happens after processing. The actions rely on ROOT's `RActionImpl` interface (`Exec`, `GetResultPtr`, `InitTask`, and `Finalize`).

Deferred range-boundary events use sorted vectors of nonzero `(bin, content)` contributions, with the piece count stored in each event record. They do not own ROOT histogram clones. Finalization transfers map nodes and their vectors into the global event map, merges matching bin indices, and releases consumed slot storage. It then performs the same native ROOT event-level fills. Sparse totals retain `float` precision for `TH1F` and `double` precision for `TH1D`. See [boundary_memory.md](boundary_memory.md) for scaling, validation, and measurements.

See [assumptions.md](assumptions.md) for correctness limits, [usage.md](usage.md) for booking examples, and [performance.md](performance.md) for measured behavior.
