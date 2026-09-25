# Examples

`py/` and `cpp/` each contain one example for the sequential, slot-local, and range-aware histogram actions. The Python examples use [`data/site_sample.json`](../data/site_sample.json) by default and compare histograms with a `fold == 0` baseline; provide `--spec` for another sample. The C++ examples use synthetic data and build with CMake.

```sh
python examples/py/single-threaded.py
python examples/py/slot-local-multi-threaded.py -j 2
python examples/py/range-multi-threaded.py -j 2

cmake --build build --target example_sequential example_slot_local example_range_aware
./build/example_sequential
./build/example_slot_local
./build/example_range_aware
```

See [`docs/usage.md`](../docs/usage.md) for the booking API and entry-range columns.
