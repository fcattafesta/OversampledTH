# OversampledTH

A ROOT RDataFrame action for creating histograms with proper uncertainty handling when performing oversampling in Flashsim.

## Testing

Run the test suite to compare performance and accuracy:

```sh
root -l -q run_all_tests.C
```

## Requirements

The tests require an oversampled dataset which can be found in EOS at:

```sh
/eos/home-f/fcattafe/oversampled_x10_dataset.root
```
