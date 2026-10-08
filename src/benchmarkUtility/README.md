# Benchmark Utility
↑ [OxySim-129](./../../README.md)

[benchmarkUtility.H](./benchmarkUtility.H)
[benchmarkUtility.C](./benchmarkUtility.C)

## Overview

A small opt-in wall-clock timer for named code sections, used only in `flameletFoam` (not `flameletCloudFoam`):

```cpp
timer.startBenchmark("ThermoUpdate");
// ...
timer.stopBenchmark("ThermoUpdate");
```

Each `stopBenchmark` takes the max duration across processors (`reduce(..., maxOp<float>())`) and accumulates it per key; `printAndResetBenchmark()` prints the running average every `benchmarkInterval` time steps and resets the accumulators.

Controlled from `controlDict`:
```
benchmarkInterval           0;      // 0 disables benchmarking entirely (default)
benchmarkPrintEveryTimeStep true;   // also print each individual duration, not just the average
```

`stopBenchmark` raises a `FatalError` if called for a key that was never `startBenchmark`'d.
