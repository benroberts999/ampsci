# Migration: IO::InputBlockLegacy → IO::InputBlock

Tracking the phased replacement of the custom `IO::InputBlockLegacy` parser with
`IO::InputBlock`, a wrapper around `nlohmann::json` (already bundled in `src/json/`).

## Why JSON?

- nlohmann/json 3.11.3 already bundled (`src/json/json.hpp`)
- `//` and `/* */` comments supported via `ignore_comments=true`
- First-class nested objects (no TOML `[[array]]` inline-table limitation)
- Excellent VSCode support
- Python `json.load()` works with no dependencies

File extension: `.jsonc` (recommended) or `.json`.

## Input file format

```jsonc
// ampsci input file
{
  "Atom": { "Z": "Cs" },
  "HartreeFock": { "core": "[Xe]", "valence": "6sp5d4f" },
  "Grid": { "r0": 1e-6, "rmax": 120.0, "num_points": 1000 },
  "Basis": { "number": 40, "states": "20spdf" },
  "Module": [
    { "type": "matrixElements", "operator": "E1", "rpa": true },
    { "type": "matrixElements", "operator": "M1", "rpa": true },
    {
      "type": "matrixElements",
      "operator": "hfs",
      "rpa": false,
      "options": { "F": "SingleParticle" }
    }
  ]
}
```

Modules are a JSON array; each entry has `"type"` giving the module name.
All other keys are passed directly to the module as `IO::InputBlock`.

## Phase status

### Phase 1 -- IO::InputBlock class [DONE]
- `src/IO/InputBlock.hpp` -- header-only JSON wrapper
- `src/IO/InputBlock.tests.cpp` -- unit tests

### Phase 2 -- Top-level parsing [DONE]
- `main.cpp`: `.json`/`.jsonc` extension detected, parsed via `InputBlock::from_file()`
- `ampsci.cpp`: JSON path still converts to `IO::InputBlock` via `to_ampsci_string()` (bridge, see Phase 4)
- `Modules::runModules2()`: dispatches `IO::InputBlock` directly to module functions

### Phase 3 -- Module functions [DONE]
All module functions now take `const IO::InputBlock &`:

- `src/Modules/basic.cpp`
- `src/Modules/Breit.cpp`
- `src/Modules/dcp.cpp`
- `src/Modules/exampleModule.cpp`
- `src/Modules/HFAnomaly.cpp`
- `src/Modules/isotopeShift.cpp`
- `src/Modules/ladder.cpp`
- `src/Modules/matrixElements.cpp`
- `src/Modules/pnc.cpp`
- `src/Modules/polarisability.cpp`
- `src/Modules/qed.cpp`
- `src/Kionisation/Module_Kionisation.cpp`

`ModuleFn` typedef updated. `runModules` reads `"Module"` JSON array.
`DiracOperator::generate()` and `FactoryFn` use `IO::InputBlock` directly.

### Phase 4 -- Main solver and operators [DONE]

All remaining code migrated to `IO::InputBlock`:

- `src/ampsci/ampsci.cpp` -- main solver driver
- `src/main.cpp` -- `.in` path parses via `IO::InputBlockLegacy` + `IO::from_legacy()`; `-o`/`-i` help queries use `IO::InputBlock` directly; `.json`/`.jsonc` path uses `IO::InputBlock::from_file()` directly
- `src/DiracOperator/GenerateOperator.hpp/.cpp` -- single `generate()` overload, single `FactoryFn`
- `src/Wavefunction/Wavefunction.hpp/.cpp` -- `radiativePotential`, `ConfigurationInteraction`, `formSigma`
- `src/Wavefunction/BSplineBasis.hpp/.cpp` -- `Parameters(IO::InputBlock)`
- `src/CI/ConfigurationInteraction.cpp` -- `configuration_interaction(IO::InputBlock)`
- `src/Potentials/RadPot.hpp/.cpp` -- `radiativePotential(IO::InputBlock)`
- All operator `*.hpp` -- `static generate(IO::InputBlock, Wavefunction)` factories

`IO::InputBlockLegacy` is retained only for `.in` file parsing (via `IO::from_legacy` bridge).

All 118 tests pass.

## IO::InputBlock API reference

```cpp
// Construction
InputBlock b{"name"};
InputBlock b{"name", json_node};
auto b = InputBlock::from_file("path/to/file.jsonc");  // parses JSONC

// Value access (same as InputBlock)
auto x = b.get("key", default_value);   // T inferred from default
auto x = b.get<T>("key");               // returns std::optional<T>

// Block access
auto sub = b.getBlock("SubBlock");      // std::optional<InputBlock>
auto sub = b.get_block("SubBlock");     // InputBlock (empty if absent)
bool ok  = b.has_block("SubBlock");
bool ok  = b.has_option("key");

// Dynamic construction (replaces add(IO::Option{key, to_string(val)}))
b.set("key", value);                    // any JSON-serialisable type
b.set_block("name", other_block);

// Validation
b.check({{"key", "description"}, {"Sub{}", "block desc"}});

// Escape hatch
const nlohmann::json &j = b.node();
```
