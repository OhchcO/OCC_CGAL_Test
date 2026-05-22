# Project Suggestions

## High Priority

1. Align the executable entry point.
   - `CMakeLists.txt` currently builds `CgalApp` from `main.cpp`.
   - `main.cpp` has its `main()` wrapped in `#if 0`.
   - The active workflow appears to live in `FileName.cpp`, which has an active
     `main()`.
   - Suggested fix: either switch the CMake target to `FileName.cpp`, or move the
     active entry into a clearly named source such as `src/main.cpp`.

2. Replace hardcoded project paths.
   - `FileName.cpp` currently uses `E:\soft\code\Project1\input\` and
     `E:\soft\code\Project1\output\`.
   - This repository has its own `input/` and `output/` folders.
   - Suggested fix: accept input/output paths through command line arguments,
     CMake definitions, or a small config file.

3. Make the input STEP file configurable.
   - `FileName.cpp` hardcodes `std::string stepfile = "02.stp";`.
   - The current repository input folder does not obviously contain `02.stp`.
   - Suggested fix: use `argv[1]` for the STEP path and default to an existing
     sample in `input/`.

4. Reduce debug-only execution branches.
   - The main loop currently processes only `if (i == 1)`, which means most
     detected split layers are skipped.
   - Suggested fix: replace this with a named debug option such as
     `--slice-index 1`, or process all layers by default.

5. Separate generated artifacts from source-controlled test data.
   - `input/` currently mixes source STEP files, logs, and generated converted
     STEP/PRT files.
   - `output/` contains generated BREP/STEP/PRT files.
   - Suggested fix: keep curated small test inputs in `input/`, move generated
     conversions/results to `output/` or a ignored `runs/` folder, and consider
     ignoring `*.prt`/debug `*.brep` if they are not canonical fixtures.

## Medium Priority

1. Split `FileName.cpp` into modules.
   - Suggested modules:
     - STEP/OCC IO helpers.
     - Z split point detection.
     - Sectioning and edge discretization.
     - 2D loop topology and nesting.
     - CGAL polygon/hull operations.
     - Cavity feature extraction.
     - Debug export helpers.

2. Normalize naming.
   - `FileName.cpp` and `conves_hull_3.cpp` are hard to understand from names.
   - Suggested names: `cavity_slicer.cpp`, `convex_hull_3d_demo.cpp`,
     `convex_hull_2d_demo.cpp`.

3. Remove duplicate includes and repeated link entries.
   - `CMakeLists.txt` repeats several OCC libraries.
   - `FileName.cpp` contains duplicate includes such as `<fstream>` and repeated
     OCC headers.
   - Suggested fix: clean after the build target is stabilized.

4. Add focused regression cases.
   - Useful cases:
     - Simple block with no cavity.
     - Open pocket.
     - Closed cavity.
     - Multiple disconnected cavities.
     - Sloped or curved wall.
     - Very small edges or fragmented section lines.

5. Replace `system("pause")`.
   - It blocks automation and CI.
   - Suggested fix: print final status and let the executable exit normally.

## Lower Priority

1. Convert source comments to a consistent UTF-8 state.
   - Many comments currently display as mojibake in PowerShell.
   - Suggested fix: after code behavior is stable, do a dedicated encoding pass
     and review with an editor configured for UTF-8.

2. Add a short `README.md`.
   - Include dependency versions, expected folder layout, configure/build
     commands, and a minimal run command.

3. Add a small command-line interface.
   - Example:
     - `CgalApp --input input/6_stp.stp --output output --slice all`
     - `CgalApp --input input/6_stp.stp --output output --slice-index 1`

4. Add CMake options for local dependency roots.
   - Example options:
     - `-DOCC_BASE_DIR=...`
     - `-DCGAL_DIR=...`
     - `-DBOOST_ROOT=...`

5. Add output naming conventions.
   - Current debug files are useful but numerous.
   - Suggested fix: use a per-run folder and include model name, slice index,
     and Z value in a consistent format.
