# Project Memory

## Project Snapshot

This repository is a Windows C++17 geometry prototype using OpenCASCADE and CGAL.
The current focus is STEP model reading, Z-layer slicing, 2D loop reconstruction,
convex hull generation, and cavity/open-pocket extraction for machining-related
analysis.

## Main Files

- `CMakeLists.txt`
  - Defines project `OCC_CGAL_Test`.
  - Uses C++17.
  - Hardcodes local dependency paths for Boost, CGAL, GMP/MPFR, and OpenCASCADE.
  - Currently builds target `CgalApp` from `main.cpp` only.

- `FileName.cpp`
  - Current main working file.
  - Contains the active `main()`.
  - Implements STEP loading, face ID mapping, Z split point detection, layer
    slicing, loop reconstruction, convex hull based open-cavity extraction,
    closed-cavity feature extraction, disconnected cavity separation, and debug
    BREP export.
  - Current hardcoded paths point to `E:\soft\code\Project1\input\` and
    `E:\soft\code\Project1\output\`, not this repository's `input` and `output`
    folders.
  - Current input file is hardcoded as `02.stp`.

- `main.cpp`
  - Older 3D convex hull/OCC conversion experiment.
  - Its `main()` is inside `#if 0`, so it is not currently an executable entry.

- `conves_hull_3.cpp`
  - 3D CGAL convex hull experiment.
  - Includes STEP/BREP conversion helpers and a disabled demo `main()`.

- `vector_convex_hull_2.cpp`
  - 2D convex hull to ordered `OneLine` loop conversion helper/demo.

- `test_2d_convex_hull.cpp`
  - Test/demo for filling an open 2D face outline with a convex hull loop.

## Data Folders

- `input/`
  - Contains multiple `.stp`, `.prt`, and `.log` files.
  - Several `*_stp.*` files are currently untracked in git.

- `output/`
  - Contains generated `.brep`, `.stp`, `.prt`, and `.log` artifacts.
  - Intended for visual/debug output from OCC workflows.

- `build/`
  - Generated CMake/Visual Studio build artifacts and ignored by git.

## Current Git State Observed

- `FileName.cpp` has existing uncommitted modifications.
- `input/10_stp_stp.prt`, `input/10_stp_stp.stp`,
  `input/6_stp_stp.prt`, and `input/6_stp_stp.stp` are untracked.
- These existing changes should be treated as user work and not reverted.

## Build Environment Assumptions

- Windows.
- Visual Studio/MSVC generator is likely used.
- OpenCASCADE path:
  `C:/OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0`
- Boost path:
  `E:/soft/boost_1_74_0`
- CGAL path:
  `E:/soft/CGAL/CGAL-5.6.2`
- GMP/MPFR libraries are linked from CGAL auxiliary paths.

## Important Workflow Notes

- Debug output is mostly written as BREP files for inspection in a BREP viewer.
- Many source comments appear mojibake in the terminal, likely due to encoding
  mismatch or previously garbled text; preserve source content unless doing an
  explicit encoding cleanup.
- Avoid changing generated model files in `input/`, `output/`, or `build/`
  unless specifically asked.
- Before compiling the current cavity workflow, align `CMakeLists.txt` with the
  intended entry source file and align hardcoded paths with this repository.
