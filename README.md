# flip

Reimplements the algorithm explained by Matthias Müller in C, using SDL3 to render it.

<img alt="preview image of fluid simulation" src="https://github.com/the2nake/flip/blob/master/docs/preview.png?raw=true"  height="600">

## [building]

```bash
cmake -H . -B build
cmake --build build
```

## [usage]

Run the program with command `./ffs` or potentially `./ffs.exe` on Windows.

## [dependencies]

- SDL3
  - package may be listed with `-devel` suffix
- C23 compliant compiler
- CMake
