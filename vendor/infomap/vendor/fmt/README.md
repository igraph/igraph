# {fmt} (vendored, header-only subset)

Vendored subset of [{fmt}](https://github.com/fmtlib/fmt) **11.2.0**.

Only the header-only formatting core is bundled:

- `include/fmt/base.h` — lightweight core (typed-argument machinery,
  `format_string`, `make_format_args`, `formatter`). This is the header
  `src/utils/format_core.h` exposes to widely-included headers like
  `utils/Log.h`. **Note:** since {fmt} 11, the old `fmt/core.h` is only a
  compatibility shim that includes the full `fmt/format.h`, so it is *not*
  vendored — `format_core.h` includes `fmt/base.h` directly to stay lightweight.
- `include/fmt/format.h` — the `fmt::format` renderer.
- `include/fmt/format-inl.h` — header-only definitions.

Deliberately **excluded**: `fmt/ostream.h`, `fmt/os.h`, `fmt/printf.h`,
`fmt/std.h`, `fmt/color.h`. These pull in `<iostream>` / `FILE*` / global
stdio streams, which the R package forbids (see
`scripts/check_cpp_stream_policy.py`). The bundled subset is verified
iostream-free.

## Usage

Do **not** include `<fmt/format.h>` (or `<fmt/base.h>`) directly. Include the
wrapper [`src/utils/format.h`](../../src/utils/format.h) for the full renderer,
or [`src/utils/format_core.h`](../../src/utils/format_core.h) for the
lightweight core; both define `FMT_HEADER_ONLY` consistently across every build
surface (CLI, Python, R, JS/Emscripten). Use `fmt::format(...)` to build
strings; never `fmt::print` (it writes to `stdout`, bypassing the
`infomap::Log` sink).

## Updating

Re-fetch the three headers and `LICENSE` from the pinned tag:

```sh
TAG=11.2.0
base="https://raw.githubusercontent.com/fmtlib/fmt/$TAG/include/fmt"
for f in base.h format.h format-inl.h; do curl -fsSL "$base/$f" -o include/fmt/$f; done
curl -fsSL "https://raw.githubusercontent.com/fmtlib/fmt/$TAG/LICENSE" -o LICENSE
```

When upgrading, re-verify the subset is iostream-free and check whether
`fmt/core.h` is still a shim (drop it if so; point `format_core.h` at whatever
the current lightweight core header is).
