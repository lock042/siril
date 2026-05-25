# Plan: drop endianness conversions in the sirilpy ↔ Siril protocol

Branch: `uv` (same branch as the rest of the uv conversion work; this
is a closely related cleanup)
Author: Adrian Knagg-Baugh
Date: 2026-05-26

## 0. Progress tracking

Every discrete unit of work has a `[ ]` checkbox. When an item is
completed, the checkbox **MUST** be updated to `[x]` in the same
commit that lands the work (same rule as the uv conversion plan).
A `[~]` is used for "in progress" / partial.

## 1. Goal and scope

### Goal

Two cleanups landing together, both motivated by the same insight
(sirilpy ↔ Siril is same-machine IPC, not cross-architecture):

1. **Remove every byte-order conversion** (network byte order /
   big-endian) from the protocol. The two sides share endianness;
   the BE conversions are pure dead weight, plus an inconsistency
   hazard (some fields already aren't converted —
   `shared_memory_info_t` and SHM payloads being the obvious
   examples).

2. **Remove the "widen to 64 bits" wire-format normalisation** in
   the `*_to_py` deep-struct serialisers (FITS keywords, imstats,
   homography, distodata, imgdata, psfstar, analysis, fits). The
   widening was originally added to make endianness conversion
   uniform across types (`siril_pythoncommands.c:125`); without
   the byte swap the widening serves no purpose, and dropping it
   shrinks the wire and lets the C macros collapse to one-liners.

After this work, every numeric field is transmitted **at its
native width, in native byte order**. The wire format becomes
"whatever the C struct already is, copied straight to the socket"
— which is what `shared_memory_info_t` and the SHM payloads do
already, so we end up with one consistent rule across the whole
protocol.

### Non-goals

- Cross-architecture IPC. We deliberately abandon any pretense of
  supporting Siril and the python child on different machines or
  architectures. Same process tree, same machine, same endianness,
  same widths; the protocol stops worrying about anything else.
- Wire-format stability across Siril versions. sirilpy is shipped
  with Siril (base venv + per-script venvs sync to the on-disk
  module version on every launch via §4.3 / §4.5 of the uv
  conversion plan), so the two sides always move in lockstep. No
  version negotiation, no fallback.
- Anything beyond endianness + widening: padding within a single
  field, alignment, field layout reordering — out of scope.

## 2. Current state

### Counts

- `src/io/siril_pythonmodule.c`: 16 BE-conversion sites
- `src/io/siril_pythoncommands.c`: 109 BE-conversion sites (mostly
  via the three `BE64` macros at lines 67–97)
- `python_module/sirilpy/connection.py`: 68 `struct.pack/unpack`
  format strings starting with `!`

### Inconsistencies in the existing code (motivation for the cleanup)

1. **`CommandHeader.length` (uint32) and `ResponseHeader.length`
   (uint32)** are BE-converted at send + parse. Python sees them via
   `struct.pack('!Bi', …)` / `struct.unpack('!BI', …)`.

2. **`incoming_image_info_t`** (width, height, channels, data_type,
   size, shm_name) — fields ARE BE-converted on parse in
   `handle_set_pixeldata_request`, `handle_plot_request`,
   `handle_save_image_file_request`, etc. (see lines 613–942 of
   `siril_pythonmodule.c`).

3. **`shared_memory_info_t`** (size, data_type, width, height,
   channels, shm_name) — fields are NOT byte-converted. The C side
   does `memcpy(&info, payload, sizeof(info))` and ships the raw
   bytes; the Python side reads them back with ctypes (which on x86
   means little-endian). This works on every supported platform
   because every supported platform is little-endian, but it is
   inconsistent with (1) and (2).

4. **SHM payloads** (the actual pixel buffers, FITS header byte
   blobs, ICC profile bytes, history strings, etc.) are not byte-
   converted either. We have no need to convert numpy `float32` /
   `uint16` arrays element-by-element — they're consumed in native
   order by both ends.

5. **Per-command structured payloads** in `siril_pythoncommands.c`
   (rectangles, vports, polygon IDs, frame indices, channel
   selectors, plot configurations, keyword dumps, …) are
   individually BE-converted via the `GUINT32_TO_BE` /
   `GUINT32_FROM_BE` / `GUINT64_TO_BE` / `GUINT64_FROM_BE` family
   and the three `BE64` macros at the top of the file. Mirrored on
   the Python side by `struct.pack('!I…')` / `struct.unpack('!I…')`.

The "BE everywhere except SHM" pattern is the worst of both worlds:
slower code paths than necessary AND an interop bug waiting to
happen the day someone adds a non-x86 build target (we'd notice
arrays full of garbage, not headers full of garbage — the SHM bug
is silently latent).

## 3. Approach

Replace every `BE`-flavoured conversion with a native operation,
AND drop the explicit widening casts at C call sites. On both
sides the change is mechanical:

| Before | After |
|---|---|
| `GUINT32_TO_BE(x)` etc. | (drop — write `x` directly) |
| `GUINT32_FROM_BE(x)` etc. | (drop — read `x` directly) |
| `COPY_BE64(val, type)` macro | new `COPY_FIELD(val, type)` macro: bounds check + `memcpy(ptr, &(type){val}, sizeof(type))` + ptr advance. NO widening, NO byte swap — `type` is now just "the on-wire width" which equals the C-struct native type. |
| `COPY_BE64((uint64_t) fit->keywords.lo, uint64_t)` | `COPY_FIELD(fit->keywords.lo, uint16_t)` (or whatever the actual C-side native type is — see §4.2 for the audit) |
| `FROM_BE64_INTO(dest, val, type)` | (eliminate macro — call sites become `memcpy(&dest, &val, sizeof(type))` or direct assignment) |
| `TO_BE64_INTO(dest, val, type)` | (eliminate macro — call sites become `dest = val`) |
| `struct.pack('!Bi', …)` (Python) | `struct.pack('=Bi', …)` |
| `struct.unpack('!BI', …)` | `struct.unpack('=BI', …)` |
| `'Q', # lo (uint64 padded)` in models.py | `'H', # lo (uint16)` — actual C native type |
| `'Q', # stackcnt (uint64)` (currently widened) | C-native width (likely `'I'` for uint32) |

Two format-character rules in Python `struct`:

1. `=` means **native byte order, standard sizes, no alignment**.
   That matches `__attribute__((packed))` on the C side, and the
   field-by-field `memcpy(ptr, &val, sizeof(val)); ptr += sizeof(val);`
   pattern that all the `*_to_py` functions use.
   Critically, `=` is NOT `@` (which uses **native sizes with
   alignment**). Using `@` would silently insert padding on some
   platforms. `=` is the right choice.

2. Each format character now has to match the **actual C-side
   native type**. The current `'Q'` (uint64) for fields that the
   C-side widens from a smaller native type becomes the smaller
   type's matching char (`'H'` for uint16, `'I'` for uint32,
   `'h'` / `'i'` for signed, `'?'` for bool 1 byte, etc.). This is
   what makes the audit in §4.2 below the most delicate part of
   the change.

The wire becomes simpler and smaller. The FITS keywords block
(currently ~16 padded numerics × 8 bytes = ~128 numeric bytes plus
strings) shrinks proportionally — exact savings depend on the
field-type audit, but a rough estimate is 40–60 bytes per keyword
block plus similar savings in imstats / homography / etc.

## 4. Work items

### 4.1 C macros in `siril_pythoncommands.c`

The three existing macros do **two** things: byte-swap to BE (what
we're removing) AND widen smaller numeric types to 8 bytes on the
wire (also what we're removing). After this change, all three
macros either collapse to trivial one-liners or shrink dramatically.

Usage audit (count of distinct types ever passed):

| Macro | `double` | `int64_t` | `uint64_t` | `size_t` | `type` (generic) | Total |
|---|---|---|---|---|---|---|
| `COPY_BE64` | 86 | 10 | 0 | 0 | 1 | **97 call sites** |
| `FROM_BE64_INTO` | 12 | 0 | 0 | 0 | 1 | **13 call sites** |
| `TO_BE64_INTO` | 10 | 0 | 0 | 6 | 1 | **17 call sites** |

All current type parameters are 8-byte. After we remove the
widening, the `type` parameter at each `COPY_*` call site becomes
"the actual C-struct native type for this field" — frequently
smaller (uint16_t, uint32_t, gint, gboolean) than the current
uniform 8 bytes. The per-call-site changes are spelled out in
§4.2 + §4.2bis below.

#### Macros to KEEP (simplified, renamed)

- [ ] Replace `COPY_BE64` macro definition (lines 67–79) with
      `COPY_FIELD` — drop the union + `GUINT64_TO_BE`. The `type`
      parameter stays (it gives the macro a uniform expression-or-
      lvalue interface — see "why the temporary" below). New body:
      ```c
      #define COPY_FIELD(val, type) \
          do { \
              type _tmp = (val); \
              if ((ptr + sizeof(_tmp)) - start_ptr > (size_t)maxlen) { \
                  siril_log_debug("Error: buffer overflow at %s\n", #val); \
                  return 1; \
              } \
              memcpy(ptr, &_tmp, sizeof(_tmp)); \
              ptr += sizeof(_tmp); \
          } while (0)
      ```
      *Why the temporary*: existing call sites include literal
      constants (`COPY_BE64(1.0, double)` at lines 129–130 etc.)
      which are rvalues — you can't take their address. The
      `type _tmp = (val)` step handles literals, casts, lvalues,
      and rvalues uniformly.
- [ ] Rename all 97 `COPY_BE64(…)` call sites in this file to
      `COPY_FIELD(…)`. Each call site's `type` argument is
      individually re-audited in §4.2bis.

#### Macros to ELIMINATE entirely

- [ ] Delete `FROM_BE64_INTO` macro (lines 81–87). At each of its
      13 call sites, replace
      `FROM_BE64_INTO(dest, val, type)` with
      `memcpy(&(dest), &(val), sizeof(type));` (or, where both are
      the same type already, plain `(dest) = (val);`).
- [ ] Delete `TO_BE64_INTO` macro (lines 89–97). At each of its
      17 call sites, replace `TO_BE64_INTO(dest, val, type)` with
      direct assignment `(dest) = (val);` — the union dance was
      entirely a byte-swap mechanism with no other side effect.

#### Wire-format comment

- [ ] Update the comment block above the macros and the line-125
      "All types shorter than 64bit are converted to 64bit types
      before endianness conversion and transmission" comment.
      Replace with a one-liner along the lines of "Each field is
      written at its native C type and width, in native byte
      order. Both ends share a machine, so neither endianness nor
      width normalisation is needed."

### 4.2bis Field-width audit (NEW — the careful part of the work)

Because the widening is going away, every `COPY_BE64((cast)X, T)`
call needs to be re-cast to use X's **actual C-struct type**, and
the matching Python format character in `models.py` (and the few
struct-format strings in `connection.py` that parse bulk-struct
returns) needs to change to match.

The risk: an off-by-one on Python format characters → silent
misalignment → wrong numbers in all later fields of the same
block. Mitigation: do the audit field-by-field, with both the C
source (`siril.h` / `fits_keywords.h` / etc.) and the Python
`_KEYWORD_FORMAT_PARTS` table side by side, and validate via a
round-trip integration test (§4.9) on every batch.

For each `*_to_py` function (the C-side serialisers) and its
matching Python deserialiser:

- [ ] **`keywords_to_py`** in `siril_pythoncommands.c:101+` ↔
      `FKeywords._KEYWORD_FORMAT_PARTS` in `models.py:174–233`.
      Field-by-field audit. The Python comments
      (`# lo (uint64 padded)`, `# binning_x (uint64)`, etc.) flag
      which fields are currently widened and need re-typing.
      Approximate work: ~25 numeric fields. The string fields
      (`{FLEN_VALUE}s`) are unchanged.
- [ ] **`imstats_to_py`** in `siril_pythoncommands.c` ↔ the
      `ImageStats` deserialiser in `models.py`. Audit pending.
- [ ] **`homography_to_py`** ↔ `Homography` deserialiser. Audit
      pending.
- [ ] **`distodata_to_py`** ↔ `DistoData` deserialiser. Audit
      pending.
- [ ] **`imgdata_to_py`** ↔ `ImgData` deserialiser. Audit pending.
- [ ] **`psfstar_to_py`** ↔ `PSFStar` deserialiser. Audit pending.
- [ ] **`analysis_to_py`** ↔ wherever the analysis tuple is parsed.
- [ ] **`fits_to_py`** ↔ `FFit` deserialiser. This is the
      composite of several sub-blocks — verify the offsets still
      add up after width changes.

Mechanical check after each: `Python.struct.calcsize(FORMAT)`
should equal the byte count emitted by the C side for the same
block. The C side has `_SIZE` constants in headers; verify these
get updated where they exist, otherwise compute via
`sizeof(structure)` style helpers.

Suggested commit per `_to_py` function — keeps each batch
small and the round-trip test surface tight.

### 4.2 Per-command BE conversions in `siril_pythoncommands.c`

Drop every `GUINT32_TO_BE` / `GUINT32_FROM_BE` / `GINT32_TO_BE` /
`GINT32_FROM_BE` / `GUINT64_TO_BE` / `GUINT64_FROM_BE` /
`GUINT16_TO_BE` / `GUINT16_FROM_BE` call. After dropping, fields
become plain integer reads/writes via `memcpy` or direct struct
field access — no helper needed.

Rough groupings to tackle one at a time (helps reviewability — each
group is a coherent batch of related lines):

- [x] `CMD_GET_CONFIG` value packing (lines 439, 446)
- [x] `CMD_PLOT` data parsing (lines 542, 587, 600)
- [x] `CMD_REQUEST_SHM` and SHM-info parsing (line 666)
- [x] `CMD_GET_DIMENSIONS` reply (lines 687–689)
- [x] `CMD_GET_SELECTION` / `CMD_SET_SELECTION` rectangles
      (lines 723–726, 772–775, 978–981, 1055–1058)
- [x] `CMD_GET_ACTIVE_VPORT` reply (line 752)
- [x] `CMD_GET_STAR_IN_SELECTION` / `CMD_GET_STATS_FOR_SELECTION`
      arguments (lines 839–844, 997, 1000)
- [x] Boolean / integer status returns (lines 1096, 1102, 1108)
- [x] Polygon ID encoding (line 1554) — done in both
      `siril_pythonmodule.c` and `gui/user_polygons.c`
- [x] All other remaining call sites — `grep` confirms zero
      `GUINT*_*_BE` / `GINT*_*_BE` calls left in
      `siril_pythoncommands.c` and `siril_pythonmodule.c`
- [x] `gui/user_polygons.c`: `g_htonl` / `g_ntohl` /
      `FROM_BE64_INTO` for polygon header and per-point doubles
      (added as part of this batch since polygon serialise is
      naturally inline payload, not bulk struct)
- [x] Python: 48 `struct.(pack|unpack)('!…')` sites in
      `connection.py` flipped to `=`. Plus polygon serialise/
      deserialise in `models.py` (8 sites) and `plot.py` (14
      sites) — these are inline payload paths, not bulk struct,
      so they belong here even though they live outside
      `connection.py`.

### 4.3 C-side header send / receive

- [ ] `send_response()` in `siril_pythonmodule.c:128` — drop the
      `GUINT32_TO_BE` wrap on `ResponseHeader.length`.
- [ ] `handle_client_communication()` (or whichever function reads
      the CommandHeader at line ~3000-ish) — drop the corresponding
      `GUINT32_FROM_BE`. (Need to verify the actual site.)
- [ ] Header type comments: update inline `// Convert to network
      byte order` style comments to reflect the new reality.

### 4.4 `incoming_image_info_t` parsing in `siril_pythonmodule.c`

- [ ] `handle_set_pixeldata_request` (lines 613–623) — drop the
      `FROM_BE` on width/height/channels/size/data_type.
- [ ] `handle_plot_request` (lines 790–800) — drop equivalents.
- [ ] `handle_save_image_file_request` (lines 937–942) — drop
      equivalents.
- [ ] Other `handle_*_request` functions taking
      `incoming_image_info_t` — sweep for missed sites.

### 4.5 `shared_memory_info_t` — confirm no change needed

- [x] **Already native.** `siril_pythonmodule.c` ships the struct
      raw via socket; Python reads it via `ctypes.Structure` (which
      is native-ordered). Verified by inspection; this item is
      pre-ticked because there is nothing to do.

### 4.6 SHM payloads — confirm no change needed

- [x] **Already native.** Pixel buffers, FITS header bytes, ICC
      profile bytes, history strings, etc. are passed raw via
      shared memory. No conversion existed before; none added.

### 4.7 `python_module/sirilpy/connection.py` (and `models.py`)

Two distinct kinds of change:

#### (a) Inline command/response payloads — endianness only

These are short, hand-coded `struct.pack/unpack` calls for things
like rectangles, vport numbers, polygon IDs, dimensions. The C
side packs them field-by-field at native width already (no
widening), so the only change here is `!` → `=`. The 68
`struct.pack('!...')` / `struct.unpack('!...')` format strings
fall into this bucket and need only an endianness flip. Group by
call site type for reviewability:

- [ ] Header send/recv (`'!Bi'`, `'!BI'`) — lines ~329, 362, 388
- [ ] Size headers (`'!Q'`) — line ~528
- [ ] Progress / float payloads (`'!f'`) — line ~1096
- [ ] Rectangles + region selectors (`'!IIII'`, `'!IIIII'`,
      `'!IIIIII'`) — many lines (~1212, 1237, 1295, 1363, 1425,
      1427, 1611, 1611, 1753, etc.)
- [ ] Single ints (`'!I'`, `'!i'`) — many lines (~1271, 1432,
      2265, 2747, 2770, 2820, 2851, 2962, 2993, 3547, 3655, 4154,
      4221)
- [ ] Two doubles (`'!2d'`) — lines ~1473, 1530
- [ ] Mixed bool + int combinations (`'!??'`, `'!I???'`, `'!???'`,
      `'!I?'`) — lines ~1600, 1741, 1742, 3267, 4191, 4444
- [ ] Multi-int with strings (`'!II'`) — lines ~2890, 2927
- [ ] Bool / int parse on response (`'!I'`, `'!i'`, `'!d'`) —
      lines ~3780, 3782, 3784, 4121
- [ ] Anything else `grep` finds after the batches above — final
      sweep to confirm `grep -E "struct\.(pack|unpack)\([^,]*['\"]\!"
      connection.py` returns nothing.
- [ ] Update inline comments that say "network byte order" to say
      "native byte order".

#### (b) Bulk-struct deserialisers — endianness AND field widths

These are the `KEYWORDS_FORMAT` / similar `ClassVar` format
strings in `models.py` that mirror the C-side `*_to_py`
serialisers covered by §4.2bis. For each, both the leading `!`
needs to flip to `=` AND every per-field format character may
need to change from a padded `'Q'` / `'q'` to the actual native
type (`'H'` for uint16, `'I'` for uint32, `'?'` for bool, etc.).

The audit and work items are owned by §4.2bis — listed here only
to clarify that some of the Python-side format-string edits live
in `models.py`, not `connection.py`. A `grep` for `'!'` returning
zero across BOTH files is the completion check.

### 4.7bis Consolidate serialize/deserialize onto model classes

While auditing the bulk-struct format strings (§4.2bis), several
classes were noticed that either (a) lack `serialize` /
`deserialize` methods entirely, or (b) have a format string
duplicated inline at the connection.py call site rather than
declared as a `ClassVar` on the class.

Existing classes already follow the pattern we want (`ImageStats`,
`FKeywords`, `BGSample`, `PSFStar`, `RegData`, `ImgData`,
`Polygon`, `ImageAnalysis`) — each owns its `_FORMAT` /
`_SIZE` `ClassVar`s and exposes a `deserialize(cls, data)` (plus
a `serialize(self)` where the write-direction is used). The
goal is to bring the stragglers up to the same shape so that the
field-type audit in §4.2bis has exactly one place per class
where the format string is declared, rather than ten.

This is not strictly required for the endianness / width work,
but doing it now (rather than later) saves an audit-on-top-of-
an-audit because the format strings change in §4.2bis anyway —
we'd be touching the same lines twice otherwise.

Work items (one commit each — keeps the change reviewable):

- [ ] **`FFit`** — currently assembled inline in
      `connection.get_image()` (around line 3201) and
      `connection.get_seq_frame()` (around line 3441), with the
      13-field `format_parts` list duplicated. Add
      `FFit._CORE_FORMAT` / `_CORE_SIZE` `ClassVar`s holding the
      core 13-field block (rx/ry/naxes/bitpix/orig_bitpix/checksum/
      mini/maxi/neg_ratio/top_down/focalkey/pixelkey/color_managed),
      a `FFit.deserialize_core(cls, data) -> values_dict` (or
      similar — returns the parsed primitives so the caller can
      assemble the full FFit with sub-blocks). Replace both
      `get_image` and `get_seq_frame` to call it.
- [ ] **`Homography`** — currently has no methods. Mirrors C's
      `homography_to_py` (lines around siril_pythoncommands.c
      where `Homography` fields get packed). Add
      `_FORMAT` / `_SIZE` ClassVars and `deserialize`. Find and
      retire any inline unpacking.
- [ ] **`DistoData`** — no methods today. Mirrors C's
      `distodata_to_py`. Add the same pair. Find call site(s).
- [ ] **`Sequence`** — no methods today. Audit C's `seq_to_py`
      (if it exists) / wherever sequence data gets serialized.
      Add `serialize` / `deserialize` to match.
- [ ] **`FPoint`** — has `serialize` but no `deserialize`. Add
      the missing half.
- [ ] **`_SharedMemoryInfo`** (in `shm.py`) — a `ctypes.Structure`
      with no parsing helpers; the 6-tuple `struct.unpack_from`
      pattern is copy-pasted three times inside
      `get_seq_frame()` alone. Add a
      `_SharedMemoryInfo.from_buffer(cls, buf, offset=0)`
      classmethod that wraps the `struct.unpack_from` + the
      kwarg-style constructor, returns `(info, new_offset)`. The
      class is already ctypes-native-ordered so no endianness
      issue, but the consolidation removes 30+ lines of
      repetition.
- [ ] **General sweep**: `grep -nE "struct\.unpack" connection.py`
      and check for any other inline parsers that duplicate
      logic owned (or that should be owned) by a model class.

Each addition follows the existing idiom on `FKeywords`:

```python
@dataclass
class Foo:
    field1: int = 0
    field2: float = 0.0
    ...

    _FORMAT_PARTS: ClassVar[List[str]] = [
        'i',  # field1 (matched against the C native type per §4.2bis)
        'd',  # field2
        ...
    ]
    _FORMAT: ClassVar[str] = '=' + ''.join(_FORMAT_PARTS)
    _SIZE: ClassVar[int] = struct.calcsize(_FORMAT)

    @classmethod
    def deserialize(cls, data: bytes) -> 'Foo':
        if len(data) != cls._SIZE:
            raise ValueError(...)
        f1, f2, ... = struct.unpack(cls._FORMAT, data)
        return cls(field1=f1, field2=f2, ...)

    def serialize(self) -> bytes:  # only when needed
        return struct.pack(cls._FORMAT, self.field1, self.field2, ...)
```

After this section is complete, the only `struct.(un)pack` calls
left in `connection.py` should be the short inline ones for
commands that don't carry a structured payload (e.g. a single
`int` reply). Anything with two or more fields belongs on a
model class.

### 4.8 Bump sirilpy version

- [ ] `python_module/pyproject.toml`: `1.1.21` → `1.1.22`.
- [ ] Document in plan §6 (rollout) that this is a wire-protocol
      change requiring matched C side. Because sirilpy is installed
      from the on-disk module on every base-venv init and per-script
      venv freshness check (uv conversion §4.3 / §4.5 with surgical
      update from Option 3), the two sides re-synchronise on the
      first Siril launch after upgrade. No old-sirilpy-meets-new-
      Siril scenario is possible in normal use.
- [ ] Scripts that explicitly check `check_module_version` for
      something below 1.1.22 are unaffected — the helper-API surface
      is unchanged, only the wire format underneath is different,
      and a script never talks to the wire directly.

### 4.9 Test plan

Mechanical checks (do all of these after each commit batch above
that lands on `master` / `uv`):

- [ ] `grep -E "GUINT[136][26]_TO_BE|GUINT[136][26]_FROM_BE|GINT[136][26]_TO_BE|GINT[136][26]_FROM_BE|GUINT64_TO_BE|GUINT64_FROM_BE" src/io/siril_python*.c` returns zero results.
- [ ] `grep -E "struct\.(pack|unpack)\([^,]*['\"]\!" python_module/sirilpy/connection.py python_module/sirilpy/models.py` returns zero results.
- [ ] `grep -E "COPY_BE64|FROM_BE64_INTO|TO_BE64_INTO" src/io/siril_python*.c` returns zero results (all macros renamed or deleted).
- [ ] After §4.7bis: every model class with a non-trivial wire
      format has both `_FORMAT` / `_SIZE` ClassVars and a
      `deserialize(cls, data)` (plus `serialize(self)` where
      send-direction is used). `connection.py` should have no
      `struct.unpack` call site with 2+ fields that bypasses a
      model class — checked by reading the diff, since this is a
      structural rather than grep-able property.
- [ ] No new compiler warnings.
- [ ] **Wire-size sanity**: after each `*_to_py` audit (§4.2bis),
      assert that `struct.calcsize(FORMAT)` in the Python class
      equals the actual size emitted by the C function. Easiest
      via the existing `_SIZE` constants where they exist; for
      new ones, log the byte count at the end of each `*_to_py`
      and compare.

Functional checks (run via `siril-cli -s -` with a small pyscript
that exercises each command class):

- [ ] **Headers**: any pyscript that runs at all (the command/
      response dispatch goes through them on every roundtrip).
- [ ] **Single-int reads/writes**: `SirilInterface.get_active_vport`,
      `SirilInterface.get_selection`, `set_selection`.
- [ ] **Dimensions / image metadata**: `get_dimensions()`,
      `is_image_loaded()`, `is_sequence_loaded()`.
- [ ] **Pixel data round-trip**: `get_pixeldata()` →
      `set_pixeldata()` with a small fits.
- [ ] **FITS header / keywords**: `get_fits_header()`,
      `get_keywords()`, `set_image_header()`. These exercise the
      `COPY_NATIVE64` / `FROM_NATIVE64_INTO` macros.
- [ ] **Float payloads**: `update_progress(0.5)`.
- [ ] **Double pairs**: `pix2wcs(x, y)`, `wcs2pix(ra, dec)`.
- [ ] **SHM transfer**: any operation that uses
      `request_shm` (large header, pixel data path, ICC profile).
      Confirms SHM-info struct + SHM payloads still work (they were
      already native, so this is just a regression check).
- [ ] **Polygons**: `add_user_polygon(...)`, `delete_user_polygon`.
      Exercises polygon-ID encoding (line 1554 in the old code).

Cross-test against a script that uses the new helpers AND old code
paths to make sure nothing was missed. A quick script that calls
`get_pixeldata()`, `set_pixeldata()`, `update_progress()`, and
`get_keywords()` covers most of the surface.

### 4.10 Commit strategy

Twelve-ish commits in this order, sized so each one is reviewable
and can be smoke-tested independently:

1. **Macros simplified + headers** — change the three C macros in
   `siril_pythoncommands.c`, the CommandHeader/ResponseHeader
   send/recv, and the two corresponding Python lines for
   `'!Bi'`/`'!BI'`. The header path is exercised by every
   roundtrip so this commit shakes out any base-level mistakes
   immediately. **[done — `384c6649d`]**
2. **`incoming_image_info_t`** — drop the `FROM_BE` calls in the
   `handle_*_request` family. Python side already constructs
   these payloads with `!IIIIQ256s`-style packs that the §4.7(a)
   batches will touch. **[done — `ee832ba93`]**
3-5. **Per-command batches 1-3 (combined)** — single mechanical
   commit covering CMD_GET_CONFIG / CMD_GET_ACTIVE_VPORT /
   CMD_GET_DIMENSIONS / is_image_loaded / is_sequence_loaded /
   selections / rectangles / regions / CMD_GET_STAR_IN_SELECTION /
   CMD_GET_STATS_FOR_SELECTION / CMD_PLOT / polygon IDs /
   `gui/user_polygons.c` SHM payload / plot.py + models.py
   inline-payload helpers / remaining stragglers. **[done — `f2d5b791a`]**
   (Originally planned as 3 separate commits but the changes are
   mechanical and symmetric enough that one diff is reviewable,
   the smoke test exercises commands across all 3 categories, and
   the split was costing more in churn than it bought in
   reviewability.)

6. **Bulk-struct: `keywords_to_py` ↔ `FKeywords`** — the
   largest single field-width audit. Rewrite the C-side casts to
   native C types, rewrite `_KEYWORD_FORMAT_PARTS` to match. Each
   field is verified individually. Round-trip test via
   `SirilInterface.get_image_keywords()` against a known image.
7. **Bulk-struct: `imstats_to_py` ↔ `ImageStats`**.
8. **Bulk-struct: `homography_to_py` ↔ `Homography`** —
   in this commit ALSO add the missing `Homography.serialize` /
   `.deserialize` (§4.7bis) so the field-format string lives on
   the class. Retire any inline unpacker.
9. **Bulk-struct: `distodata_to_py` + `imgdata_to_py` ↔
    `DistoData`, `ImgData`** — same: add `DistoData.serialize` /
    `.deserialize` (§4.7bis), retire inline parsers.
10. **Bulk-struct: `psfstar_to_py` + `analysis_to_py` ↔
    `PSFStar`, analysis tuple**.
11. **Bulk-struct: `fits_to_py` ↔ `FFit`** — the composite that
    embeds several of the above. Verify the per-block offsets
    still add up. Add the missing `FFit.deserialize_core(...)`
    classmethod (§4.7bis), replace both `get_image()` and
    `get_seq_frame()` to call it instead of the duplicated
    inline `format_parts` list. Largest integration test.

11a. **Class consolidation cleanup** — add the remaining missing
    methods that aren't naturally landed in the per-block
    commits above: `Sequence.serialize` / `.deserialize`,
    `FPoint.deserialize` (complement existing serialize),
    `_SharedMemoryInfo.from_buffer` (retire the 3× copy-pasted
    unpack pattern in `get_seq_frame`). General sweep of
    `connection.py` `struct.unpack` sites — anything with 2+
    fields gets moved onto its model class.

12. **Final sweep + version bump** — `grep` confirms zero
    remaining `!` format chars, zero `BE64`/`GUINT*_*_BE`
    references, and zero inline 2+ field `struct.unpack` calls
    in `connection.py` for which a model class exists. Bump
    sirilpy to 1.1.22; update the uv conversion plan's
    compatibility note.

The single biggest source of risk is steps 6–11 (the bulk-struct
audits). If any of those gets the field types wrong, the
deserialised values silently come out as nonsense from that field
onward. The `_SIZE` self-check in §4.9 catches block-size
mismatches, but it can't detect a uint32 read as uint16 (the
block size happens to match by coincidence). Field-by-field
round-trip values against known inputs is the only reliable
check; budget time for it.

## 5. Risks and mitigations

1. **Mismatched protocol if a stale sirilpy is loaded** — e.g. user
   has Siril 1.5 → 1.5.x upgrade, and a per-script venv's sirilpy
   was 1.1.21 but the new Siril expects 1.1.22 wire format.

   **Mitigation**: the per-script-venv freshness check (uv
   conversion §4.5, Option 3) already triggers a surgical sirilpy
   reinstall on every script launch when the recorded
   `sirilpy_version` differs from `com.python_version`. So upgrading
   Siril → all per-script venvs get the new sirilpy on first use.
   Base venv: handled at startup by `uv_install_sirilpy`.
   **No additional work needed** — the existing infrastructure
   already enforces this lockstep.

2. **Missed a conversion site** — silent data corruption (garbage
   integers) on whatever endpoint is missed.

   **Mitigation**: the per-batch grep checks in §4.9 catch the C
   side; the Python `!` grep catches Python. Plus the functional
   smoke tests touch every category. After each commit, do both
   greps + run the smoke test.

3. **`struct.pack('=...')` size mismatch with `__attribute__((packed))`
   on the C side** — would corrupt headers.

   **Mitigation**: `=` in Python `struct` means *standard sizes, no
   alignment*, which is exactly what the C side does (memcpy of
   each field followed by `ptr += sizeof(field)`). Sizes are
   identical. Verified by the very first commit's smoke test (any
   pyscript that runs at all uses the header path).

4. **Test gap: macOS aarch64** — also little-endian, but a
   different ABI from x86_64.

   **Mitigation**: `=` and `__attribute__((packed))` are about byte
   layout, not ABI. They don't differ across LE architectures. If
   the change works on x86_64 Linux it works on aarch64 macOS.
   Manual smoke test by the macOS maintainer when they pick this
   up post-merge.

5. **A new contributor adds a `GUINT32_TO_BE` later** — the
   inconsistency creeps back in.

   **Mitigation**: the C-side grep returning zero is a one-line
   CI check we could add. Out of scope for this PR but worth
   noting as a follow-up.

## 6. Summary of code touch-points

| Done  | File                                            | Change                                                |
| ----- | ----------------------------------------------- | ----------------------------------------------------- |
| `[ ]` | `src/io/siril_pythoncommands.c`                 | rename `COPY_BE64` → `COPY_FIELD` (keep simplified, ~97 call sites; each one re-typed to native C width per §4.2bis); delete `FROM_BE64_INTO` (~13 sites → memcpy); delete `TO_BE64_INTO` (~17 sites → `dest = val`); drop ~76 remaining inline `GUINTxx_*_BE` calls |
| `[ ]` | `src/io/siril_pythonmodule.c`                   | drop ~16 inline `GUINTxx_*_BE` calls (CommandHeader, ResponseHeader, incoming_image_info_t) |
| `[ ]` | `python_module/sirilpy/connection.py`           | swap all 68 `struct.(pack\|unpack)('!…')` → `('=…')` (inline payloads — width unchanged) |
| `[ ]` | `python_module/sirilpy/models.py`               | rewrite `KEYWORDS_FORMAT` (and analogous `_FORMAT` ClassVars for ImageStats / Homography / DistoData / ImgData / PSFStar / FFit) to use native C widths; `=` prefix; per-field char audit per §4.2bis. **Also**: add the missing serialize/deserialize methods on `FFit`, `Homography`, `DistoData`, `Sequence`, `FPoint` (§4.7bis), retiring the inline parsers in `connection.py` |
| `[ ]` | `python_module/sirilpy/shm.py`                  | add `_SharedMemoryInfo.from_buffer(buf, offset=0)` classmethod to retire the 3–4× copy-pasted unpack pattern in `connection.py:get_seq_frame` (§4.7bis) |
| `[ ]` | `python_module/pyproject.toml`                  | version bump 1.1.21 → 1.1.22 |
| `[—]` | `src/io/siril_pythonmodule.h`                   | (no change — struct layouts unchanged) |

End of plan.
