# Siril headless CLI Docker image

A small Debian-based image containing the command-line `siril-cli` only. It is
built with meson `-Dgtk=false`, so **no GTK, GtkSourceView or X11** is linked or
installed — it is purpose-built for headless/batch/server use (CI pipelines,
automated processing, projects like astrolab). All optional dependencies are
enabled (libraw, TIFF/JPEG/PNG/HEIF/JXL, exiv2, ffms2/ffmpeg, libcurl, libgit2,
sqlite, XISF, healpix/htmesh, Python scripting); only the test framework
(criterion) is left out.

> This image is **not** a desktop/GUI delivery mechanism. For interactive use,
> install the native packages, Flatpak or AppImage. Docker is unlike Flatpak
> here: there is no sandbox interposed between Siril and your files — the
> container simply sees whatever directories you bind-mount, with normal POSIX
> access inside them.

## Obtaining the image

The CI does **not** push to a registry. The `docker-cli-image` job builds the
image, smoke-tests it, and saves it as a downloadable artifact,
`siril-cli-image.tar.gz`. Download that artifact from the pipeline and load it:

```sh
docker load -i siril-cli-image.tar.gz      # -> siril-cli:<ref>  (e.g. siril-cli:1.4.4)
docker image ls siril-cli
```

Official release images are published to Docker Hub manually from the
release-tag artifact (`docker load`, `docker tag`, `docker push`).

## Running

The entrypoint is `siril-cli`, so arguments are passed straight through. Mount
the directory holding your data/scripts at `/data` (the working directory).
Replace `siril-cli:<ref>` below with the tag you loaded:

```sh
# Run a script over data in the current directory
docker run --rm -v "$PWD:/data" --user "$(id -u):$(id -g)" \
    siril-cli:<ref> -s /data/process.ssf

# Quick check that it works headlessly
docker run --rm siril-cli:<ref> --version
```

### File ownership

Files written by the container are owned by its user. The image defaults to
uid/gid `1000`; if your host user differs, add `--user "$(id -u):$(id -g)"`
(as above) so outputs land owned by you rather than by root or uid 1000.

### Configuration & persistent state

Siril stores all of its state under `$HOME` using the standard XDG locations:

| Path                                  | Contents                                  |
| ------------------------------------- | ----------------------------------------- |
| `~/.config/siril/`                    | settings (`siril.config`), equipment, etc.|
| `~/.local/share/siril/`               | application data                          |
| `~/.local/share/siril-scripts/`       | scripts pulled from the script repository |
| `~/.local/share/siril-spcc-database/` | the SPCC photometric database             |

Inside the container `$HOME` is `/home/siril`. **By default this lives in the
container's writable layer, so with `--rm` all settings, downloaded catalogues,
the SPCC database and scripts are discarded when the container exits** — only
the data in your bind-mounted directory persists. That is usually fine for
one-shot batch processing, but to keep state between runs, persist `$HOME`:

```sh
# Named volume — simplest; survives across runs, isolated from the host tree
docker volume create siril-state
docker run --rm \
    -v siril-state:/home/siril \
    -v "$PWD:/data" \
    siril-cli:<ref> -s /data/process.ssf
```

Or fold the config into your data tree by pointing the XDG dirs at the mount
(handy if you want the settings to live alongside the project):

```sh
docker run --rm -v "$PWD:/data" \
    -e XDG_CONFIG_HOME=/data/.siril-config \
    -e XDG_DATA_HOME=/data/.siril-data \
    --user "$(id -u):$(id -g)" \
    siril-cli:<ref> -s /data/process.ssf
```

The image pins `HOME=/home/siril` and makes it world-writable, so it works even
with an arbitrary `--user <uid>` that has no entry in the container's
`/etc/passwd` (without this, Docker leaves `HOME=/` and Siril cannot create
`~/.config`). The cache (`XDG_CACHE_HOME`) is redirected to `/tmp` so it never
accumulates in a persisted home volume.

## Building locally

The build needs the meson subprojects, so check out submodules first, then
build from the repository root:

```sh
git submodule update --init --recursive
docker build -f build/docker/Dockerfile -t siril-cli .
docker run --rm siril-cli --version
```

## How it stays small

A multi-stage build compiles Siril in a full builder image, installs it into a
staging tree, then resolves (via `ldd` + `dpkg -S`) the exact set of Debian
runtime library packages the binary actually links against. Only those — plus
`ca-certificates`/`glib-networking` (TLS) and `python3`/`python3-venv` (Siril
Python scripting), which are loaded/spawned rather than linked — are installed
into a fresh `debian:testing-slim` runtime stage. No `-dev` packages, no GTK
stack, no build toolchain ship in the final image.
