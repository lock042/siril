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

## Pulling

Images are published to the project GitLab Container Registry:

```sh
docker pull registry.gitlab.com/free-astro/siril:latest   # newest release
docker pull registry.gitlab.com/free-astro/siril:1.4       # latest 1.4.x
docker pull registry.gitlab.com/free-astro/siril:1.4.4     # exact version
docker pull registry.gitlab.com/free-astro/siril:edge      # tip of master
```

## Running

The entrypoint is `siril-cli`, so arguments are passed straight through. Mount
the directory holding your data/scripts at `/data` (the working directory):

```sh
# Run a script over data in the current directory
docker run --rm -v "$PWD:/data" --user "$(id -u):$(id -g)" \
    registry.gitlab.com/free-astro/siril -s /data/process.ssf

# Quick check that it works headlessly
docker run --rm registry.gitlab.com/free-astro/siril --version
```

### File ownership

Files written by the container are owned by its user. The image defaults to
uid/gid `1000`; if your host user differs, add `--user "$(id -u):$(id -g)"`
(as above) so outputs land owned by you rather than by root or uid 1000.

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
