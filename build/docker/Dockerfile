# syntax=docker/dockerfile:1
###############################################################################
# Siril headless CLI image
#
# Produces a small Debian-based image containing ONLY the command-line
# `siril-cli` binary, built with meson `-Dgtk=false`: no GTK, no GtkSourceView
# and no X11 are linked or installed.  All optional dependencies are enabled
# (libraw, TIFF/JPEG/PNG/HEIF/JXL, exiv2, ffms2/ffmpeg, libcurl, libgit2,
# sqlite, XISF, healpix/htmesh, …); only the test framework (criterion) is left
# out.
#
# Build (from the repository root, with submodules checked out):
#   git submodule update --init --recursive
#   docker build -f build/docker/Dockerfile -t siril-cli .
#
# Run:
#   docker run --rm -v "$PWD:/data" --user "$(id -u):$(id -g)" \
#       siril-cli -s /data/script.ssf
#
# Formal release procedure (publishing to Docker Hub)
# ---------------------------------------------------
# CI never pushes to a registry; release images are published by hand from the
# artifact of the release-tag pipeline.  To cut a release:
#
#   1. Push the release tag (e.g. 1.4.4).  This runs the `docker-cli-image` CI
#      job, which builds + smoke-tests the image and uploads it as the artifact
#      `siril-cli-image.tar.gz` (the image inside is tagged `siril-cli:<tag>`).
#   2. Download that artifact from the tagged pipeline, then load, retag and
#      push it (replace `siril` with the actual Docker Hub namespace):
#
#        VER=1.4.4                       # the release tag
#        REPO=siril/siril                # Docker Hub <namespace>/<repo>
#        docker load -i siril-cli-image.tar.gz          # -> siril-cli:$VER
#        docker tag siril-cli:$VER $REPO:$VER
#        docker tag siril-cli:$VER $REPO:${VER%.*}       # e.g. 1.4
#        docker tag siril-cli:$VER $REPO:latest
#        docker login -u <dockerhub-user>
#        docker push $REPO:$VER
#        docker push $REPO:${VER%.*}
#        docker push $REPO:latest
#
#   Only move `latest` for a stable release, not for pre-releases/RCs.
###############################################################################

###############################################################################
# Stage 1 — builder
###############################################################################
FROM debian:testing AS build

ENV DEBIAN_FRONTEND=noninteractive

# Build toolchain + all (non-GUI) build dependencies.  Deliberately omits
# libgtk-4-dev / libgtksourceview-5-dev (GUI) and libcriterion-dev (tests).
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake ninja-build pkg-config git ca-certificates \
        python3 python3-pip \
        libglib2.0-dev libjson-glib-dev \
        libcairo2-dev libpango1.0-dev libgdk-pixbuf-2.0-dev \
        libgsl-dev libfftw3-dev liblcms2-dev \
        wcslib-dev libcfitsio-dev \
        libopencv-dev \
        libraw-dev libtiff-dev libjpeg-dev libpng-dev libheif-dev libjxl-dev \
        libexiv2-dev \
        libffms2-dev libavformat-dev libavcodec-dev libavutil-dev \
        libswscale-dev libswresample-dev \
        libcurl4-gnutls-dev libgit2-dev libsqlite3-dev \
        libzstd-dev liblz4-dev libpugixml-dev libxisf-dev \
        libhealpix-cxx-dev \
    && rm -rf /var/lib/apt/lists/*

# Match the version pins the project uses in CI (avoid known-bad meson/ninja).
RUN pip install --break-system-packages "meson!=1.10.0" "ninja!=1.11.1"

WORKDIR /src
COPY . /src

# Configure a headless, release build with every optional feature explicitly
# enabled.  Passing each option explicitly makes the build fail loudly if an
# optional dependency is missing, rather than silently dropping the feature.
RUN meson setup _build \
        --prefix=/usr \
        --buildtype=release \
        -Dgtk=false \
        -Dcriterion=false \
        -Dopenmp=true \
        -Dlibcurl=true -Dsqlite=true -Dexiv2=true -Dlibgit2=true \
        -Dlibraw=true -Dlibtiff=true -Dlibjpeg=true -Dlibpng=true \
        -Dlibheif=true -Dlibjxl=true -DlibXISF=true \
        -Dffms2=true -Dffmpeg=true \
    && ninja -C _build

# Install into a staging tree, then compute the exact set of Debian runtime
# packages the installed binaries/libraries actually link against.  This keeps
# the runtime image minimal and version-correct without hard-coding fragile
# `libfoo<sover>` package names.  Criterion is never linked into siril-cli, so
# it is naturally excluded.
RUN set -eux; \
    DESTDIR=/install ninja -C _build install; \
    mkdir -p /install-meta; \
    find /install -type f \( -perm -u+x -o -name '*.so*' \) -print0 \
      | xargs -0 -r ldd 2>/dev/null \
      | awk '/=> \// && $3 ~ /^\/(usr\/)?lib/ { print $3 }' \
      | sort -u \
      | xargs -r dpkg -S 2>/dev/null \
      | cut -d: -f1 | tr ',' '\n' | sort -u \
      > /install-meta/run-deps.txt; \
    echo "Runtime packages:"; cat /install-meta/run-deps.txt

###############################################################################
# Stage 2 — runtime
###############################################################################
FROM debian:testing-slim AS runtime

ENV DEBIAN_FRONTEND=noninteractive \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

# Bring in the built tree (prefix /usr) and the computed runtime-dep list.
COPY --from=build /install /
COPY --from=build /install-meta/run-deps.txt /tmp/run-deps.txt

# Install the auto-detected runtime libraries plus a few things that are
# loaded/spawned rather than linked, so they don't show up in ldd:
#   ca-certificates + glib-networking — TLS for libcurl / gio (online catalogues)
#   python3 + python3-venv           — Siril's Python scripting (pyscript)
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates glib-networking \
        python3 python3-venv \
        $(cat /tmp/run-deps.txt) \
    && rm -rf /var/lib/apt/lists/* /tmp/run-deps.txt \
    && ldconfig

# Run as a non-root user with uid/gid 1000 so files written to a bind-mounted
# /data are owned by the typical host user.  Override with `--user $(id -u):$(id -g)`
# when your host uid differs.
#
# Siril keeps all its state under $HOME via the XDG dirs (~/.config/siril,
# ~/.local/share/siril, ~/.local/share/siril-scripts, the SPCC database, …).
# Pin HOME and make it world-writable so the image still works when run with
# `--user <arbitrary-uid>`: such a uid has no /etc/passwd entry, so Docker
# would otherwise leave HOME=/ and Siril could not create ~/.config.  The
# cache is redirected to /tmp so it never bloats a persisted home volume.
# Mount a volume at /home/siril (or set XDG_*_HOME) to persist state across runs.
ENV HOME=/home/siril \
    XDG_CACHE_HOME=/tmp/siril-cache
RUN groupadd -g 1000 siril && useradd -u 1000 -g 1000 -m -s /bin/bash siril \
    && chmod 0777 /home/siril
USER siril
WORKDIR /data

# `docker run <image> -s script.ssf` runs the headless CLI directly.
ENTRYPOINT ["siril-cli"]
CMD ["--version"]
