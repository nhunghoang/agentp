# AGENT-P runtime image
#
# Build from the repo root (plink2 and regenie are x86_64-only, so on
# Apple Silicon the --platform flag is required):
#   docker build --platform linux/amd64 -t agentp-img .
#
# Run the bundled example:
#   docker run --platform linux/amd64 -it agentp-img
#   cd /usr/home/example && python run_example.py

FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# System deps, R, and the R packages FUSION needs (from Ubuntu so nothing
# is compiled from CRAN). libgomp1 is required by the regenie binary.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bash \
        build-essential \
        ca-certificates \
        curl \
        g++-10 \
        gcc-10 \
        libgomp1 \
        python3 \
        python3-dev \
        python3-pip \
        r-base-core \
        r-base-dev \
        r-cran-optparse \
        r-cran-rcolorbrewer \
        r-cran-rcpp \
        r-cran-rcppeigen \
        unzip \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# agentp shells out with brace expansion ({2..22}), which needs bash as sh.
# It also calls `python`, not `python3`.
RUN ln -sf /bin/bash /bin/sh \
    && ln -sf /usr/bin/python3 /usr/bin/python

# bgenix, built from a pinned bgen release
ARG BGEN_VERSION=v1.1.7
RUN curl -fsSL "http://code.enkre.net/bgen/tarball/release/${BGEN_VERSION}" -o /tmp/bgen.tgz \
    && tar -xzf /tmp/bgen.tgz -C /opt \
    && rm /tmp/bgen.tgz \
    && cd "/opt/${BGEN_VERSION}" \
    && CC=gcc-10 CXX=g++-10 python3 waf configure \
    && python3 waf \
    && ln -s "/opt/${BGEN_VERSION}/build/apps/bgenix" /usr/local/bin/bgenix

# regenie, pinned prebuilt release binary
ARG REGENIE_VERSION=4.1.3
RUN curl -fsSL "https://github.com/rgcgithub/regenie/releases/download/v${REGENIE_VERSION}/regenie_v${REGENIE_VERSION}.gz_x86_64_Linux.zip" -o /tmp/regenie.zip \
    && unzip -q /tmp/regenie.zip -d /tmp/regenie \
    && install -m 755 "/tmp/regenie/regenie_v${REGENIE_VERSION}.gz_x86_64_Linux" /usr/local/bin/regenie \
    && rm -rf /tmp/regenie /tmp/regenie.zip

# plink2
RUN curl -fsSL https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20241203.zip -o /tmp/plink2.zip \
    && unzip -q /tmp/plink2.zip -d /usr/local/bin \
    && chmod +x /usr/local/bin/plink2 \
    && rm /tmp/plink2.zip

# Install agentp from this checkout (not PyPI). Editable so the bundled
# external/ scripts and models/ resolve to /opt/agentp/agentp/...
# setuptools<81 keeps pkg_resources, which agentp/Static.py imports.
WORKDIR /opt/agentp
COPY setup.py MANIFEST.in README.md LICENSE ./
COPY agentp ./agentp
RUN python3 -m pip install --no-cache-dir --upgrade pip "setuptools<81" wheel \
    && python3 -m pip install --no-cache-dir -e .

# plink2R (required by FUSION) from the bundled source
RUN R CMD INSTALL agentp/external/plink2R/plink2R

# Fail the build if the models are Git LFS pointers instead of real files,
# and warn about any model set that is empty.
RUN if grep -rlsI --include='*.db' --include='*.gz' --include='*.RDat' \
            '^version https://git-lfs' agentp/models | head -n 5 | grep -q .; then \
        echo 'ERROR: agentp/models contains Git LFS pointers. Run `git lfs pull` before building.' >&2; \
        exit 1; \
    fi \
    && for d in JTI/weights JTI/covariance PrediXcan/weights PrediXcan/covariance \
                UTMOST_dzhou/weights UTMOST_dzhou/covariance FUSION/WEIGHTS FUSION/LDREF; do \
        [ -n "$(ls -A "agentp/models/$d" 2>/dev/null)" ] \
            || echo "WARNING: agentp/models/$d is empty; that model set will not work." >&2; \
    done

# Example data
COPY example /usr/home/example

# Smoke test: every tool agentp shells out to must be callable
RUN python -c "import agentp" \
    && plink2 --version \
    && ! (ldd /usr/local/bin/regenie | grep 'not found') \
    && regenie --version \
    && (bgenix -help 2>&1 | head -n 1) \
    && Rscript -e "suppressMessages({library(plink2R); library(optparse)})"

WORKDIR /usr/home
CMD ["/bin/bash"]
