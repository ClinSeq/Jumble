FROM rocker/r-ver:4.4.2

LABEL org.opencontainers.image.title="Jumble" \
      org.opencontainers.image.description="Copy number analysis of short read sequencing data" \
      org.opencontainers.image.source="https://github.com/ClinSeq/Jumble" \
      org.opencontainers.image.licenses="GPL-3.0-or-later"

# System libraries for Rsamtools/bamsignals (htslib) and curl/xml based Bioconductor deps
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/Jumble

# 1. Install CRAN/Bioconductor dependencies from DESCRIPTION only,
#    so this expensive layer stays cached when only R code changes.
COPY DESCRIPTION .
RUN Rscript -e 'install.packages("pak", repos = sprintf("https://r-lib.github.io/p/pak/stable/%s/%s/%s", .Platform$pkgType, R.Version()$os, R.Version()$arch))' \
    && Rscript -e 'pak::local_install_deps(dependencies = "hard")' \
    && rm -rf /tmp/* /root/.cache

# 2. Install the Jumble package itself
COPY . .
RUN R CMD INSTALL --no-multiarch --without-keep.source . \
    && Rscript -e 'library(Jumble)' \
    && chmod +x /usr/local/lib/R/site-library/Jumble/scripts/*.R

# Command-line wrappers (jumble-run.R, jumble-count.R, ...) on PATH
ENV PATH="/usr/local/lib/R/site-library/Jumble/scripts:${PATH}"


CMD ["R"]
