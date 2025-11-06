# powerNMA Docker Image - Production-Ready Network Meta-Analysis
# Ultimate Edition with all Phase 1-11 features

FROM rocker/shiny:4.3.0

# Metadata
LABEL maintainer="powerNMA Development Team"
LABEL description="Ultimate Network Meta-Analysis Platform with ML, AI, and Cloud-Ready Architecture"
LABEL version="1.0.0-ultimate"

# Set working directory
WORKDIR /srv/shiny-server/powernma

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgit2-dev \
    libsqlite3-dev \
    libpq-dev \
    python3-pip \
    python3-dev \
    git \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Install Ollama (for AI features)
RUN wget https://ollama.ai/install.sh -O install_ollama.sh \
    && chmod +x install_ollama.sh \
    && ./install_ollama.sh \
    && rm install_ollama.sh

# Install R package dependencies
RUN R -e "install.packages(c( \
    'shiny', 'shinydashboard', 'shinyWidgets', 'shinyjs', 'shinyBS', \
    'DT', 'plotly', 'ggplot2', 'tidyverse', \
    'netmeta', 'meta', 'metafor', 'gemtc', 'multinma', \
    'igraph', 'network', 'visNetwork', \
    'rmarkdown', 'knitr', 'officer', 'officedown', \
    'jsonlite', 'httr', 'plumber', \
    'RSQLite', 'DBI', 'pool', \
    'future', 'promises', 'parallel', 'foreach', 'doParallel', \
    'randomForest', 'gbm', 'xgboost', 'nnet', 'caret', \
    'tm', 'e1071', 'text2vec', \
    'devtools', 'roxygen2', 'testthat', \
    'readxl', 'writexl', 'openxlsx', \
    'htmlwidgets', 'htmltools', \
    'lubridate', 'stringr', 'purrr', 'forcats' \
    ), repos='https://cloud.r-project.org/')"

# Copy powerNMA package
COPY . /tmp/powerNMA

# Install powerNMA package
RUN R -e "devtools::install('/tmp/powerNMA', dependencies = TRUE)"

# Create necessary directories
RUN mkdir -p /srv/shiny-server/powernma/data \
    && mkdir -p /srv/shiny-server/powernma/results \
    && mkdir -p /srv/shiny-server/powernma/logs \
    && mkdir -p /srv/shiny-server/powernma/cache \
    && mkdir -p /srv/shiny-server/powernma/uploads \
    && chown -R shiny:shiny /srv/shiny-server/powernma

# Copy Shiny app files
RUN R -e "library(powerNMA); \
    cat('library(powerNMA)\nlaunch_powernma_app(host=\"0.0.0.0\", port=3838)', \
    file='/srv/shiny-server/powernma/app.R')"

# Copy configuration files
COPY inst/docker/shiny-server.conf /etc/shiny-server/shiny-server.conf
COPY inst/docker/startup.sh /usr/local/bin/startup.sh
RUN chmod +x /usr/local/bin/startup.sh

# Expose ports
EXPOSE 3838 8000 11434

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD curl -f http://localhost:3838/ || exit 1

# Set environment variables
ENV SHINY_LOG_STDERR=1 \
    SHINY_LOG_LEVEL=INFO \
    POWERNMA_DATA_DIR=/srv/shiny-server/powernma/data \
    POWERNMA_RESULTS_DIR=/srv/shiny-server/powernma/results \
    POWERNMA_MAX_UPLOAD_MB=500 \
    POWERNMA_N_WORKERS=4

# Volume for persistent data
VOLUME ["/srv/shiny-server/powernma/data", "/srv/shiny-server/powernma/results"]

# Run startup script
CMD ["/usr/local/bin/startup.sh"]
