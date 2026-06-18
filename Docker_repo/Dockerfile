# Use Python 3.12.7 base image
FROM python:3.12.7

# Set the environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    CONDA_DIR=/opt/conda \
    PATH="/tools/sratoolkit/bin:$PATH" \
    MPLCONFIGDIR=/tmp/matplotlib \
    FONTCONFIG_PATH=/etc/fonts \
    FONTCONFIG_FILE=/etc/fonts/fonts.conf \
    FONTCONFIG_CACHEDIR=/tmp/fontconfig \
    XDG_CACHE_HOME=/tmp/.cache \
    MPLBACKEND=Agg \
    PYTHONUNBUFFERED=1

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    unzip \
    pigz \
    git \
    curl \
    bzip2 \
    libglib2.0-0 \
    libgomp1 \
    libxext6 \
    libsm6 \
    libxrender1 \
    libtiff-dev \
    cpanminus \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p $CONDA_DIR && \
    rm /tmp/miniconda.sh && \
    $CONDA_DIR/bin/conda clean --all --yes

ENV PATH="$CONDA_DIR/bin:$PATH"

# Accept conda ToS
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Install zstandard in base conda env to suppress warnings
RUN conda install -n base -c conda-forge zstandard -y

# Create /tools directory
RUN mkdir -p /tools

# Create /tools directory and SRA cache mount point
RUN mkdir -p /tools /data/sra_cache && \
    chmod 777 /data/sra_cache

# Create dedicated conda env for your Python script + tools
RUN conda create -p /tools/xpcollect_env python=3.12.7 -y && \
    conda install -p /tools/xpcollect_env -c conda-forge -c bioconda \
        busco=6.0.0 \
        kallisto \
        perl-sort-naturally \
        perl-clone \
        perl-http-message \
        perl-libwww-perl \
    -y && \
    conda clean --all --yes

# Copy and run external tool installation script
COPY install_tools.sh /app/install_tools.sh
RUN chmod +x /app/install_tools.sh && /app/install_tools.sh

# Create directories with proper permissions for non-root users
RUN mkdir -p /tmp/matplotlib /tmp/.cache /tmp/.local /tmp/fontconfig && \
    chmod -R 777 /tmp


# Create working directory
WORKDIR /app

# Copy Python dependencies and Python scripts
COPY requirements.txt .
COPY Xpression_collector_docker.py .


# Install Python dependencies inside xpcollect_env
RUN /tools/xpcollect_env/bin/pip install --no-cache-dir -r requirements.txt

# Make /app writable for non-root users (BUSCO needs to write logs)
RUN chmod 777 /app

# Add conda env to PATH (this makes busco/agat/python available globally)
# Add conda env to PATH but keep Diamond priority
ENV PATH="/tools/xpcollect_env/bin:$PATH" \
    PERL5LIB="/tools/xpcollect_env/lib/perl5/site_perl:/tools/xpcollect_env/lib/perl5"

# Verify installations
RUN /tools/xpcollect_env/bin/busco --help > /dev/null 2>&1 && echo "BUSCO installation verified" || echo "BUSCO installation failed"
RUN /tools/xpcollect_env/bin/kallisto --help > /dev/null 2>&1 && echo "Kallisto installation verified" || echo "Kallisto installation failed"

# Use bash shell by default
SHELL ["/bin/bash", "-c"]

# Default command runs your Python script using the xpcollect_env Python
ENTRYPOINT ["/tools/xpcollect_env/bin/python", "Xpression_collector_docker.py"]

