#!/bin/bash
set -e

# Create /tools directory if not exists
mkdir -p /tools
cd /tools

echo "Installing bioinformatics tools..."


# Install SRA Toolkit
SRA_VERSION="3.4.1"

if [ ! -d "/tools/sratoolkit" ]; then
		echo "Installing SRA Toolkit ${SRA_VERSION}..."
		wget -q --show-progress \
				--output-document /tools/sratoolkit.tar.gz \
				"https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/sratoolkit.${SRA_VERSION}-ubuntu64.tar.gz"
		tar -xzf /tools/sratoolkit.tar.gz -C /tools/
		# Rename versioned dir to stable name (e.g. sratoolkit.3.4.1-ubuntu64 → sratoolkit)
		mv /tools/sratoolkit.${SRA_VERSION}-ubuntu64 /tools/sratoolkit
		rm -f /tools/sratoolkit.tar.gz

		# Configure cache directory non-interactively (no vdb-config -i needed)
		mkdir -p /root/.ncbi /data/sra_cache
		cat > /root/.ncbi/user-settings.mkfg << EOF
/repository/user/main/public/enabled = "true"
/repository/user/main/public/root = "/data/sra_cache"
EOF
		echo "✓ SRA Toolkit ${SRA_VERSION} installed, cache → /data/sra_cache"
else
		echo "✓ SRA Toolkit already installed"
fi


echo "All bioinformatics tools installed successfully!"

# Verify installations
echo "Verifying tool installations..."
ls -la /tools/
echo ""
/tools/sratoolkit/bin/fastq-dump --version 2>&1 | head -1
echo "Installation verification complete."