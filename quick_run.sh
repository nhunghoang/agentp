#!/bin/bash
set -euo pipefail

wget --no-clobber https://github.com/nhunghoang/agentp/blob/main/install.sh
wget --no-clobber https://github.com/nhunghoang/agentp/blob/main/run_example.sh

chmod +x install.sh
chmod +x run_example.sh

sudo bash install.sh
bash run_example.sh