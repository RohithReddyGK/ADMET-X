#!/bin/bash
set -e  # Exit immediately if a command fails

# ---------------- Install Miniconda if not already installed ----------------
if [ ! -d "$HOME/miniconda" ]; then
  echo "Installing Miniconda..."
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
  bash miniconda.sh -b -p $HOME/miniconda
  rm miniconda.sh
fi

export PATH="$HOME/miniconda/bin:$PATH"
source "$HOME/miniconda/etc/profile.d/conda.sh"

# ---------------- Create / Update Conda environment ----------------
ENV_NAME="admet_env"

if conda env list | grep -q "$ENV_NAME"; then
  echo "Updating existing Conda environment '$ENV_NAME'..."
  conda env update -f environment.yml --prune
else
  echo "Creating Conda environment '$ENV_NAME'..."
  conda env create -f environment.yml
fi

# Activate environment
echo "Activating Conda environment '$ENV_NAME'..."
conda activate $ENV_NAME

# ---------------- Set environment variables ----------------
export FLASK_APP=app.py      # Main backend file
export FLASK_ENV=production
export PORT=${PORT:-5000}    # Use Render-assigned PORT or fallback to 5000

# ---------------- Start Flask server ----------------
echo "Starting Flask server on port $PORT..."
flask run --host=0.0.0.0 --port=$PORT
