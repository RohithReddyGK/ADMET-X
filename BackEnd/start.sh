#!/bin/bash
set -e

# ---------------- Install Miniconda ----------------
if [ ! -d "$HOME/miniconda" ]; then
  echo "Installing Miniconda..."
  wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
  bash miniconda.sh -b -p $HOME/miniconda
  rm miniconda.sh
fi

export PATH="$HOME/miniconda/bin:$PATH"
source "$HOME/miniconda/etc/profile.d/conda.sh"

# ---------------- Create / Update Environment ----------------
ENV_NAME="admet_env"
if conda env list | grep -q "$ENV_NAME"; then
  echo "Updating environment..."
  conda env update -f environment.yml --prune
else
  echo "Creating environment..."
  conda env create -f environment.yml
fi

# Activate environment
conda activate "$ENV_NAME"

# ---------------- Set Environment Variables ----------------
export FLASK_APP=app.py
export FLASK_ENV=production
export PORT=${PORT:-5000}

# ---------------- Start Gunicorn ----------------
echo "Starting Gunicorn on port $PORT..."
exec gunicorn --bind 0.0.0.0:$PORT app:app
