#!/bin/bash
set -e

# Install Miniconda if not already installed
if [ ! -d "$HOME/miniconda" ]; then
  echo "Installing Miniconda..."
  wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
  bash miniconda.sh -b -p $HOME/miniconda
  rm miniconda.sh
fi

# Add conda to PATH
export PATH="$HOME/miniconda/bin:$PATH"
source "$HOME/miniconda/etc/profile.d/conda.sh"

# Create/update env
ENV_NAME="admet_env"
if conda env list | grep -q "$ENV_NAME"; then
  echo "Updating environment..."
  conda env update -f environment.yml --prune
else
  echo "Creating environment..."
  conda env create -f environment.yml
fi

conda activate $ENV_NAME

# Environment variables
export FLASK_APP=app.py
export FLASK_ENV=production
export PORT=${PORT:-5000}

# Start backend with Gunicorn
echo "Starting Gunicorn on port $PORT..."
exec gunicorn -w 4 -b 0.0.0.0:$PORT app:app
