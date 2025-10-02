#!/bin/bash
set -e  # Exit immediately if a command fails

# Activate Conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate admet_env

# Set environment variables
export FLASK_APP=app.py
export FLASK_ENV=production
export PORT=${PORT:-5000}

echo "Starting Gunicorn on port $PORT..."
exec gunicorn -w 4 -b 0.0.0.0:$PORT app:app
