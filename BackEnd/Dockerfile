# Use Miniconda base image
FROM continuumio/miniconda3

# Set working directory
WORKDIR /app

# Copy environment.yml and create environment
COPY environment.yml .

RUN conda env create -f environment.yml

# Activate the environment automatically
SHELL ["conda", "run", "-n", "admet_env", "/bin/bash", "-c"]

# Copy the rest of the project
COPY . .

# Expose port 8080 (Fly.io expects this)
EXPOSE 8080

# Run the Flask app
CMD ["conda", "run", "-n", "admet_env", "python", "app.py"]
