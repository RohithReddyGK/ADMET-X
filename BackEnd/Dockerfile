# Use Miniconda base image
FROM continuumio/miniconda3

# Set working directory
WORKDIR /app

# Copy environment.yml and create environment
COPY environment.yml .
RUN conda env create -f environment.yml

# Copy the rest of the project
COPY . .

# Expose port 8080
EXPOSE 8080

# Run Flask app in the conda environment
CMD ["conda", "run", "-n", "admet_env", "python", "app.py"]
