# ---------- Use Miniconda base ----------
FROM continuumio/miniconda3

# Set working directory
WORKDIR /app

# ---------- Copy environment file first ----------
COPY environment.yml .
<<<<<<< HEAD
RUN conda env create -f environment.yml

# Copy the rest of the project
COPY . .

# Expose port 8080
EXPOSE 8080

# Run Flask app in the conda environment
CMD ["conda", "run", "-n", "admet_env", "python", "app.py"]
=======

# Create conda environment
RUN conda env create -f environment.yml && conda clean -afy

# Activate conda environment by default
SHELL ["conda", "run", "-n", "admet_env", "/bin/bash", "-c"]

# ---------- Copy project files ----------
COPY app.py .
COPY utils/ ./utils/
COPY descriptors/ ./descriptors/
COPY Models/ ./Models/

# Expose port
EXPOSE 8080

# ---------- Run the app ----------
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "admet_env", "gunicorn", "--bind", "0.0.0.0:8080", "app:app", "--workers", "1", "--threads", "4"]
>>>>>>> af9ce8401c07dbbf694915af6756665fd57727a9
