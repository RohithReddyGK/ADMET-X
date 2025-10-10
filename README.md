<p align="center">
  <h1>An AI-Driven Platform for In-Silico ADMET Prediction</h1>
</p>

<p align="center">
  <img width="200" height="200" alt="ADMET-AI Logo" src="https://github.com/user-attachments/assets/af4851f7-a076-486d-9cce-06365cfe0bb6" />
</p>

<p align="center">
  <!-- License & Core -->
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"></a>
  <a href="https://www.python.org/"><img src="https://img.shields.io/badge/Python-3.12-blue.svg" alt="Python"></a>
  <a href="https://reactjs.org/"><img src="https://img.shields.io/badge/React-18.2.0-blue.svg" alt="React"></a>
  <a href="https://vitejs.dev/"><img src="https://img.shields.io/badge/Vite-React-orange.svg" alt="Vite + React"></a>
  <a href="https://www.docker.com/"><img src="https://img.shields.io/badge/Docker-Container-blue.svg" alt="Docker"></a>
  <a href="https://fly.io/"><img src="https://img.shields.io/badge/Deployment-Fly.io-purple.svg" alt="Fly.io"></a>
</p>
<p align="center">
  <!-- Libraries/Tools -->
  <a href="https://www.rdkit.org/"><img src="https://img.shields.io/badge/RDKit-Chemistry-green.svg" alt="RDKit"></a>
  <a href="https://flask.palletsprojects.com/"><img src="https://img.shields.io/badge/Flask-Backend-orange.svg" alt="Flask"></a>
  <a href="https://framer.com/motion/"><img src="https://img.shields.io/badge/FramerMotion-Animation-pink.svg" alt="Framer Motion"></a>
  <a href="https://lottiefiles.com/"><img src="https://img.shields.io/badge/Lottie-Animations-blue.svg" alt="Lottie"></a>
</p>

**ADMET-AI** is a comprehensive, AI-powered platform for predicting **Absorption, Distribution, Metabolism, Excretion, and Toxicity (ADMET)** properties of drug molecules based on their SMILES representations. It leverages machine learning models to provide accurate drug-likeness assessments and visual insights through radar plots and molecular drawings.

---

## 🚀 Project Overview

ADMET-AI aims to accelerate early-stage drug discovery by providing in-silico predictions of pharmacokinetics and pharmacodynamics, reducing the need for extensive laboratory experiments.  

Key highlights:

- Batch and single SMILES prediction.
- Downlodable report via CSV.
- Interactive radar plots for ADMET visualization.
- Toxicity prediction integrated with traditional ADME properties.
- Support for drawing molecules and uploading files.
- Ready-to-deploy with Docker and Fly.io

---

## 💡 Features

- **SMILES Input Options**: Text input, file upload (.txt/.csv), drawing molecules, or example molecules.
- **ADMET Predictions**:
  - Absorption: Bioavailability, Caco2, HIA, etc.
  - Distribution: BBB, PPBR, etc.
  - Metabolism: CYP450 enzyme, etc.
  - Excretion: Clearance, Half-Life, etc.
  - Toxicity: AMES, Carcinogenicity, hERG, LD50, and more.
- **Interactive Visualization**: Molecule images, radar plots, and color-coded property status.
- **Export Results**: Download predictions as a CSV file.
- **Local & Online Deployment**: Works both locally and via Fly.io deployment.

---

## 🛠 Tech Stack

- **Backend**: Python, Flask, Joblib, Flask-CORS
- **Frontend**: Vite-React, TailwindCSS, Framer-Motion, Lottie
- **AI/ML Models**: Custom trained models for ADMET prediction
- **Deployment**: Docker, Fly.io and Vercel
- **Utilities**: RDKit, ChemUtils, Plotting modules

---

## ⚙️ Installation (Local Setup)

### **1. Clone the repository**
```bash
git clone https://github.com/yourusername/ADMET.git
cd ADMET/BackEnd
```

### **2. Setup Conda Environment**
```bash
conda env create -f environment.yml
conda activate admet_env
```

### **3. Run BackEnd**
```bash
python app.py
```
### **Running localhosts**
```bash
Backend will run at: http://127.0.0.1:8080
Test endpoint: http://127.0.0.1:8080/ → should return "Backend is running."
```

### **4. Setup FrontEnd**
```bash
cd ../FrontEnd
npm install
npm run dev
```

---

## 🌐 Deployment

### **1. Docker**
```bash
docker build -t admet-ai .
docker run -p 8080:8080 admet-ai
```

### **2. Fly.io**
```bash
"C:\Users\Rohith Reddy G K\.fly\bin\flyctl.exe" launch
```

---

## 🖥 Usage

- Input molecule SMILES via text, file, draw, or example.
- Click Predict.
- View interactive ADMET radar plots and molecular images.
- Optionally, download all results as CSV.

## 🧪 Contribution

Contributions are welcome! To contribute:
-Fork the repository
-Create a branch: git checkout -b feature-name
-Make changes and commit: git commit -m "Add new feature"
-Push to branch: git push origin feature-name
-Open a Pull Request



