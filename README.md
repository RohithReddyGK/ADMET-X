<div align="center">

  <!-- Title -->
  <h1 style="font-size: 10em;">An AI-Driven Platform for In-Silico ADMET Prediction</h1>

  <!-- Logo -->
  <img width="200" height="200" alt="ADMET-X Logo" src="https://github.com/user-attachments/assets/5fe01a51-fec8-441e-a973-11cfdee642a1"/>


  <!-- Badges -->
  <p>
    <a href="https://admet-x.vercel.app/"><img src="https://img.shields.io/badge/Live-Demo-green" alt="LIVE DEMO"></a>
    <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"></a>
    <a href="https://www.python.org/"><img src="https://img.shields.io/badge/Python-3.12-blue.svg" alt="Python"></a>
    <a href="https://reactjs.org/"><img src="https://img.shields.io/badge/React-18.2.0-blue.svg" alt="React"></a>
    <a href="https://vitejs.dev/"><img src="https://img.shields.io/badge/Vite-React-orange.svg" alt="Vite + React"></a>
    <a href="https://www.docker.com/"><img src="https://img.shields.io/badge/Docker-Container-blue.svg" alt="Docker"></a>
    <a href="https://fly.io/"><img src="https://img.shields.io/badge/Deployment-Fly.io-purple.svg" alt="Fly.io"></a>
    <a href="https://www.rdkit.org/"><img src="https://img.shields.io/badge/RDKit-Chemistry-green.svg" alt="RDKit"></a>
    <a href="https://flask.palletsprojects.com/"><img src="https://img.shields.io/badge/Flask-Backend-orange.svg" alt="Flask"></a>
    <a href="https://framer.com/motion/"><img src="https://img.shields.io/badge/FramerMotion-Animation-pink.svg" alt="Framer Motion"></a>
    <a href="https://lottiefiles.com/"><img src="https://img.shields.io/badge/Lottie-Animations-blue.svg" alt="Lottie"></a>
  </p>

</div>

**ADMET-X** is a comprehensive, AI-powered platform for predicting **Absorption, Distribution, Metabolism, Excretion, and Toxicity (ADMET)** properties of drug molecules based on their SMILES representations. It leverages machine learning models to provide accurate drug-likeness assessments and visual insights through radar plots and molecular drawings.

---

## 🚀 Project Overview

ADMET-X aims to accelerate early-stage drug discovery by providing in-silico predictions of pharmacokinetics and pharmacodynamics, reducing the need for extensive laboratory experiments.  

Key highlights:

- Batch and single SMILES prediction.
- Downlodable report via CSV.
- Interactive radar plots for ADMET visualization.
- Toxicity prediction integrated with traditional ADME properties.
- Support for drawing molecules and uploading files.
- Ready-to-deploy with Docker and Fly.io.
- Used Vercel for FrontEnd deployment.

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

## 📊 Datasets

- **TDC**(Therapeutics Data Commons): [TDC](https://tdcommons.ai/single_pred_tasks/overview/) - Used for Model -> training, testing and validating.
- **PubChem**: [PubChem](https://pubchem.ncbi.nlm.nih.gov/) - Used for Batch Predictions.

---

## 🛠 Tech Stack

- **Backend**: Python, Flask, Joblib, Flask-CORS
- **Frontend**: Vite-React, TailwindCSS, Framer-Motion, Lottie
- **AI/ML Models**: Custom trained models for ADMET prediction
- **Deployment**: Docker, Fly.io and Vercel
- **Utilities**: RDKit, ChemUtils, Plotting modules

---

### 📁 Project Files & Navigation
```bash
ADMET-X/
├── BackEnd/                # Flask backend code, app.py, utils, Models folder
├── FrontEnd/               # React frontend code
├── Model_predictions/      # Predicted ADMET results (optional storage)
├── Model_training/         # Scripts and notebooks for training ML models
├── Test_Model/             # Unit tests or test scripts for models
├── admet_data/             # Raw datasets used for training/testing
├── LICENSE                 # MIT License file
├── README.md               # Project README with badges, instructions, contributors
├── SECURITY.md             # Security policy and responsible disclosure
├── example.py              # Example script demonstrating usage
├── package-lock.json       # Frontend dependency lock file
├── package.json            # Frontend dependency definitions
└── paths.py                # Paths configuration for project directories/files
```

---

## ⚙️ Installation (Local Setup)

### **1. Clone the repository**
```bash
git clone https://github.com/yourusername/ADMET-X.git
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

---

### 👥 Contributors
| Name                 | GitHub                                     | LinkedIn                                               |
| ---------------------| ------------------------------------------ | ------------------------------------------------------ |
| Sheik Arshad Ibrahim | [GitHub](https://github.com/arshadibrahim882) | [LinkedIn](https://www.linkedin.com/in/sheik-arshad-ibrahim-33a16829a/) |
| Rohith Reddy G K     | [GitHub](https://github.com/RohithReddyGK)  | [LinkedIn](https://www.linkedin.com/in/rohithreddygk/)  |
| Sayed Jahangir Ali   | [GitHub](https://github.com/Jahangir-ali-74)      |  -  |
| Thirumurugan M       | [GitHub](https://github.com/thirumuruganmeganath-ops)      |  -     |

---

## 🧪 Contribution

Contributions are welcome! To contribute:
- Fork the repository
- Create a branch: git checkout -b feature-name
- Make changes and commit: git commit -m "Add new feature"
- Push to branch: git push origin feature-name
- Open a Pull Request

---

## 🌟If you liked our project, git it a ⭐.
