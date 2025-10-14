# Security Policy for ADMET-X

## 🚨 Reporting a Vulnerability

If you discover a security vulnerability in ADMET-X, please **report it responsibly**:

- **Email**: gkrohithreddy@gmail.com   
- **Subject Line**: `[ADMET-X Security] Vulnerability Report`

Please include:

- Steps to reproduce the issue
- Expected vs actual behavior
- Any relevant screenshots or logs
- Environment details (OS, Python version, frontend version)

We will acknowledge your report within **48 hours** and aim to provide a fix or guidance within **7 business days**.

---

## 🛡 Security Overview

ADMET-X is designed with security and privacy in mind. Key considerations include:

- **Data Handling**: SMILES input is processed **locally or on the deployed backend**. No data is shared externally.  
- **Dependencies**: All Python and JavaScript dependencies are regularly updated to prevent known vulnerabilities.  
- **Authentication**: The platform currently does not require user authentication. If added in future versions, it will follow secure practices.  
- **Deployment Security**: Recommended deployment using **Docker** and **Fly.io**, which provide container isolation and HTTPS enforcement.

---

## 🔐 Recommendations for Secure Usage

- Avoid uploading **sensitive or proprietary SMILES data** unless running on a private deployment.  
- Always use the **latest version** of the backend and frontend.  
- If self-hosting, configure HTTPS and firewall rules.  
- Regularly update all dependencies (`pip`/`npm`) to patch known vulnerabilities.

---

## 📝 Dependencies Security

- **Backend**:
  - Python 3.12
  - Flask 2.x
  - RDKit 2025
  - Joblib
  - Flask-CORS

- **Frontend**:
  - React 18.2
  - Vite 4.x
  - TailwindCSS 3.x
  - Framer-Motion
  - Lottie

> We recommend checking [Dependabot](https://docs.github.com/en/code-security/supply-chain-security/keeping-your-dependencies-updated-automatically/about-dependabot-version-updates) to monitor dependencies for vulnerabilities.

---

## 📅 Security Updates

All critical or high-severity security issues will be patched **promptly** and communicated via GitHub releases and email notifications to registered users (if applicable).  

---

## 🏆 Responsible Disclosure

By reporting security vulnerabilities responsibly, you help us protect the platform and all users. All reports will be credited in the `Contributors` section of the project.

---

*Thank you for helping make ADMET-X safe and secure!*
