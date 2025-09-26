# PAMPA_NCATS.py

import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors
import numpy as np
from rdkit.Chem import rdFingerprintGenerator
from sklearn.metrics import accuracy_score, classification_report, roc_auc_score
import xgboost as xgb
import joblib

# ---------- Morgan Fingerprints (New RDKit Standard) ----------
morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def compute_morgan_fingerprint(smiles: str, radius=2, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    arr = np.zeros((n_bits,), dtype=int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

# ---------- Rule-based Filters ----------
def check_lipinski(mol):
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    h_donors = Descriptors.NumHDonors(mol)
    h_acceptors = Descriptors.NumHAcceptors(mol)
    return "Yes" if (mw <= 500 and logp <= 5 and h_donors <= 5 and h_acceptors <= 10) else "No"

def check_ghose(mol):
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    atoms = mol.GetNumAtoms()
    return "Yes" if (160 <= mw <= 480 and -0.4 <= logp <= 5.6 and 20 <= atoms <= 70) else "No"

def check_veber(mol):
    tpsa = Descriptors.TPSA(mol)
    rot_bonds = Descriptors.NumRotatableBonds(mol)
    return "Yes" if (tpsa <= 140 and rot_bonds <= 10) else "No"

def check_egan(mol):
    tpsa = Descriptors.TPSA(mol)
    logp = Descriptors.MolLogP(mol)
    return "Yes" if (tpsa <= 131 and logp <= 5.88) else "No"

def check_muegge(mol):
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    h_donors = Descriptors.NumHDonors(mol)
    h_acceptors = Descriptors.NumHAcceptors(mol)
    rot_bonds = Descriptors.NumRotatableBonds(mol)
    tpsa = Descriptors.TPSA(mol)
    return "Yes" if (200 <= mw <= 600 and -2 <= logp <= 5 and h_donors <= 5 and
                     h_acceptors <= 10 and rot_bonds <= 15 and tpsa <= 150) else "No"

# ---------- Train Model ----------
def train_model(train_csv: str, smiles_col="Drug", label_col="Y", save_model_path=None):
    df = pd.read_csv(train_csv)
    df = df.dropna(subset=[smiles_col, label_col])

    fps = [compute_morgan_fingerprint(s) for s in df[smiles_col].astype(str)]
    fps = np.array([fp for fp in fps if fp is not None])
    y = df[label_col].values[:len(fps)]

    model = xgb.XGBClassifier(
        n_estimators=500,
        max_depth=8,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.8,
        random_state=42,
        n_jobs=-1,
        use_label_encoder=False,
        eval_metric="logloss"
    )
    model.fit(fps, y)

    if save_model_path:
        joblib.dump({"model": model, "smiles_col": smiles_col, "label_col": label_col}, save_model_path)
        print(f"Model trained and saved to {save_model_path}")

    return model

# ---------- Evaluate Model ----------
def evaluate_model(model, test_csv, smiles_col="Drug", label_col="Y"):
    df = pd.read_csv(test_csv)
    df = df.dropna(subset=[smiles_col, label_col])

    fps = [compute_morgan_fingerprint(s) for s in df[smiles_col].astype(str)]
    fps = np.array([fp for fp in fps if fp is not None])
    y_true = df[label_col].values[:len(fps)]

    y_pred = model.predict(fps)
    y_prob = model.predict_proba(fps)

    print("Evaluation on Test Set (Classification)")
    print(f"Accuracy: {accuracy_score(y_true, y_pred):.4f}")
    print(f"ROC-AUC: {roc_auc_score(y_true, y_prob[:,1]):.4f}")
    print("Classification Report:")
    print(classification_report(y_true, y_pred))

# ---------- Predict and Save ----------
def predict_and_save(model, valid_csv, out_csv, smiles_col="Drug"):
    df = pd.read_csv(valid_csv)
    df = df.dropna(subset=[smiles_col])

    fps = [compute_morgan_fingerprint(s) for s in df[smiles_col].astype(str)]
    fps = np.array([fp for fp in fps if fp is not None])

    y_pred_class = model.predict(fps)
    y_pred_prob = model.predict_proba(fps)

    df["Predicted_Class"] = y_pred_class
    df["Prob_Class_0"] = y_pred_prob[:,0]
    df["Prob_Class_1"] = y_pred_prob[:,1]
    df["Units"] = "-"

    # Add molecular descriptors
    molwt, logp, donors, acceptors, tpsa, rot, atoms = [], [], [], [], [], [], []
    lipinski, ghose, veber, egan, muegge = [], [], [], [], []

    for smi in df[smiles_col]:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            molwt.append(np.nan)
            logp.append(np.nan)
            donors.append(np.nan)
            acceptors.append(np.nan)
            tpsa.append(np.nan)
            rot.append(np.nan)
            atoms.append(np.nan)
            lipinski.append("NA")
            ghose.append("NA")
            veber.append("NA")
            egan.append("NA")
            muegge.append("NA")
        else:
            molwt.append(Descriptors.MolWt(mol))
            logp.append(Descriptors.MolLogP(mol))
            donors.append(Descriptors.NumHDonors(mol))
            acceptors.append(Descriptors.NumHAcceptors(mol))
            tpsa.append(Descriptors.TPSA(mol))
            rot.append(Descriptors.NumRotatableBonds(mol))
            atoms.append(mol.GetNumAtoms())
            lipinski.append(check_lipinski(mol))
            ghose.append(check_ghose(mol))
            veber.append(check_veber(mol))
            egan.append(check_egan(mol))
            muegge.append(check_muegge(mol))

    df["MolWt"] = molwt
    df["LogP"] = logp
    df["NumHDonors"] = donors
    df["NumHAcceptors"] = acceptors
    df["TPSA"] = tpsa
    df["NumRotatableBonds"] = rot
    df["AtomCount"] = atoms
    df["Lipinski"] = lipinski
    df["Ghose"] = ghose
    df["Veber"] = veber
    df["Egan"] = egan
    df["Muegge"] = muegge

    df.to_csv(out_csv, index=False)
    print(f"Predictions saved to {out_csv}")
    return df

# ---------- Main ----------
if __name__ == "__main__":
    train_path = r"C:\Users\Rohith Reddy G K\Dropbox\ADMET\admet_data\Absorption\PAMPA_NCATS\train.csv"
    test_path = r"C:\Users\Rohith Reddy G K\Dropbox\ADMET\admet_data\Absorption\PAMPA_NCATS\test.csv"
    valid_path = r"C:\Users\Rohith Reddy G K\Dropbox\ADMET\admet_data\Absorption\PAMPA_NCATS\valid.csv"
    save_pred_path = r"C:\Users\Rohith Reddy G K\Dropbox\ADMET\Predictions\PAMPA_NCATS.csv"
    save_model_path = r"C:\Users\Rohith Reddy G K\Dropbox\ADMET\Models\PAMPA_NCATS.joblib"

    model = train_model(train_path, save_model_path=save_model_path)
    evaluate_model(model, test_path)
    df_valid = predict_and_save(model, valid_path, save_pred_path)
    print(df_valid.head())
