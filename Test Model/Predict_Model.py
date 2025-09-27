# Predict_Model.py

import pandas as pd
import numpy as np
import joblib
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors

# ---------- Morgan Fingerprint Generator ----------
def compute_morgan_fingerprint(smiles, radius=2, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    arr = np.zeros((n_bits,), dtype=int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

# ---------- Descriptor Generator (for classifier models) ----------
def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return [np.nan]*7  # Number of descriptors used in classifier.
    return [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        Descriptors.NumHDonors(mol),
        Descriptors.NumHAcceptors(mol),
        Descriptors.TPSA(mol),
        Descriptors.NumRotatableBonds(mol),
        mol.GetNumAtoms()
    ]

# ---------- Main Prediction Function ----------
def predict_dataset(model_path, dataset_path):
    # Load model
    model_data = joblib.load(model_path)
    model = model_data["model"]
    n_feat = getattr(model, "n_features_in_", None)

    # Load dataset
    df = pd.read_csv(dataset_path)

    # ---------- Generate features and make predictions ----------
    if n_feat == 2048:
        # Regression model
        fps = [compute_morgan_fingerprint(s) for s in df["Drug"]]
        X = np.array([fp for fp in fps if fp is not None])
        df["Prediction"] = model.predict(X)
        df["Units"] = "continuous"

    elif n_feat == 7:
        # Classifier model
        desc = df["Drug"].apply(compute_descriptors)
        X = np.vstack(desc.values)
        df["Predicted_Class"] = model.predict(X)
        y_pred_prob = model.predict_proba(X)
        df["Prob_Class_0"] = y_pred_prob[:,0]
        df["Prob_Class_1"] = y_pred_prob[:,1]
        df["Units"] = "-"
    else:
        raise ValueError(f"Model features ({n_feat}) do not match known types.")

    # ---------- Apply 5 Scientists' Rules for ALL models ----------
    lipinski, ghose, veber, egan, muegge = [], [], [], [], []
    for smi in df["Drug"]:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            lipinski.append("NA"); ghose.append("NA"); veber.append("NA")
            egan.append("NA"); muegge.append("NA")
        else:
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            h_donors = Descriptors.NumHDonors(mol)
            h_acceptors = Descriptors.NumHAcceptors(mol)
            tpsa = Descriptors.TPSA(mol)
            rot_bonds = Descriptors.NumRotatableBonds(mol)
            atoms = mol.GetNumAtoms()

            lipinski.append("Yes" if (mw <= 500 and logp <=5 and h_donors<=5 and h_acceptors<=10) else "No")
            ghose.append("Yes" if (160<=mw<=480 and -0.4<=logp<=5.6 and 20<=atoms<=70) else "No")
            veber.append("Yes" if (tpsa<=140 and rot_bonds<=10) else "No")
            egan.append("Yes" if (tpsa<=131 and logp<=5.88) else "No")
            muegge.append("Yes" if (200<=mw<=600 and -2<=logp<=5 and h_donors<=5 and h_acceptors<=10 and rot_bonds<=15 and tpsa<=150) else "No")

    df["Lipinski"] = lipinski
    df["Ghose"] = ghose
    df["Veber"] = veber
    df["Egan"] = egan
    df["Muegge"] = muegge

    # ---------- Display in terminal ----------
    print("\n--- Predictions ---")
    print(df.head(20))  # Shows first 20 rows.
    return df

# ---------- Run Example ----------
if __name__ == "__main__":

    model_path = r"C:\Users\Rohith Reddy G K\Dropbox\ADMET\Model_training\Absorption\Solubility_AqSolDB.joblib"
    dataset_path = r"C:\Users\Rohith Reddy G K\Dropbox\ADMET\admet_data\Absorption\Solubility_AqSolDB\valid.csv"

    predict_dataset(model_path, dataset_path)
