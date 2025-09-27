# Solubility_AqSolDB.py

import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem import rdFingerprintGenerator
import numpy as np
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
import xgboost as xgb
import joblib

# ---------- Morgan Fingerprints ----------
morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def compute_morgan_fingerprint(smiles: str, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = morgan_gen.GetFingerprint(mol)
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
    fps = [compute_morgan_fingerprint(s) for s in df[smiles_col].astype(str)]
    valid_idx = [i for i, fp in enumerate(fps) if fp is not None]
    fps = np.array([fps[i] for i in valid_idx])
    df = df.iloc[valid_idx]
    y = df[label_col].values

    model = xgb.XGBRegressor(
        n_estimators=500,
        max_depth=8,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.8,
        random_state=42,
        n_jobs=-1
    )
    model.fit(fps, y)

    if save_model_path:
        joblib.dump({"model": model, "smiles_col": smiles_col, "label_col": label_col}, save_model_path)
        print(f"Model trained and saved to {save_model_path}")

    return model

# ---------- Evaluate Model ----------
def evaluate_model(model, test_csv, smiles_col="Drug", label_col="Y"):
    df = pd.read_csv(test_csv)
    fps = [compute_morgan_fingerprint(s) for s in df[smiles_col].astype(str)]
    valid_idx = [i for i, fp in enumerate(fps) if fp is not None]
    fps = np.array([fps[i] for i in valid_idx])
    df = df.iloc[valid_idx]
    y = df[label_col].values

    y_pred = model.predict(fps)

    print("Evaluation on Test Set (Regression)")
    print(f"R² Score: {r2_score(y, y_pred):.4f}")
    print(f"MAE: {mean_absolute_error(y, y_pred):.4f}")
    print(f"MSE: {mean_squared_error(y, y_pred):.4f}")

# ---------- Predict and Save ----------
def predict_and_save(model, valid_csv, out_csv, smiles_col="Drug"):
    df = pd.read_csv(valid_csv)
    fps = [compute_morgan_fingerprint(s) for s in df[smiles_col].astype(str)]
    valid_idx = [i for i, fp in enumerate(fps) if fp is not None]
    fps = np.array([fps[i] for i in valid_idx])
    df = df.iloc[valid_idx]

    df["Prediction"] = model.predict(fps)
    df["Units"] = "continuous"

    # Compute descriptors for the final CSV
    df["MolWt"] = df[smiles_col].apply(lambda x: Descriptors.MolWt(Chem.MolFromSmiles(x)) if Chem.MolFromSmiles(x) else np.nan)
    df["LogP"] = df[smiles_col].apply(lambda x: Descriptors.MolLogP(Chem.MolFromSmiles(x)) if Chem.MolFromSmiles(x) else np.nan)
    df["NumHDonors"] = df[smiles_col].apply(lambda x: Descriptors.NumHDonors(Chem.MolFromSmiles(x)) if Chem.MolFromSmiles(x) else np.nan)
    df["NumHAcceptors"] = df[smiles_col].apply(lambda x: Descriptors.NumHAcceptors(Chem.MolFromSmiles(x)) if Chem.MolFromSmiles(x) else np.nan)
    df["TPSA"] = df[smiles_col].apply(lambda x: Descriptors.TPSA(Chem.MolFromSmiles(x)) if Chem.MolFromSmiles(x) else np.nan)
    df["NumRotatableBonds"] = df[smiles_col].apply(lambda x: Descriptors.NumRotatableBonds(Chem.MolFromSmiles(x)) if Chem.MolFromSmiles(x) else np.nan)
    df["AtomCount"] = df[smiles_col].apply(lambda x: Chem.MolFromSmiles(x).GetNumAtoms() if Chem.MolFromSmiles(x) else np.nan)

    # Rule filters
    lipinski, ghose, veber, egan, muegge = [], [], [], [], []
    for smi in df[smiles_col]:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            lipinski.append("NA")
            ghose.append("NA")
            veber.append("NA")
            egan.append("NA")
            muegge.append("NA")
        else:
            lipinski.append(check_lipinski(mol))
            ghose.append(check_ghose(mol))
            veber.append(check_veber(mol))
            egan.append(check_egan(mol))
            muegge.append(check_muegge(mol))

    df["Lipinski"] = lipinski
    df["Ghose"] = ghose
    df["Veber"] = veber
    df["Egan"] = egan
    df["Muegge"] = muegge

    # Predicted class: Yes if >=3 rules passed
    df["Predicted_Class"] = df[["Lipinski","Ghose","Veber","Egan","Muegge"]].apply(
        lambda row: "Yes" if list(row).count("Yes") >= 3 else "No", axis=1
    )

    # Probabilities placeholder (since regression, just normalized prediction)
    df["Prob_Class_0"] = 1 - df["Prediction"] / df["Prediction"].max()
    df["Prob_Class_1"] = df["Prediction"] / df["Prediction"].max()

    # Reorder columns to match your required format
    columns_order = ["Drug_ID","Drug","Y","MolWt","LogP","NumHDonors","NumHAcceptors",
                     "TPSA","NumRotatableBonds","AtomCount","Prob_Class_0","Prob_Class_1",
                     "Predicted_Class","Lipinski","Ghose","Veber","Egan","Muegge","Units"]
    df = df[columns_order]

    df.to_csv(out_csv, index=False)
    print(f"Predictions saved to {out_csv}")
    return df

# ---------- Main ----------
if __name__ == "__main__":
    train_path = r"C:\Users\Rohith Reddy G K\Dropbox\ADMET\admet_data\Absorption\Solubility_AqSolDB\train.csv"
    test_path = r"C:\Users\Rohith Reddy G K\Dropbox\ADMET\admet_data\Absorption\Solubility_AqSolDB\test.csv"
    valid_path = r"C:\Users\Rohith Reddy G K\Dropbox\ADMET\admet_data\Absorption\Solubility_AqSolDB\valid.csv"
    save_pred_path = r"C:\Users\Rohith Reddy G K\Dropbox\ADMET\Model_predictions\Absorption\Solubility_AqSolDB.csv"
    save_model_path = r"C:\Users\Rohith Reddy G K\Dropbox\ADMET\Model_training\Absorption\Solubility_AqSolDB.joblib"

    model = train_model(train_path, save_model_path=save_model_path)
    evaluate_model(model, test_path)
    df_valid = predict_and_save(model, valid_path, save_pred_path)
    print(df_valid.head())
