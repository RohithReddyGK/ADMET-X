# HydrationFreeEnergy_FreeSolv.py

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, Lipinski
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import classification_report, roc_auc_score, r2_score, mean_absolute_error, mean_squared_error
import joblib   
import numpy as np
from paths import train_path, valid_path, test_path, prediction_path, model_path

# ---------- Helper: Compute Molecular Descriptors ----------
def compute_descriptors(smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "MolWt": float(Descriptors.MolWt(mol)),
        "LogP": float(Crippen.MolLogP(mol)),
        "NumHDonors": int(Lipinski.NumHDonors(mol)),
        "NumHAcceptors": int(Lipinski.NumHAcceptors(mol)),
        "TPSA": float(rdMolDescriptors.CalcTPSA(mol)),
        "NumRotatableBonds": int(rdMolDescriptors.CalcNumRotatableBonds(mol)),
        "AtomCount": int(mol.GetNumAtoms())
    }

# ---------- Druglikeness Evaluation ----------
def evaluate_druglikeness_from_desc(desc: dict) -> dict:
    if desc is None:
        return {"Lipinski": "Invalid", "Ghose": "Invalid", "Veber": "Invalid", "Egan": "Invalid", "Muegge": "Invalid"}

    violations = sum([
        desc["MolWt"] > 500,
        desc["LogP"] > 5,
        desc["NumHDonors"] > 5,
        desc["NumHAcceptors"] > 10
    ])
    lipinski = "Yes" if violations <= 1 else "No"
    ghose = "Yes" if (160 <= desc["MolWt"] <= 480 and -0.4 <= desc["LogP"] <= 5.6 and 40 <= desc["AtomCount"] <= 70) else "No"
    veber = "Yes" if (desc["NumRotatableBonds"] <= 10 and desc["TPSA"] <= 140) else "No"
    egan = "Yes" if (desc["LogP"] <= 5.88 and desc["TPSA"] <= 131) else "No"
    muegge = "Yes" if (200 <= desc["MolWt"] <= 600 and desc["LogP"] <= 5 and desc["NumRotatableBonds"] <= 15 and desc["TPSA"] <= 150) else "No"
    return {"Lipinski": lipinski, "Ghose": ghose, "Veber": veber, "Egan": egan, "Muegge": muegge}

# ---------- Auto-detect task type ----------
def detect_task_type(df: pd.DataFrame, target_col="Y"):
    if "bioavailability" in str(df.columns).lower() or "bioavailability" in str(df).lower():
        task_type = "classification"
        units = "-"
    elif "solubility" in str(df.columns).lower() or "solubility" in str(df).lower():
        task_type = "regression"
        units = "log(mol/L)"
    else:
        task_type = "classification"
        units = "-"
    print(f"Detected task: {task_type}, Units: {units}")
    return task_type, units

# ---------- Train Model ----------
def train_model(train_csv: str, smiles_col="Drug", label_col="Y",
                save_model_path=None):
    df = pd.read_csv(train_csv)

    task_type, units = detect_task_type(df, label_col)

    # Compute descriptors
    desc_list = [compute_descriptors(s) for s in df[smiles_col].astype(str)]
    desc_df = pd.DataFrame(desc_list)
    df = pd.concat([df, desc_df], axis=1).dropna()

    feature_cols = ["MolWt", "LogP", "NumHDonors", "NumHAcceptors", "TPSA", "NumRotatableBonds", "AtomCount"]
    X, y = df[feature_cols].values, df[label_col].values

    # Choose model based on task type
    if task_type == "classification":
        model = RandomForestClassifier(n_estimators=200, random_state=42)
        y = y.astype(int)
    else:
        model = RandomForestRegressor(n_estimators=200, random_state=42)

    model.fit(X, y)
    joblib.dump({"model": model, "feature_cols": feature_cols, "task_type": task_type, "units": units}, save_model_path)
    print(f"Model trained and saved to {save_model_path}")
    return model, feature_cols, task_type, units

# ---------- Evaluate Model ----------
def evaluate_model(model, feature_cols, test_csv, smiles_col="Drug", label_col="Y",
                   task_type="classification"):
    df = pd.read_csv(test_csv)
    desc_list = [compute_descriptors(s) for s in df[smiles_col].astype(str)]
    desc_df = pd.DataFrame(desc_list)
    df = pd.concat([df, desc_df], axis=1).dropna()

    X, y = df[feature_cols].values, df[label_col].values
    if task_type == "classification":
        y = y.astype(int)
        y_pred = model.predict(X)
        y_prob = model.predict_proba(X)[:, 1]
        print("Evaluation on Test Set (Classification)")
        print(classification_report(y, y_pred))
        print("ROC-AUC:", roc_auc_score(y, y_prob))
    else:
        y_pred = model.predict(X)
        print("Evaluation on Test Set (Regression)")
        print(f"R² Score: {r2_score(y, y_pred):.4f}")
        print(f"MAE: {mean_absolute_error(y, y_pred):.4f}")
        print(f"MSE: {mean_squared_error(y, y_pred):.4f}")

# ---------- Predict and Save ----------
def predict_and_save(model, feature_cols, valid_csv, out_csv, smiles_col="Drug",
                     task_type="classification", units="-"):
    df = pd.read_csv(valid_csv)
    desc_list = [compute_descriptors(s) for s in df[smiles_col].astype(str)]
    desc_df = pd.DataFrame(desc_list)
    df = pd.concat([df, desc_df], axis=1).dropna()

    X = df[feature_cols].values
    if task_type == "classification":
        df["Prediction"] = model.predict_proba(X)[:, 1]
    else:
        df["Prediction"] = model.predict(X)

    # Add druglikeness
    druglike = df[smiles_col].astype(str).apply(lambda s: evaluate_druglikeness_from_desc(compute_descriptors(s)))
    druglike_df = pd.DataFrame(druglike.tolist(), index=df.index)
    df = pd.concat([df, druglike_df], axis=1)

    df["Units"] = units
    df.to_csv(out_csv, index=False)
    print(f"Predictions saved to {out_csv}")
    return df

# ---------- Example Usage ----------
if __name__ == "__main__":
    category = "Absorption"
    model_name = "HydrationFreeEnergy_FreeSolv"

    train_csv = train_path(category, model_name)
    valid_csv = valid_path(category, model_name)
    test_csv  = test_path(category, model_name)
    out_csv   = prediction_path(category, model_name)
    save_model_file = model_path(category, model_name)

    model, feature_cols, task_type, units = train_model(train_csv, save_model_path=save_model_file)
    evaluate_model(model, feature_cols, test_csv, task_type=task_type)
    df_valid = predict_and_save_binary(model, feature_cols, valid_csv, out_csv, task_type=task_type, units=units)
    print(df_valid.head())