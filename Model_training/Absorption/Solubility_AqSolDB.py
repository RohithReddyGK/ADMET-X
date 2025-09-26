# Solubility_AqSolDB.py (Robust Version)
'''
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, Lipinski
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import classification_report, roc_auc_score, r2_score, mean_absolute_error, mean_squared_error
import joblib
import numpy as np

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
        return "classification", "-"
    elif "solubility" in str(df.columns).lower() or "solubility" in str(df).lower():
        return "regression", "log(mol/L)"
    else:
        return "classification", "-"

# ---------- Descriptor Computation (Robust) ----------
def compute_descriptors_for_dataframe(df, smiles_col="Drug"):
    desc_list = []
    valid_indices = []
    invalid_smiles = []

    for i, s in enumerate(df[smiles_col].astype(str)):
        desc = compute_descriptors(s)
        if desc is not None:
            desc_list.append(desc)
            valid_indices.append(i)
        else:
            invalid_smiles.append(s)

    if invalid_smiles:
        print(f"Skipped {len(invalid_smiles)} invalid SMILES. Examples: {invalid_smiles[:5]}")

    desc_df = pd.DataFrame(desc_list)
    df_valid = df.iloc[valid_indices].reset_index(drop=True)
    df_final = pd.concat([df_valid, desc_df], axis=1)
    return df_final

# ---------- Train Model ----------
def train_model(train_csv: str, smiles_col="Drug", label_col="Y",
                save_model_path=r"D:\VS Code Editor\Major Project\ADMET\Model_training\Absorption\Solubility_AqSolDB.joblib"):

    df = pd.read_csv(train_csv)
    task_type, units = detect_task_type(df, label_col)

    df = compute_descriptors_for_dataframe(df, smiles_col)

    feature_cols = ["MolWt", "LogP", "NumHDonors", "NumHAcceptors", "TPSA", "NumRotatableBonds", "AtomCount"]
    X, y = df[feature_cols].values, df[label_col].values

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
def evaluate_model(model, feature_cols, test_csv, smiles_col="Drug", label_col="Y", task_type="classification"):
    df = pd.read_csv(test_csv)
    df = compute_descriptors_for_dataframe(df, smiles_col)
    X, y = df[feature_cols].values, df[label_col].values

    if task_type == "classification":
        y = y.astype(int)
        y_pred = model.predict(X)
        if len(np.unique(y)) == 2:
            # Binary classification: ROC-AUC
            y_prob = model.predict_proba(X)[:, 1]
            print("Evaluation on Test Set (Binary Classification)")
            print(classification_report(y, y_pred))
            print("ROC-AUC:", roc_auc_score(y, y_prob))
        else:
            # Multi-class classification: ROC-AUC requires multi_class parameter
            print("Evaluation on Test Set (Multi-class Classification)")
            print(classification_report(y, y_pred))
    else:
        y_pred = model.predict(X)
        print("Evaluation on Test Set (Regression)")
        print(f"R² Score: {r2_score(y, y_pred):.4f}")
        print(f"MAE: {mean_absolute_error(y, y_pred):.4f}")
        print(f"MSE: {mean_squared_error(y, y_pred):.4f}")


# ---------- Predict and Save ----------
def predict_and_save(model, feature_cols, valid_csv, out_csv, smiles_col="Drug", task_type="classification", units="-"):
    df = pd.read_csv(valid_csv)
    df = compute_descriptors_for_dataframe(df, smiles_col)
    X = df[feature_cols].values

    if task_type == "classification":
        df["Prediction"] = model.predict_proba(X)[:, 1]
    else:
        df["Prediction"] = model.predict(X)

    druglike = df[smiles_col].astype(str).apply(lambda s: evaluate_druglikeness_from_desc(compute_descriptors(s)))
    druglike_df = pd.DataFrame(druglike.tolist(), index=df.index)
    df = pd.concat([df, druglike_df], axis=1)

    df["Units"] = units
    df.to_csv(out_csv, index=False)
    print(f"Predictions saved to {out_csv}")
    return df

# ---------- Main ----------
if __name__ == "__main__":
    train_csv = r"D:\VS Code Editor\Major Project\ADMET\admet_data\Absorption\Solubility_AqSolDB\train.csv"
    test_csv  = r"D:\VS Code Editor\Major Project\ADMET\admet_data\Absorption\Solubility_AqSolDB\test.csv"
    valid_csv = r"D:\VS Code Editor\Major Project\ADMET\admet_data\Absorption\Solubility_AqSolDB\valid.csv"
    save_model_path = r"D:\VS Code Editor\Major Project\ADMET\Model_training\Absorption\Solubility_AqSolDB.joblib"
    out_csv = r"D:\VS Code Editor\Major Project\ADMET\Model_predictions\Absorption\Solubility_AqSolDB.csv"

    model, feature_cols, task_type, units = train_model(train_csv, save_model_path=save_model_path)
    evaluate_model(model, feature_cols, test_csv, task_type=task_type)
    df_valid = predict_and_save(model, feature_cols, valid_csv, out_csv, task_type=task_type, units=units)
    print(df_valid.head())
    '''