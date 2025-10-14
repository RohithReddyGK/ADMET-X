# herg_central.py

import os
import numpy as np
import pandas as pd
import joblib
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, Crippen, rdMolDescriptors, Lipinski
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from xgboost import XGBClassifier, XGBRegressor
from sklearn.metrics import classification_report, roc_auc_score, mean_squared_error, r2_score
from paths import train_path, valid_path, test_path, prediction_path, model_path

# ---------- Fingerprint size ----------
FP_SIZE = 2048

# ---------- Feature builders ----------
def descriptors_from_mol(mol):
    return [
        float(Descriptors.MolWt(mol)),
        float(Crippen.MolLogP(mol)),
        int(Lipinski.NumHDonors(mol)),
        int(Lipinski.NumHAcceptors(mol)),
        float(rdMolDescriptors.CalcTPSA(mol)),
        int(rdMolDescriptors.CalcNumRotatableBonds(mol)),
        int(mol.GetNumAtoms())
    ]

gen = GetMorganGenerator(radius=2, fpSize=FP_SIZE)

def fingerprint_from_mol(mol):
    fp = gen.GetFingerprint(mol)
    arr = np.zeros((FP_SIZE,), dtype=int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

def build_features_and_desc(df, smiles_col="Drug", use_descriptors=True, use_fingerprints=True):
    feats, desc_rows, valid_idx = [], {}, []
    for idx, smi in df[smiles_col].astype(str).items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        parts = []
        if use_descriptors:
            d = descriptors_from_mol(mol)
            parts.append(np.array(d, dtype=float))
            desc_rows[idx] = d
        if use_fingerprints:
            arr = fingerprint_from_mol(mol)
            parts.append(arr.astype(float))
        feat = np.concatenate(parts) if parts else np.array([], dtype=float)
        feats.append(feat)
        valid_idx.append(idx)
    X = np.vstack(feats) if feats else np.empty((0,0))
    desc_df = pd.DataFrame.from_dict(
        desc_rows, orient="index",
        columns=["MolWt","LogP","NumHDonors","NumHAcceptors","TPSA","NumRotatableBonds","AtomCount"]
    ) if use_descriptors else None
    return X, desc_df, valid_idx

# ---------- Druglikeness rules ----------
def lipinski_rule(mol):
    return (Descriptors.MolWt(mol) <= 500 and Crippen.MolLogP(mol) <= 5 and
            Lipinski.NumHDonors(mol) <= 5 and Lipinski.NumHAcceptors(mol) <= 10)

def ghose_rule(mol):
    mw, logp, atoms = Descriptors.MolWt(mol), Crippen.MolLogP(mol), mol.GetNumAtoms()
    return (160 <= mw <= 480 and -0.4 <= logp <= 5.6 and 20 <= atoms <= 70)

def veber_rule(mol):
    return (rdMolDescriptors.CalcTPSA(mol) <= 140 and rdMolDescriptors.CalcNumRotatableBonds(mol) <= 10)

def egan_rule(mol):
    return (rdMolDescriptors.CalcTPSA(mol) <= 131 and Crippen.MolLogP(mol) <= 5.88)

def muegge_rule(mol):
    mw, logp = Descriptors.MolWt(mol), Crippen.MolLogP(mol)
    return (200 <= mw <= 600 and -2 <= logp <= 5 and
            Lipinski.NumHAcceptors(mol) <= 10 and Lipinski.NumHDonors(mol) <= 5 and
            rdMolDescriptors.CalcTPSA(mol) <= 150)

# ---------- Scientific units mapping ----------
ADMET_UNITS = {
    "hERG": "μM",
    "Caco-2": "nm/s",
    "Ames": "-",
    "BBB": "-",
    "VDss": "L/kg",
    "CL": "mL/min/kg",
    "Half-life": "h",
    "F20": "%",
    "F30": "%",
    # Add more datasets here as needed
}

# ---------- Train model ----------
def train_model(train_csv, valid_csv, smiles_col="Drug", label_col="Y",
                save_model_path=None, use_descriptors=True, use_fingerprints=True, model_type="auto"):

    train_df = pd.read_csv(train_csv)
    valid_df = pd.read_csv(valid_csv)

    X_train, _, idx_train = build_features_and_desc(train_df, smiles_col, use_descriptors, use_fingerprints)
    y_train = train_df.loc[idx_train, label_col].values

    X_valid, _, idx_valid = build_features_and_desc(valid_df, smiles_col, use_descriptors, use_fingerprints)
    y_valid = valid_df.loc[idx_valid, label_col].values

    dataset_name = os.path.basename(train_csv).split("_")[0]

    # Determine classification vs regression
    is_classification = False
    if np.issubdtype(y_train.dtype, np.integer):
        unique_vals = np.unique(y_train)
        is_classification = len(unique_vals) <= 10
    else:
        # If floats close to 0/1 → classify as classification
        if set(np.unique(y_train)) <= {0.0,1.0} or np.all(np.isclose(y_train, y_train.astype(int))):
            is_classification = True
            y_train = y_train.astype(int)
            y_valid = y_valid.astype(int)
        else:
            is_classification = False

    print(f"Task detected: {'Classification' if is_classification else 'Regression'} ({dataset_name})")

    if model_type=="auto":
        model_type = "xgb"

    # Instantiate model
    if is_classification:
        if model_type=="rf":
            model = RandomForestClassifier(n_estimators=300, random_state=42, class_weight="balanced", n_jobs=-1)
        else:  # XGB
            pos, neg = np.sum(y_train==1), np.sum(y_train==0)
            scale_pos_weight = float(neg) / max(1.0, float(pos))
            model = XGBClassifier(
                n_estimators=500, learning_rate=0.05, max_depth=6,
                subsample=0.8, colsample_bytree=0.8,
                scale_pos_weight=scale_pos_weight,
                use_label_encoder=False, eval_metric="logloss",
                random_state=42, n_jobs=-1
            )
    else:
        if model_type=="rf":
            model = RandomForestRegressor(n_estimators=300, random_state=42, n_jobs=-1)
        else:
            model = XGBRegressor(
                n_estimators=500, learning_rate=0.05, max_depth=6,
                subsample=0.8, colsample_bytree=0.8,
                random_state=42, n_jobs=-1
            )

    print(f"Training using {model_type}...")
    model.fit(X_train, y_train)

    # Validation results
    if X_valid.shape[0] > 0:
        yv_pred = model.predict(X_valid)
        if is_classification:
            try:
                yv_prob = model.predict_proba(X_valid)[:,1]
                auc = roc_auc_score(y_valid, yv_prob)
            except Exception:
                auc = None
            print("\nValidation results:")
            print(classification_report(y_valid, yv_pred))
            print("ROC-AUC (valid):", auc)
        else:
            mse, r2 = mean_squared_error(y_valid, yv_pred), r2_score(y_valid, yv_pred)
            print("\nValidation results:")
            print("MSE:", mse)
            print("R²:", r2)

    meta = {
        "model": model,
        "use_descriptors": use_descriptors,
        "use_fingerprints": use_fingerprints,
        "feature_order": (
            ["MolWt","LogP","NumHDonors","NumHAcceptors","TPSA","NumRotatableBonds","AtomCount"] if use_descriptors else []
        ) + (["FP_%d"%i for i in range(FP_SIZE)] if use_fingerprints else []),
        "model_type": model_type,
        "is_classification": is_classification,
        "dataset_name": dataset_name
    }

    if save_model_path:
        os.makedirs(os.path.dirname(save_model_path), exist_ok=True)
        joblib.dump(meta, save_model_path)
        print(f"\nModel saved to: {save_model_path}")

    return meta

# ---------- Evaluate and save predictions ----------
def evaluate_and_save_predictions(model_meta, test_csv, prediction_csv=None, smiles_col="Drug", label_col="Y"):
    df_test = pd.read_csv(test_csv)
    use_desc, use_fp = model_meta["use_descriptors"], model_meta["use_fingerprints"]
    model, is_classification = model_meta["model"], model_meta["is_classification"]
    dataset_name = model_meta.get("dataset_name", "ADMET")

    X_test, desc_test_df, idx_test = build_features_and_desc(df_test, smiles_col, use_desc, use_fp)
    result = df_test.copy()

    desc_cols = ["MolWt","LogP","NumHDonors","NumHAcceptors","TPSA","NumRotatableBonds","AtomCount"]
    for c in desc_cols: result[c] = np.nan
    if desc_test_df is not None and not desc_test_df.empty:
        for c in desc_test_df.columns:
            result.loc[desc_test_df.index, c] = desc_test_df[c].values

    # Predictions
    unit = ADMET_UNITS.get(dataset_name, "-")
    if X_test.shape[0] > 0:
        y_pred = model.predict(X_test)
        if is_classification:
            result["Predicted_Class"] = np.nan
            try: probs = model.predict_proba(X_test)
            except Exception: probs = np.full((X_test.shape[0],2), np.nan)
            classes = list(getattr(model, "classes_", [0,1]))
            for cls in classes: result[f"Prob_Class_{cls}"] = np.nan
            for i, idx in enumerate(idx_test):
                result.at[idx, "Predicted_Class"] = y_pred[i]
                if probs.ndim==2:
                    for j, cls in enumerate(classes):
                        result.at[idx, f"Prob_Class_{cls}"] = probs[i,j]
        else:
            result["Predicted_Value"] = np.nan
            for i, idx in enumerate(idx_test):
                result.at[idx, "Predicted_Value"] = y_pred[i]
    else:
        if is_classification: result["Predicted_Class"] = np.nan
        else: result["Predicted_Value"] = np.nan

    # Druglikeness rules
    for rule_name, rule_func in zip(
        ["Lipinski","Ghose","Veber","Egan","Muegge"],
        [lipinski_rule, ghose_rule, veber_rule, egan_rule, muegge_rule]
    ):
        result[rule_name] = [
            "Invalid" if Chem.MolFromSmiles(smi) is None else ("Yes" if rule_func(Chem.MolFromSmiles(smi)) else "No")
            for smi in result[smiles_col].astype(str)
        ]

    result["Units"] = unit

    # Save CSV
    if prediction_csv:
        os.makedirs(os.path.dirname(prediction_csv), exist_ok=True)
        result.to_csv(prediction_csv, index=False)
        print(f"\nPredictions saved to: {prediction_csv}")
        print(result.head(5))

    return result

# ---------- Example run ----------
if __name__ == "__main__":
    category = "Toxicity"
    model_name = "herg_central"

    train_csv = train_path(category, model_name)
    valid_csv = valid_path(category, model_name)
    test_csv  = test_path(category, model_name)
    save_model_file = model_path(category, model_name)
    prediction_csv = prediction_path(category, model_name)

    meta = train_model(train_csv, valid_csv, use_descriptors=True, use_fingerprints=True, model_type="auto",
                       save_model_path=save_model_file)

    evaluate_and_save_predictions(meta, test_csv, prediction_csv)
