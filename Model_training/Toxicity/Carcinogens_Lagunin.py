# Carcinogens_Lagunin.py

import os
import numpy as np
import pandas as pd
import joblib
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, Crippen, rdMolDescriptors, Lipinski
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.metrics import classification_report, roc_auc_score
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
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

# Morgan fingerprint generator
gen = GetMorganGenerator(radius=2, fpSize=FP_SIZE)

def fingerprint_from_mol(mol):
    fp = gen.GetFingerprint(mol)
    arr = np.zeros((FP_SIZE,), dtype=int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

def build_features_and_desc(df, smiles_col="Drug", use_descriptors=True, use_fingerprints=True):
    feats = []
    desc_rows = {}
    valid_idx = []

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

    X = np.vstack(feats) if feats else np.empty((0, 0))

    desc_df = pd.DataFrame.from_dict(desc_rows, orient="index",
                                     columns=["MolWt", "LogP", "NumHDonors", "NumHAcceptors",
                                              "TPSA", "NumRotatableBonds", "AtomCount"]) if use_descriptors else None

    return X, desc_df, valid_idx

# ---------- Druglikeness rules ----------
def lipinski_rule(mol):
    return (Descriptors.MolWt(mol) <= 500 and Crippen.MolLogP(mol) <= 5 and
            Lipinski.NumHDonors(mol) <= 5 and Lipinski.NumHAcceptors(mol) <= 10)

def ghose_rule(mol):
    mw = Descriptors.MolWt(mol); logp = Crippen.MolLogP(mol); atoms = mol.GetNumAtoms()
    return (160 <= mw <= 480 and -0.4 <= logp <= 5.6 and 20 <= atoms <= 70)

def veber_rule(mol):
    return (rdMolDescriptors.CalcTPSA(mol) <= 140 and rdMolDescriptors.CalcNumRotatableBonds(mol) <= 10)

def egan_rule(mol):
    return (rdMolDescriptors.CalcTPSA(mol) <= 131 and Crippen.MolLogP(mol) <= 5.88)

def muegge_rule(mol):
    mw = Descriptors.MolWt(mol); logp = Crippen.MolLogP(mol)
    return (200 <= mw <= 600 and -2 <= logp <= 5 and
            Lipinski.NumHAcceptors(mol) <= 10 and Lipinski.NumHDonors(mol) <= 5 and
            rdMolDescriptors.CalcTPSA(mol) <= 150)

# ---------- Train model ----------
def train_model(train_csv, valid_csv, smiles_col="Drug", label_col="Y",
                save_model_path=None, use_descriptors=True, use_fingerprints=False, model_type="rf"):

    train_df = pd.read_csv(train_csv)
    valid_df = pd.read_csv(valid_csv)

    X_train, desc_train_df, idx_train = build_features_and_desc(train_df, smiles_col, use_descriptors, use_fingerprints)
    y_train = train_df.loc[idx_train, label_col].astype(int).values

    X_valid, desc_valid_df, idx_valid = build_features_and_desc(valid_df, smiles_col, use_descriptors, use_fingerprints)
    y_valid = valid_df.loc[idx_valid, label_col].astype(int).values

    pos, neg = np.sum(y_train==1), np.sum(y_train==0)

    if model_type == "rf":
        model = RandomForestClassifier(n_estimators=300, random_state=42, class_weight="balanced", n_jobs=-1)
    elif model_type == "xgb":
        scale_pos_weight = float(neg) / max(1.0, float(pos))
        model = XGBClassifier(n_estimators=500, learning_rate=0.05, max_depth=6,
                              subsample=0.8, colsample_bytree=0.8,
                              scale_pos_weight=scale_pos_weight,
                              use_label_encoder=False, eval_metric="logloss",
                              random_state=42, n_jobs=-1)
    else:
        raise ValueError("model_type must be 'rf' or 'xgb'")

    print(f"Training using {model_type}...")
    model.fit(X_train, y_train)

    if X_valid.shape[0] > 0:
        yv_pred = model.predict(X_valid)
        try:
            yv_prob = model.predict_proba(X_valid)[:, 1]
            auc = roc_auc_score(y_valid, yv_prob)
        except Exception:
            yv_prob = None
            auc = None
        print("\nValidation results:")
        print(classification_report(y_valid, yv_pred))
        print("ROC-AUC (valid):", auc)

    meta = {
        "model": model,
        "use_descriptors": use_descriptors,
        "use_fingerprints": use_fingerprints,
        "feature_order": (["MolWt","LogP","NumHDonors","NumHAcceptors","TPSA","NumRotatableBonds","AtomCount"] if use_descriptors else []) +
                         (["FP_%d"%i for i in range(FP_SIZE)] if use_fingerprints else []),
        "model_type": model_type
    }

    if save_model_path:
        os.makedirs(os.path.dirname(save_model_path), exist_ok=True)
        joblib.dump(meta, save_model_path)
        print(f"\nModel saved to: {save_model_path}")

    return meta

# ---------- Evaluate and save predictions ----------
def evaluate_and_save_predictions(model_meta, test_csv, prediction_csv=None, smiles_col="Drug", label_col="Y"):
    df_test = pd.read_csv(test_csv)
    use_descriptors = model_meta.get("use_descriptors", True)
    use_fingerprints = model_meta.get("use_fingerprints", False)
    model = model_meta["model"]

    X_test, desc_test_df, idx_test = build_features_and_desc(df_test, smiles_col, use_descriptors, use_fingerprints)

    result = df_test.copy()
    desc_cols = ["MolWt","LogP","NumHDonors","NumHAcceptors","TPSA","NumRotatableBonds","AtomCount"]
    for c in desc_cols:
        result[c] = np.nan

    if desc_test_df is not None and not desc_test_df.empty:
        for c in desc_test_df.columns:
            result.loc[desc_test_df.index, c] = desc_test_df[c].values

    if X_test.shape[0] > 0:
        y_pred = model.predict(X_test)
        try:
            probs = model.predict_proba(X_test)
        except Exception:
            probs = np.full((X_test.shape[0], 2), np.nan)

        result["Predicted_Class"] = np.nan
        classes = list(getattr(model, "classes_", [0,1]))
        for cls in classes:
            result[f"Prob_Class_{cls}"] = np.nan

        for i, idx in enumerate(idx_test):
            result.at[idx, "Predicted_Class"] = y_pred[i]
            if probs.ndim==2:
                for j, cls in enumerate(classes):
                    result.at[idx, f"Prob_Class_{cls}"] = probs[i,j]

    else:
        result["Predicted_Class"] = np.nan
        result["Prob_Class_0"] = np.nan
        result["Prob_Class_1"] = np.nan

    # Add druglikeness rules
    for rule_name, rule_func in zip(["Lipinski","Ghose","Veber","Egan","Muegge"],
                                    [lipinski_rule, ghose_rule, veber_rule, egan_rule, muegge_rule]):
        result[rule_name] = [
            "Invalid" if Chem.MolFromSmiles(smi) is None else ("Yes" if rule_func(Chem.MolFromSmiles(smi)) else "No")
            for smi in result[smiles_col].astype(str)
        ]

    result["Units"] = "-"

    # Column order
    ordered_cols = ["Drug_ID","Drug",label_col]+desc_cols+["Prob_Class_0","Prob_Class_1","Predicted_Class",
                                                          "Lipinski","Ghose","Veber","Egan","Muegge","Units"]
    final_cols = [c for c in ordered_cols if c in result.columns] + [c for c in result.columns if c not in ordered_cols]

    if prediction_csv:
        os.makedirs(os.path.dirname(prediction_csv), exist_ok=True)
        result.to_csv(prediction_csv, index=False, columns=final_cols)
        print(f"\nPredictions saved to: {prediction_csv}")
        print(result[final_cols].head(5))

    return result

# ---------- Example run ----------
if __name__ == "__main__":
    category = "Toxicity"
    model_name = "Carcinogens_Lagunin"

    train_csv = train_path(category, model_name)
    valid_csv = valid_path(category, model_name)
    test_csv  = test_path(category, model_name)
    save_model_file = model_path(category, model_name)
    prediction_csv = prediction_path(category, model_name)

    meta = train_model(train_csv, valid_csv, use_descriptors=True, use_fingerprints=True, model_type="xgb",
                       save_model_path=save_model_file)

    evaluate_and_save_predictions(meta, test_csv, prediction_csv)