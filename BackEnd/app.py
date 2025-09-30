from flask import Flask, request, jsonify
from flask_cors import CORS
import os
import joblib
from rdkit.Chem import Descriptors
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import numpy as np
from utils.generate_images import mol_to_image, radar_plot_image

# ---------- CONFIG ----------
MODELS_DIR = "Models"  # Folder where all ADMET models are stored
SCIENCE = ["Lipinski", "Ghose", "Veber", "Egan", "Muegge"]

PROPERTY_UNITS = {
    # Absorption
    "Bioavailability_Ma": "-",
    "Caco2_Wang": "log(10^-6 cm/s)",
    "HIA_Hou": "-",
    "HydrationFreeEnergy_FreeSolv": "kcal/mol",
    "Lipophilicity_AstraZeneca": "logP",
    "PAMPA_NCATS": "cm/s",
    "Pgp_Broccatelli": "-",
    "Solubility_AqSolDB": "log(mol/L)",
    # Distribution
    "BBB_MAartins": "-",
    "PPBR_AZ": "%",
    "VDss_Lombardo": "L/kg",
    # Metabolism
    "CYP1A2_Veith": "-",
    "CYP2C19_Veith": "-",
    "CYP2C9_Substrate_CarbonMangels": "-",
    "CYP2C9_Veith": "-",
    "CYP2D6_Substrate_CarbonMangels": "-",
    "CYP2D6_Veith": "-",
    "CYP3A4_Substrate_CarbonMangels": "-",
    "CYP3A4_Veith": "-",
    # Excretion
    "Clearance_Hepatocyte_AZ": "mL/min/kg",
    "Half_Life_Obach": "h",
    # Toxicity
    "AMES": "-",
    "Carcinogens_Lagunin": "-",
    "ClinTox": "-",
    "DILI": "-",
    "hERG": "-",
    "herg_central": "-",
    "hERG_Karim": "-",
    "LD50_Zhu": "mg/kg",
    "Skin_Reaction": "-",
    "Tox21": "-",
    "Toxcast": "-",
}

category_props = {
    "Absorption": ["Bioavailability_Ma","Caco2_Wang","HIA_Hou","HydrationFreeEnergy_FreeSolv",
                   "Lipophilicity_AstraZeneca","PAMPA_NCATS","Pgp_Broccatelli","Solubility_AqSolDB"],
    "Distribution": ["BBB_MAartins","PPBR_AZ","VDss_Lombardo"],
    "Metabolism": ["CYP1A2_Veith","CYP2C19_Veith","CYP2C9_Substrate_CarbonMangels","CYP2C9_Veith",
                   "CYP2D6_Substrate_CarbonMangels","CYP2D6_Veith","CYP3A4_Substrate_CarbonMangels","CYP3A4_Veith"],
    "Excretion": ["Clearance_Hepatocyte_AZ","Half_Life_Obach"],
    "Toxicity": ["AMES","Carcinogens_Lagunin","ClinTox","DILI","hERG","herg_central","hERG_Karim",
                 "LD50_Zhu","Skin_Reaction","Tox21","Toxcast"]
}

# ---------- LOAD MODELS ----------
ADMET_MODELS = {}
for category in category_props.keys():
    category_path = os.path.join(MODELS_DIR, category)
    if os.path.isdir(category_path):
        ADMET_MODELS[category] = {}
        for model_file in os.listdir(category_path):
            if model_file.endswith(".joblib"):
                model_name = model_file.replace(".joblib", "")
                model_path = os.path.join(category_path, model_file)
                ADMET_MODELS[category][model_name] = joblib.load(model_path)

# ---------- FLASK APP ----------
app = Flask(__name__)
CORS(app)  # Allow all origins for testing

# ---------- Helper Functions ----------
def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "MolWt": round(Descriptors.MolWt(mol), 2),
        "LogP": round(Descriptors.MolLogP(mol), 2),
        "H-Donors": Descriptors.NumHDonors(mol),
        "H-Acceptors": Descriptors.NumHAcceptors(mol),
        "TPSA": round(Descriptors.TPSA(mol), 2),
        "RotBonds": Descriptors.NumRotatableBonds(mol),
        "Atoms": mol.GetNumAtoms(),
    }

def descriptors_to_array(descriptors):
    return [
        descriptors["MolWt"],
        descriptors["LogP"],
        descriptors["H-Donors"],
        descriptors["H-Acceptors"],
        descriptors["TPSA"],
        descriptors["RotBonds"],
        descriptors["Atoms"],
    ]

def compute_morgan_fingerprint(smiles, radius=2, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    arr = np.zeros((n_bits,), dtype=int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

# ---------- ROUTES ----------
@app.route("/", methods=["GET"])
def home():
    return "Flask backend is running!"

@app.route("/predict", methods=["POST"])
def predict():
    try:
        data = request.get_json()
        smiles = data.get("smiles", "").strip()
        if not smiles:
            return jsonify({"error": "SMILES string required"}), 400

        descriptors = compute_descriptors(smiles)
        if descriptors is None:
            mol_image = mol_to_image("")  # blank image
            radar_image = radar_plot_image({
                "MolWt": 0, "LogP": 0, "H-Donors": 0,
                "H-Acceptors": 0, "TPSA": 0, "RotBonds": 0, "Atoms": 0
            })
            return jsonify({
                "smiles": smiles,
                "molImage": mol_image,
                "radarImage": radar_image,
                "descriptors": {},
                "ADMET": {},
                "error": "Invalid SMILES string"
            }), 400

        # Generate molecule image
        mol_image = mol_to_image(smiles)
        radar_image = radar_plot_image(descriptors)

        # ---------- ADMET Predictions ----------
        admet_results = {}

        for category, props in category_props.items():
            admet_results[category] = []
            for prop in props:
                model_dict = ADMET_MODELS.get(category, {}).get(prop)
                prediction = None

                units = PROPERTY_UNITS.get(prop, "-")
                
                if model_dict:
                    # Access the actual model inside the dict
                    model = model_dict['model'] if isinstance(model_dict, dict) else model_dict
    
                    try:
                        n_feat = getattr(model, "n_features_in_", None)
                        if n_feat == 7:
                            feature_array = descriptors_to_array(descriptors)
                        elif n_feat == 2048:
                            fp = compute_morgan_fingerprint(smiles)
                            feature_array = fp.tolist() if fp is not None else [0]*2048
                        elif n_feat == 2055:
                            fp = compute_morgan_fingerprint(smiles)
                            desc_array = descriptors_to_array(descriptors)
                            feature_array = desc_array + (fp.tolist() if fp is not None else [0]*2048)
                        else:
                            feature_array = descriptors_to_array(descriptors)

                        prediction = float(model.predict([feature_array])[0])
                    except Exception as e:
                        print(f"Prediction error for {prop}: {e}")
                        prediction = None

                druglikeness = [
                    {"scientist": s, "value": "Yes" if prediction is not None and prediction > 0.5 else "No"}
                    for s in SCIENCE
                ]

                admet_results[category].append({
                    "property": prop,
                    "prediction": prediction,
                    "druglikeness": druglikeness,
                    "units": units
                })

        response = {
            "smiles": smiles,
            "molImage": mol_image,
            "radarImage": radar_image,
            "descriptors": descriptors,
            "ADMET": admet_results
        }

        return jsonify(response)

    except Exception as e:
        print(f"Unexpected error: {e}")
        return jsonify({"error": "Server error occurred"}), 500

# ---------- RUN ----------
if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000, debug=True)
