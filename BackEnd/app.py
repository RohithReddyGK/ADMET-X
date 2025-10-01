from flask import Flask, request, jsonify 
from flask_cors import CORS
import os
import joblib

from utils.chem_utils import compute_descriptors, mol_to_image
from utils.plotting import radar_plot_image
from utils.prediction_utils import safe_predict

# ---------- CONFIG ----------
MODELS_DIR = "Models"
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
    "BBB_Martins": "-",
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
    "Absorption": [
        "Bioavailability_Ma","Caco2_Wang","HIA_Hou","HydrationFreeEnergy_FreeSolv",
        "Lipophilicity_AstraZeneca","PAMPA_NCATS","Pgp_Broccatelli","Solubility_AqSolDB"
    ],
    "Distribution": ["BBB_Martins","PPBR_AZ","VDss_Lombardo"],
    "Metabolism": [
        "CYP1A2_Veith","CYP2C19_Veith","CYP2C9_Substrate_CarbonMangels","CYP2C9_Veith",
        "CYP2D6_Substrate_CarbonMangels","CYP2D6_Veith","CYP3A4_Substrate_CarbonMangels","CYP3A4_Veith"
    ],
    "Excretion": ["Clearance_Hepatocyte_AZ","Half_Life_Obach"],
    "Toxicity": [
        "AMES","Carcinogens_Lagunin","ClinTox","DILI","hERG","herg_central","hERG_Karim",
        "LD50_Zhu","Skin_Reaction","Tox21","Toxcast"
    ]
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
                try:
                    ADMET_MODELS[category][model_name] = joblib.load(model_path)
                except Exception as e:
                    print(f"⚠️ Failed to load {model_file}: {e}")

# ---------- THRESHOLD FUNCTION ----------
def categorize_property(prop, value):
    """
    Assign color status (green/yellow/red/gray) based on known thresholds.
    """
    if value is None:
        return "gray"

    try:
        v = float(value)
    except Exception:
        return "gray"

    # ---- Absorption ----
    if prop == "Bioavailability_Ma":
        return "green" if v >= 0.5 else "red"
    elif prop == "Caco2_Wang":
        return "green" if v > -5.15 else ("yellow" if -6 <= v <= -5.15 else "red")
    elif prop == "HIA_Hou":
        return "green" if v >= 0.5 else "red"
    elif prop == "HydrationFreeEnergy_FreeSolv":
        return "green" if v < -5 else ("yellow" if -5 <= v <= 0 else "red")
    elif prop == "Lipophilicity_AstraZeneca":
        return "green" if -0.7 <= v <= 5.0 else "red"
    elif prop == "PAMPA_NCATS":
        return "green" if v >= 1e-5 else ("yellow" if 1e-6 <= v < 1e-5 else "red")
    elif prop == "Pgp_Broccatelli":
        return "green" if v >= 0.5 else "red"
    elif prop == "Solubility_AqSolDB":
        return "green" if v >= -2 else ("yellow" if -4 <= v < -2 else "red")

    # ---- Distribution ----
    elif prop == "BBB_Martins":
        return "green" if v >= 0.5 else "red"
    elif prop == "PPBR_AZ":
        return "red" if v < 50 else ("yellow" if 50 <= v <= 90 else "green")
    elif prop == "VDss_Lombardo":
        return "red" if v < 0.71 else ("yellow" if 0.71 <= v <= 5 else "green")

    # ---- Metabolism ----
    elif prop in [
        "CYP1A2_Veith","CYP2C19_Veith","CYP2C9_Substrate_CarbonMangels",
        "CYP2C9_Veith","CYP2D6_Substrate_CarbonMangels","CYP2D6_Veith",
        "CYP3A4_Substrate_CarbonMangels","CYP3A4_Veith"
    ]:
        return "green" if v >= 0.5 else "red"

    # ---- Excretion ----
    elif prop == "Clearance_Hepatocyte_AZ":
        return "red" if v < 5 else ("yellow" if 5 <= v <= 30 else "green")
    elif prop == "Half_Life_Obach":
        return "green" if v >= 8 else ("yellow" if 2 <= v < 8 else "red")

    # ---- Toxicity ----
    elif prop in ["AMES","Carcinogens_Lagunin","ClinTox","DILI","hERG_Karim","Skin_Reaction","Tox21","Toxcast"]:
        return "red" if v >= 0.5 else "green"
    elif prop in ["hERG","herg_central"]:
        return "red" if v <= 1 else ("yellow" if 1 < v <= 10 else "green")
    elif prop == "LD50_Zhu":
        return "green" if v > 2000 else ("yellow" if 50 <= v <= 2000 else "red")

    return "gray"

# ---------- FLASK APP ----------
app = Flask(__name__)
CORS(app)

@app.route("/", methods=["GET"])
def home():
    return "Flask backend running!"

# ---------- PREDICT ----------
@app.route("/predict", methods=["POST"])
def predict():
    try:
        data = request.get_json(silent=True) or {}
        raw_smiles = data.get("smiles", [])
        if isinstance(raw_smiles, str):
            smiles_list = [raw_smiles.strip()]
        elif isinstance(raw_smiles, list):
            smiles_list = [s.strip() for s in raw_smiles if isinstance(s, str) and s.strip()]
        else:
            return jsonify({"molecules": []})

        molecules = []
        for smi in smiles_list:
            descriptors = compute_descriptors(smi) or {
                "MolWt": 0, "LogP": 0, "H-Donors": 0, "H-Acceptors": 0,
                "TPSA": 0, "RotBonds": 0, "Atoms": 0
            }

            mol_image = mol_to_image(smi) or ""
            radar_image = radar_plot_image(descriptors) or ""

            admet_results = {}
            for category, props in category_props.items():
                admet_results[category] = []
                for prop in props:
                    pred = safe_predict(prop, descriptors, smi, ADMET_MODELS)
                    value = "N/A"
                    if pred is not None:
                        try:
                            value = round(float(pred), 4)
                        except Exception:
                            value = str(pred)

                    status = categorize_property(prop, pred)
                    units = PROPERTY_UNITS.get(prop, "-")
                    druglikeness = [
                        {"scientist": s, "value": "Yes" if pred is not None and float(pred) > 0.5 else "No"}
                        for s in SCIENCE
                    ]

                    admet_results[category].append({
                        "property": prop,
                        "prediction": value,
                        "status": status,
                        "units": units,
                        "druglikeness": druglikeness
                    })

            molecules.append({
                "smiles": smi,
                "molImage": mol_image,
                "radarImage": radar_image,
                "descriptors": descriptors,
                "ADMET": admet_results
            })

        return jsonify({"molecules": molecules})

    except Exception as e:
        print(f"/predict unexpected error: {e}")
        return jsonify({"molecules": []}), 500

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000, debug=True)
