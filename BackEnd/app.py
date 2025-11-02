from flask import Flask, request, jsonify
from flask_cors import CORS
import os
import sys
import joblib
import traceback

# Deployment-safe paths 
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Ensure modules are importable
sys.path.append(BASE_DIR)
sys.path.append(os.path.join(BASE_DIR, "descriptors"))
sys.path.append(os.path.join(BASE_DIR, "utils"))

from utils.chem_utils import compute_descriptors, mol_to_image
from utils.plotting import radar_plot_image
from utils.prediction_utils import safe_predict

# CONFIG
MODELS_DIR = os.path.join(os.getcwd(), "Models")
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

# LOAD MODELS
ADMET_MODELS = {}
print("[INFO] Starting to load ADMET models...")
for category, props in category_props.items():
    category_path = os.path.join(MODELS_DIR, category)
    ADMET_MODELS[category] = {}
    if os.path.isdir(category_path):
        for model_file in os.listdir(category_path):
            if model_file.endswith(".joblib"):
                model_name = model_file.replace(".joblib", "")
                model_path = os.path.join(category_path, model_file)
                try:
                    ADMET_MODELS[category][model_name] = joblib.load(model_path)
                    print(f"✅ Loaded {category}/{model_file}")
                except Exception as e:
                    print(f"⚠️ Failed to load {category}/{model_file}: {e}")
                    traceback.print_exc()
    else:
        print(f"⚠️ Category folder not found: {category_path}")
print(f"[INFO] Model loading finished. Total categories: {len(ADMET_MODELS)}")

# PROPERTY STATUS
def categorize_property(prop, value):
    try:
        v = float(value)
    except:
        return "gray"
    
    # Absorption
    if prop == "Bioavailability_Ma": return "green" if v>=0.5 else "red"
    if prop == "Caco2_Wang": return "green" if v>-5.15 else ("yellow" if -6<=v<=-5.15 else "red")
    if prop == "HIA_Hou": return "green" if v>=0.5 else "red"
    if prop == "HydrationFreeEnergy_FreeSolv": return "green" if v<-5 else ("yellow" if -5<=v<=0 else "red")
    if prop == "Lipophilicity_AstraZeneca": return "green" if -0.7<=v<=5.0 else "red"
    if prop == "PAMPA_NCATS": return "green" if v>=1e-5 else ("yellow" if 1e-6<=v<1e-5 else "red")
    if prop == "Pgp_Broccatelli": return "green" if v>=0.5 else "red"
    if prop == "Solubility_AqSolDB": return "green" if v>=-2 else ("yellow" if -4<=v<-2 else "red")
    
    # Distribution
    if prop == "BBB_Martins": return "green" if v>=0.5 else "red"
    if prop == "PPBR_AZ": return "red" if v<50 else ("yellow" if 50<=v<=90 else "green")
    if prop == "VDss_Lombardo": return "red" if v<0.71 else ("yellow" if 0.71<=v<=5 else "green")
    
    # Metabolism
    if prop in ["CYP1A2_Veith","CYP2C19_Veith","CYP2C9_Substrate_CarbonMangels","CYP2C9_Veith",
                "CYP2D6_Substrate_CarbonMangels","CYP2D6_Veith","CYP3A4_Substrate_CarbonMangels","CYP3A4_Veith"]:
        return "green" if v>=0.5 else "red"
    
    # Excretion
    if prop=="Clearance_Hepatocyte_AZ": return "red" if v<5 else ("yellow" if 5<=v<=30 else "green")
    if prop=="Half_Life_Obach": return "green" if v>=8 else ("yellow" if 2<=v<8 else "red")
    
    # Toxicity
    if prop in ["AMES","Carcinogens_Lagunin","ClinTox","DILI","hERG_Karim","Skin_Reaction","Tox21","Toxcast"]:
        return "red" if v>=0.5 else "green"
    if prop in ["hERG","herg_central"]: return "red" if v<=1 else ("yellow" if 1<v<=10 else "green")
    if prop=="LD50_Zhu": return "green" if v>2000 else ("yellow" if 50<=v<=2000 else "red")
    
    return "gray"

# FLASK APP
app = Flask(__name__)
CORS(app, resources={r"/*": {"origins":"*"}})

@app.route("/", methods=["GET"])
def status():
    return jsonify({"status":"Backend is running"}), 200

@app.route("/predict", methods=["POST"])
def predict():
    try:
        data = request.get_json(force=True) or {}
        raw_smiles = data.get("smiles", [])
        if isinstance(raw_smiles,str): smiles_list=[raw_smiles.strip()]
        elif isinstance(raw_smiles,list): smiles_list=[s.strip() for s in raw_smiles if s.strip()]
        else: return jsonify({"molecules":[]})
        
        molecules=[]
        for smi in smiles_list:
            # Descriptors
            descriptors = compute_descriptors(smi) or {
                "MolWt":0,"LogP":0,"H-Donors":0,"H-Acceptors":0,
                "TPSA":0,"RotBonds":0,"Atoms":0
            }
            mol_img = mol_to_image(smi) or ""
            radar_img = radar_plot_image(descriptors) or ""

            # ADMET Predictions 
            admet_results = {}
            for category, props in category_props.items():
                admet_results[category] = []
                for prop in props:
                    try:
                        pred = safe_predict(prop, descriptors, smi, ADMET_MODELS)
                        if pred is None:
                            print(f"[WARN] Prediction for {prop} returned None")
                    except Exception as e:
                        print(f"[ERROR] safe_predict failed for {prop}: {e}")
                        traceback.print_exc()
                        pred = None

                    value = "N/A"
                    if pred is not None:
                        try:
                            value = round(float(pred), 4)
                        except:
                            value = str(pred)
                    status_val = categorize_property(prop, pred)
                    units = PROPERTY_UNITS.get(prop, "-")
                    druglikeness = [
                        {"scientist": s, "value": "Yes" if pred is not None and float(pred) > 0.5 else "No"}
                        for s in SCIENCE
                    ]
                    admet_results[category].append({
                        "property": prop,
                        "prediction": value,
                        "status": status_val,
                        "units": units,
                        "druglikeness": druglikeness
                    })

            molecules.append({
                "smiles": smi,
                "molImage": mol_img,
                "radarImage": radar_img,
                "descriptors": descriptors,
                "ADMET": admet_results
            })

        return jsonify({"molecules": molecules})

    except Exception as e:
        print(f"[ERROR] /predict failed: {e}")
        traceback.print_exc()
        return jsonify({"molecules": []}), 500

# START SERVER
if __name__ == "__main__":
    print("Starting Flask server...")
    port = int(os.environ.get("PORT", 8080))
    app.run(host="0.0.0.0", port=port, debug=True, threaded=True)