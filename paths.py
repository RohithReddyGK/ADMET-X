import os

# ---------------- Project Root ----------------
PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))

# ---------------- Dataset paths ----------------
def train_path(category, model_name):
    return os.path.join(PROJECT_ROOT, "admet_data", category, model_name, "train.csv")

def valid_path(category, model_name):
    return os.path.join(PROJECT_ROOT, "admet_data", category, model_name, "valid.csv")

def test_path(category, model_name):
    return os.path.join(PROJECT_ROOT, "admet_data", category, model_name, "test.csv")

# ---------------- Model & Prediction paths ----------------
def model_path(category, model_name):
    return os.path.join(PROJECT_ROOT, "BackEnd", "Models", category, model_name + ".joblib")

def prediction_path(category, model_name):
    return os.path.join(PROJECT_ROOT, "Model_predictions", category, model_name + ".csv")
