from .chem_utils import descriptors_to_array, compute_morgan_fingerprint
import numpy as np

def safe_predict(prop, descriptors, smiles, ADMET_MODELS):
    """
    Safely predict a property using the preloaded ADMET_MODELS.

    Arguments:
        prop: str, property name
        descriptors: dict, molecular descriptors
        smiles: str, SMILES string
        ADMET_MODELS: dict, preloaded models from app.py

    Returns:
        float or None: predicted value, None if prediction fails
    """
    try:
        # Find the model for this property
        model_dict = ADMET_MODELS.get(next((cat for cat, props in ADMET_MODELS.items() if prop in props), None), {}).get(prop)

        if not model_dict:
            return None  # model missing

        # Access actual sklearn model (handle if nested in dict)
        model = model_dict['model'] if isinstance(model_dict, dict) else model_dict

        n_feat = getattr(model, "n_features_in_", None)
        if n_feat is None:
            return None

        # Build feature array depending on model input size
        if n_feat == 7:
            features = descriptors_to_array(descriptors)
        elif n_feat == 2048:
            fp = compute_morgan_fingerprint(smiles)
            features = fp.tolist() if fp is not None else [0]*2048
        elif n_feat == 2055:
            fp = compute_morgan_fingerprint(smiles)
            desc_array = descriptors_to_array(descriptors)
            features = desc_array + (fp.tolist() if fp is not None else [0]*2048)
        else:
            # fallback
            features = descriptors_to_array(descriptors)

        # Ensure correct length
        if len(features) != n_feat:
            if len(features) < n_feat:
                features += [0]*(n_feat - len(features))
            else:
                features = features[:n_feat]

        prediction = float(model.predict([features])[0])
        return prediction

    except Exception as e:
        # Never crash, log optionally
        print(f"safe_predict error for {prop} / {smiles}: {e}")
        return None
