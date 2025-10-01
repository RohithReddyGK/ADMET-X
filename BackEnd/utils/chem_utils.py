from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
import numpy as np
from rdkit import DataStructs
from rdkit.Chem import Draw
import io
import base64

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

def descriptors_to_array(desc):
    return [
        desc["MolWt"], desc["LogP"], desc["H-Donors"], desc["H-Acceptors"],
        desc["TPSA"], desc["RotBonds"], desc["Atoms"]
    ]

def compute_morgan_fingerprint(smiles, radius=2, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    arr = np.zeros((n_bits,), dtype=int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

def mol_to_image(smiles, size=(300, 300)):
    """
    Generate a molecule image as base64 PNG if available,
    otherwise fall back to inline SVG.
    Always returns a string safe for frontend display.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        # placeholder if SMILES parsing fails
        mol = Chem.MolFromSmiles("C")

    try:
        # Try PNG rendering (requires Cairo/Pillow support)
        img = Draw.MolToImage(mol, size=size)
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        img_str = base64.b64encode(buf.getvalue()).decode("utf-8")
        return f"data:image/png;base64,{img_str}"
    except Exception:
        # Fallback to SVG rendering (always supported)
        svg = Draw.MolsToGridImage([mol], molsPerRow=1, subImgSize=size, useSVG=True)
        # Encode SVG as base64 for consistency with frontend <img src>
        svg_bytes = svg.encode("utf-8")
        svg_str = base64.b64encode(svg_bytes).decode("utf-8")
        return f"data:image/svg+xml;base64,{svg_str}"