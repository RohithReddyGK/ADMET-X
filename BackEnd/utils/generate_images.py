from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import base64
from io import BytesIO
from PIL import Image
import cairosvg 

# ---------------- Molecule Image ----------------
def mol_to_image(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        img = Image.new("RGB", (300, 300), (255, 255, 255))
        buffer = BytesIO()
        img.save(buffer, format="PNG")
        img_str = base64.b64encode(buffer.getvalue()).decode("utf-8")
        return f"data:image/png;base64,{img_str}"

    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    png_bytes = cairosvg.svg2png(bytestring=svg.encode("utf-8"))
    img_str = base64.b64encode(png_bytes).decode("utf-8")
    return f"data:image/png;base64,{img_str}"

# Fixed descriptor ranges for normalization
FIXED_RANGES = {
    "Atoms": (5, 60),
    "H-Acceptors": (0, 12),
    "H-Donors": (0, 10),
    "LogP": (-2, 7),
    "MolWt": (50, 800),
    "RotBonds": (0, 15),
    "TPSA": (0, 200),
}

# ---------------- Radar Plot ----------------
def radar_plot_multiple(molecule_descriptors_list, molecule_names=None, dark_mode=False):
    """
    Radar chart for single or multiple molecules using fixed descriptor ranges.
    Returns base64 PNG image.
    """
    # Ensure input is a list
    if isinstance(molecule_descriptors_list, dict):
        molecule_descriptors_list = [molecule_descriptors_list]

    if not molecule_descriptors_list:
        # Return blank image
        img = np.ones((300, 300, 3), dtype=np.uint8) * 255
        buffer = BytesIO()
        plt.imsave(buffer, img)
        img_str = base64.b64encode(buffer.getvalue()).decode("utf-8")
        return f"data:image/png;base64,{img_str}"

    try:
        # Use first molecule to get labels
        labels = list(molecule_descriptors_list[0].keys())
        num_vars = len(labels)

        # Fill missing descriptors with 0
        for desc in molecule_descriptors_list:
            for label in labels:
                if label not in desc or desc[label] is None:
                    desc[label] = 0

        angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
        angles += angles[:1]

        fig, ax = plt.subplots(figsize=(6,6), subplot_kw=dict(polar=True))

        if molecule_names is None:
            molecule_names = [f"Molecule {i+1}" for i in range(len(molecule_descriptors_list))]

        cmap = plt.cm.get_cmap('tab10')
        num_colors = len(molecule_descriptors_list)

        for idx, desc in enumerate(molecule_descriptors_list):
            values = []
            for label in labels:
                val = desc[label]
                min_val, max_val = FIXED_RANGES.get(label, (0, 1))
                norm_val = (val - min_val) / (max_val - min_val)
                norm_val = np.clip(norm_val, 0, 1)  # keep in [0,1]
                values.append(norm_val)
            values += values[:1]  # close loop

            color = cmap(idx / max(1, num_colors-1))
            ax.plot(angles, values, 'o-', linewidth=2, color=color, label=molecule_names[idx])
            ax.fill(angles, values, alpha=0.25, color=color)

        # Tick labels
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(labels, color='white' if dark_mode else 'black', fontsize=12)

        # Y-axis
        ax.set_yticks(np.linspace(0, 1, 5))
        ax.set_yticklabels([f"{round(y,2)}" for y in np.linspace(0,1,5)],
                           color='white' if dark_mode else 'black')
        ax.set_ylim(0, 1)
        ax.set_facecolor('black' if dark_mode else 'white')
        fig.patch.set_facecolor('black' if dark_mode else 'white')
        ax.legend(loc='upper right', bbox_to_anchor=(1.2,1.1))

        # Save to base64 PNG
        buffer = BytesIO()
        plt.savefig(buffer, format='PNG', bbox_inches='tight', transparent=True)
        plt.close(fig)
        img_str = base64.b64encode(buffer.getvalue()).decode("utf-8")
        return f"data:image/png;base64,{img_str}"

    except Exception as e:
        print(f"Radar plot error: {e}")
        img = np.ones((300,300,3), dtype=np.uint8) * 255
        buffer = BytesIO()
        plt.imsave(buffer, img)
        img_str = base64.b64encode(buffer.getvalue()).decode("utf-8")
        return f"data:image/png;base64,{img_str}"

# Keep old name for backward compatibility
radar_plot_image = radar_plot_multiple