# utils/plotting.py
import matplotlib
matplotlib.use("Agg")  # Non-GUI backend
import matplotlib.pyplot as plt
import numpy as np
import io
import base64

def radar_plot_image(descriptors):
    """
    Convert descriptor values to a radar/spider plot PNG (base64 string)
    """
    labels = ["MolWt", "LogP", "H-Donors", "H-Acceptors", "TPSA", "RotBonds", "Atoms"]
    values = [
        descriptors.get("MolWt", 0)/500,
        (descriptors.get("LogP", 0)+5)/10,
        descriptors.get("H-Donors", 0)/5,
        descriptors.get("H-Acceptors", 0)/10,
        descriptors.get("TPSA", 0)/150,
        descriptors.get("RotBonds", 0)/10,
        descriptors.get("Atoms", 0)/50,
    ]

    # Complete loop
    angles = np.linspace(0, 2 * np.pi, len(labels), endpoint=False).tolist()
    values += values[:1]
    angles += angles[:1]

    # Increase figure size and DPI for clarity
    fig, ax = plt.subplots(figsize=(5,5), dpi=120, subplot_kw=dict(polar=True))
    ax.plot(angles, values, 'o-', linewidth=2)
    ax.fill(angles, values, alpha=0.25)
    ax.set_thetagrids(np.degrees(angles[:-1]), labels)
    ax.set_ylim(0,1)

    buf = io.BytesIO()
    plt.tight_layout()
    fig.savefig(buf, format="PNG", transparent=True)
    plt.close(fig)
    buf.seek(0)

    img_str = base64.b64encode(buf.read()).decode("utf-8")
    return f"data:image/png;base64,{img_str}"
