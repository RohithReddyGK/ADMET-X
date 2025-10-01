import React, { useState, useEffect } from "react";
import { Collapse, Button, message, Tag } from "antd";
import Plot from "react-plotly.js";
import Conclusion from "./Conclusion";

import absorptionIcon from "../assets/icons/Absorption.png";
import distributionIcon from "../assets/icons/Distribution.png";
import metabolismIcon from "../assets/icons/Metabolism.png";
import excretionIcon from "../assets/icons/Excretion.png";
import toxicityIcon from "../assets/icons/Toxicity.png";

const { Panel } = Collapse;

const PredictionPanel = ({ results, smilesInput, setSmilesInput }) => {
  const initialMolecules = results?.molecules || [];
  const [molecules, setMolecules] = useState(initialMolecules);

  useEffect(() => {
    setMolecules(results?.molecules || []);
  }, [results]);

  const handleCopy = (smi, name) => {
    const textToCopy = name ? `${smi} ${name}` : smi;
    navigator.clipboard.writeText(textToCopy);
    message.success("SMILES copied!");
  };

  const handleRemove = (idx) => {
    const updated = [...molecules];
    const removedMol = updated.splice(idx, 1)[0];
    setMolecules(updated);

    if (setSmilesInput && removedMol) {
      const inputLines = smilesInput.split("\n");
      inputLines.splice(idx, 1);
      setSmilesInput(inputLines.join("\n"));
    }
  };

  const downloadCSV = (molData) => {
    if (!molData) return;
    const smi = molData.smiles || "";
    const name = molData.name || "";

    const csvRows = ["SMILES,Drug Name,Category,Property,Prediction,Units,Status,Druglikeness"];

    if (molData.ADMET) {
      Object.keys(molData.ADMET).forEach((cat) => {
        const catData = molData.ADMET[cat];
        if (Array.isArray(catData) && catData.length > 0) {
          catData.forEach((row) => {
            const prediction = row.prediction ?? "";
            const units = row.units || "";
            const status = row.status || "";
            const druglikeness = Array.isArray(row.druglikeness)
              ? row.druglikeness.map((d) => `${d.scientist}: ${d.value ?? "N/A"}`).join("; ")
              : "";
            csvRows.push(`${smi},${name},${cat},${row.property},${prediction},${units},${status},${druglikeness}`);
          });
        }
      });
    }

    const csvContent = "data:text/csv;charset=utf-8," + csvRows.join("\n");
    const encodedUri = encodeURI(csvContent);
    const link = document.createElement("a");
    link.setAttribute("download", `${smi.replace(/[^a-z0-9]/gi, "_")}_prediction.csv`);
    link.setAttribute("href", encodedUri);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  };

  if (!molecules.length) {
    return (
      <div className="mt-10 p-6 bg-gray-100 dark:bg-gray-900 rounded-lg shadow">
        <h2 className="text-2xl font-bold mb-4 text-gray-900 dark:text-white">Prediction Results</h2>
        <p className="text-gray-700 dark:text-gray-300">No molecules to display.</p>
      </div>
    );
  }

  const statusTag = (status) => {
    switch (status) {
      case "green": return <Tag color="green">Good</Tag>;
      case "yellow": return <Tag color="gold">Moderate</Tag>;
      case "red": return <Tag color="red">Poor</Tag>;
      default: return <Tag color="default">N/A</Tag>;
    }
  };

  const getCategoryIcon = (cat) => {
    switch (cat.toLowerCase()) {
      case "absorption": return absorptionIcon;
      case "distribution": return distributionIcon;
      case "metabolism": return metabolismIcon;
      case "excretion": return excretionIcon;
      case "toxicity": return toxicityIcon;
      default: return null;
    }
  };

  return (
    <div className="mt-10 p-4 md:p-6 bg-gray-100 dark:bg-gray-900 rounded-lg shadow overflow-x-hidden">
      <h2 className="text-2xl font-bold mb-4 text-gray-900 dark:text-white">Prediction Results</h2>
      <Collapse accordion={false}>
        {molecules.map((mol, idx) => {
          const smi = mol.smiles || smilesInput.split("\n")[idx] || "";
          const name = mol.name || "";

          return (
            <Panel
              header={
                <div className="flex flex-col sm:flex-row sm:justify-between w-full gap-2">
                  {/* SMILES & Name */}
                  <div
                    className="break-words text-gray-900 dark:text-white"
                    style={{ wordBreak: "break-word", overflowWrap: "anywhere" }}
                    title={`${smi}${name ? ` (${name})` : ""}`}
                  >
                    {smi}{name ? ` (${name})` : ""}
                  </div>

                  {/* Buttons */}
                  <div className="flex flex-wrap gap-2 mt-2 sm:mt-0">
                    <Button size="small" onClick={(e) => { e.stopPropagation(); handleCopy(smi, name); }}>
                      Copy Smile Name
                    </Button>
                    <Button size="small" danger onClick={(e) => { e.stopPropagation(); handleRemove(idx); }}>
                      Remove
                    </Button>
                    <Button size="small" type="primary" onClick={(e) => { e.stopPropagation(); downloadCSV(mol); }}>
                      Download CSV
                    </Button>
                  </div>
                </div>
              }
              key={idx}
            >
              {/* Images + Radar Plot Side by Side */}
              <div className="flex flex-col md:flex-row gap-4 mb-6 items-center">
                {mol.molImage && (
                  <img
                    src={mol.molImage}
                    alt={smi}
                    className="w-full md:w-1/2 rounded object-contain"
                    style={{ maxHeight: "300px" }}
                  />
                )}
                {mol.radarImage && (
                  <img
                    src={mol.radarImage}
                    alt="radar"
                    className="w-full md:w-1/2 rounded object-contain"
                    style={{ maxHeight: "300px" }}
                  />
                )}
              </div>

              {/* Spider Plot */}
              {mol.spider_data?.labels?.length > 0 && (
                <div className="mb-4">
                  <Plot
                    data={[{
                      type: "scatterpolar",
                      r: mol.spider_data.values,
                      theta: mol.spider_data.labels,
                      fill: "toself",
                    }]}
                    layout={{
                      polar: { radialaxis: { visible: true, range: [0, 1] } },
                      showlegend: false,
                      margin: { t: 0, b: 0, l: 0, r: 0 },
                      autosize: true,
                    }}
                    style={{ width: "100%", height: 400 }}
                    useResizeHandler={true}
                  />
                </div>
              )}

              {/* ADMET Tables */}
              <div className="mt-4 space-y-4">
                {mol.ADMET && Object.keys(mol.ADMET).map((cat) => (
                  <div key={cat} className="overflow-x-auto">
                    <h4 className="text-xl font-bold text-blue-600 bg-white dark:bg-gray-800 px-2 py-1 rounded flex items-center gap-2">
                      {getCategoryIcon(cat) && (
                        <img src={getCategoryIcon(cat)} alt={cat} className="w-10 h-10 inline-block" />
                      )}
                      {cat}
                    </h4>
                    {Array.isArray(mol.ADMET[cat]) && mol.ADMET[cat].length > 0 ? (
                      <table className="min-w-full border text-sm text-gray-900 dark:text-white bg-white dark:bg-gray-800">
                        <thead>
                          <tr className="bg-blue-200 dark:bg-gray-700">
                            <th className="px-2 py-1 border">Property</th>
                            <th className="px-2 py-1 border">Prediction</th>
                            <th className="px-2 py-1 border">Units</th>
                            <th className="px-2 py-1 border">Status</th>
                            <th className="px-2 py-1 border">Druglikeness</th>
                          </tr>
                        </thead>
                        <tbody>
                          {mol.ADMET[cat].map((row, i) => (
                            <tr key={i} className={i % 2 === 0 ? "bg-blue-50 dark:bg-gray-900" : "bg-white dark:bg-gray-800"}>
                              <td className="px-2 py-1 border">{row.property}</td>
                              <td className="px-2 py-1 border">{row.prediction?.toFixed?.(3) ?? "N/A"}</td>
                              <td className="px-2 py-1 border">{row.units || "-"}</td>
                              <td className="px-2 py-1 border">{statusTag(row.status)}</td>
                              <td className="px-2 py-1 border">
                                {Array.isArray(row.druglikeness)
                                  ? row.druglikeness.map(d => `${d.scientist}: ${d.value ?? "N/A"}`).join("; ")
                                  : "-"}
                              </td>
                            </tr>
                          ))}
                        </tbody>
                      </table>
                    ) : (
                      <div className="italic text-gray-600 dark:text-gray-300 p-2">No data</div>
                    )}
                  </div>
                ))}
              </div>
              <Conclusion admet={mol.ADMET} />
            </Panel>
          );
        })}
      </Collapse>
    </div>
  );
};

export default PredictionPanel;
