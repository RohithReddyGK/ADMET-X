import React, { useState, useEffect } from "react";
import { Collapse, Button, message, Tag } from "antd";
import Plot from "react-plotly.js";
import Conclusion from "./Conclusion";
import TooltipCell from "./ToolTipCell";

// Import your icons
import molecularDescriptorIcon from "../assets/icons/Molecular_Descriptors.png";
import absorptionIcon from "../assets/icons/Absorption.png";
import distributionIcon from "../assets/icons/Distribution.png";
import metabolismIcon from "../assets/icons/Metabolism.png";
import excretionIcon from "../assets/icons/Excretion.png";
import toxicityIcon from "../assets/icons/Toxicity.png";

const { Panel } = Collapse;

const categoryIcons = {
  Absorption: absorptionIcon,
  Distribution: distributionIcon,
  Metabolism: metabolismIcon,
  Excretion: excretionIcon,
  Toxicity: toxicityIcon,
};

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

    // Add Molecular Descriptors
    if (molData.descriptors) {
      Object.entries(molData.descriptors).forEach(([desc, value]) => {
        csvRows.push(`${smi},${name},Molecular Descriptor,${desc},${value},-,-,-`);
      });
    }

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
                  <div
                    className="break-words text-gray-900 dark:text-white"
                    style={{ wordBreak: "break-word", overflowWrap: "anywhere" }}
                    title={`${smi}${name ? ` (${name})` : ""}`}
                  >
                    {smi}{name ? ` (${name})` : ""}
                  </div>
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
              {/* Images + Radar Plot */}
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

              {/* Molecular Descriptors Table */}
              {mol.descriptors && Object.keys(mol.descriptors).length > 0 && (
                <div className="mt-4 overflow-x-auto">
                  <h4 className="flex items-center gap-2 text-2xl font-bold text-gray-800 dark:text-white bg-white dark:bg-gray-600 px-2 py-1 rounded">
                    <img src={molecularDescriptorIcon} alt="Molecular Descriptors" className="w-12 h-12" />
                    Molecular Descriptors
                  </h4>
                  <table className="min-w-full border text-sm text-gray-900 dark:text-gray-100 bg-blue-50 dark:bg-gray-900 rounded-lg overflow-hidden">
                    <thead>
                      <tr className="bg-blue-200 dark:bg-gray-800 text-gray-900 dark:text-gray-100 font-semibold">
                        <th className="px-2 py-1 border border-blue-300 dark:border-gray-700">Descriptor</th>
                        <th className="px-2 py-1 border border-blue-300 dark:border-gray-700">Value</th>
                      </tr>
                    </thead>
                    <tbody>
                      {Object.entries(mol.descriptors).map(([desc, value], i) => (
                        <tr
                          key={i}
                          className={`${i % 2 === 0
                            ? "bg-blue-100 dark:bg-gray-800"
                            : "bg-blue-50 dark:bg-gray-700"
                            }`}
                        >
                          <td className="px-2 py-1 border border-blue-300 dark:border-gray-700">{desc}</td>
                          <td className="px-2 py-1 border border-blue-300 dark:border-gray-700">{value}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              )}

              {/* ADMET Tables with Icons */}
              <div className="mt-4 space-y-4">
                {mol.ADMET && Object.keys(mol.ADMET).map((cat) => (
                  <div key={cat} className="overflow-x-auto relative z-10 overflow-visible">
                    <h4 className="flex items-center gap-2 text-2xl font-bold text-gray-800 dark:text-white bg-white dark:bg-gray-600 px-2 py-1 rounded flex item-center gap 2">
                      {categoryIcons[cat] && (
                        <img src={categoryIcons[cat]} alt={cat} className="w-10 h-10" />
                      )}
                      {cat}
                    </h4>
                    {Array.isArray(mol.ADMET[cat]) && mol.ADMET[cat].length > 0 ? (
                      <table className="min-w-full border text-sm text-gray-900 dark:text-gray-100 bg-blue-50 dark:bg-gray-900 rounded-lg overflow-hidden">
                        <thead>
                          <tr className="bg-blue-200 dark:bg-gray-800 text-gray-900 dark:text-gray-100 font-semibold">
                            <th className="px-2 py-1 border border-blue-300 dark:border-gray-700">Property</th>
                            <th className="px-2 py-1 border border-blue-300 dark:border-gray-700">Prediction</th>
                            <th className="px-2 py-1 border border-blue-300 dark:border-gray-700">Units</th>
                            <th className="px-2 py-1 border border-blue-300 dark:border-gray-700">Status</th>
                            <th className="px-2 py-1 border border-blue-300 dark:border-gray-700">Druglikeness</th>
                          </tr>
                        </thead>
                        <tbody>
                          {mol.ADMET[cat].map((row, i) => (
                            <tr
                              key={i}
                              className={`${i % 2 === 0
                                ? "bg-blue-100 dark:bg-gray-800" // softened for even rows
                                : "bg-blue-50 dark:bg-gray-700"  // softened for odd rows
                                }`}
                            >
                              <td className="px-2 py-1 border border-blue-300 dark:border-gray-700">
                                <TooltipCell name={row.property} />
                              </td>
                              <td className="px-2 py-1 border border-blue-300 dark:border-gray-700">
                                {row.prediction?.toFixed?.(3) ?? "N/A"}
                              </td>
                              <td className="px-2 py-1 border border-blue-300 dark:border-gray-700">{row.units || "-"}</td>
                              <td className="px-2 py-1 border border-blue-300 dark:border-gray-700">{statusTag(row.status)}</td>
                              <td className="px-2 py-1 border border-blue-300 dark:border-gray-700">
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
              <Conclusion admet={mol.ADMET} smiles={smi} drugName={name} />
            </Panel>
          );
        })}
      </Collapse>
    </div>
  );
};

export default PredictionPanel;