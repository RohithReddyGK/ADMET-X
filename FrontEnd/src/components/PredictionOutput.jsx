import React from "react";
import absorptionIcon from "../assets/icons/Absorption.png";
import distributionIcon from "../assets/icons/Distribution.png";
import metabolismIcon from "../assets/icons/Metabolism.png";
import excretionIcon from "../assets/icons/Excretion.png";
import toxicityIcon from "../assets/icons/Toxicity.png";
import descriptorIcon from "../assets/icons/Molecular_Descriptors.png";
import { FiDownload } from "react-icons/fi";

const SCIENTISTS = ["Lipinski", "Ghose", "Veber", "Egan", "Muegge"];

const PredictionOutput = ({ result }) => {
  if (!result) return null;

  const { molImage, radarImage, descriptors, ADMET } = result;


  // ---------------- CSV Download ----------------
  const downloadReport = () => {
    const rows = [];

    // Title
    rows.push(["ADMET & Molecular Descriptors Report"]);
    rows.push([]);

    // Molecular Descriptors
    rows.push(["Molecular Descriptors"]);
    rows.push(["Descriptor", "Value"]);
    Object.entries(descriptors).forEach(([key, value]) => {
      rows.push([key, value]);
    });
    rows.push([]);

    // ADMET Properties
    Object.entries(ADMET).forEach(([category, props]) => {
      rows.push([`${category} Properties`]);
      rows.push(["Model", "Prediction", ...SCIENTISTS, "Units"]);
      props.forEach((model) => {
        const row = [
          model.property,
          model.prediction?.toFixed(3) ?? "-",
          ...model.druglikeness.map((d) => d.value),
          model.units,
        ];
        rows.push(row);
      });
      rows.push([]);
    });

    // Convert to CSV string
    const csvContent = rows.map((r) => r.join(",")).join("\n");

    // Trigger download
    const blob = new Blob([csvContent], { type: "text/csv" });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `ADMET_Report_${result.smiles}.csv`;
    document.body.appendChild(a);
    a.click();
    a.remove();
    window.URL.revokeObjectURL(url);
  };

  // Helper to render descriptors table
  const renderDescriptors = () => (
    <div className="overflow-x-auto">
      <table className="table-auto w-full border border-blue-200 rounded-lg overflow-hidden">
        <thead>
          <tr className="bg-blue-100 dark:bg-blue-800">
            <th className="border px-3 py-1 text-left">Descriptor</th>
            <th className="border px-3 py-1 text-left">Value</th>
          </tr>
        </thead>
        <tbody>
          {Object.entries(descriptors).map(([key, value]) => (
            <tr
              key={key}
              className="odd:bg-white even:bg-blue-50 dark:odd:bg-gray-700 dark:even:bg-gray-800"
            >
              <td className="border px-3 py-1">{key}</td>
              <td className="border px-3 py-1">{value}</td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );

  // Helper to render property tables (Absorption, Distribution, etc.)
  const renderPropertyTable = (propName) => {
    const propData = ADMET[propName];
    if (!propData) return null;

    return (
      <div className="overflow-x-auto">
        <table className="table-auto w-full border border-blue-200 rounded-lg overflow-hidden text-sm">
          <thead>
            <tr className="bg-blue-100 dark:bg-blue-800">
              <th className="border px-2 py-1">Model</th>
              <th className="border px-2 py-1">Prediction</th>
              {SCIENTISTS.map((sci) => (
                <th key={sci} className="border px-2 py-1 text-center">{sci}</th>
              ))}
              <th className="border px-2 py-1">Units</th>
            </tr>
          </thead>
          <tbody>
            {propData.map((model) => (
              <tr
                key={model.property}
                className="odd:bg-white even:bg-blue-50 dark:odd:bg-gray-700 dark:even:bg-gray-800"
              >
                <td className="border px-2 py-1">{model.property}</td>
                <td className="border px-2 py-1">{model.prediction?.toFixed(3) || "-"}</td>
                {model.druglikeness.map((d) => (
                  <td key={d.scientist} className="border px-2 py-1 text-center">{d.value}</td>
                ))}
                <td className="border px-2 py-1">{model.units}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    );
  };

  return (
    <div className="max-w-7xl mx-auto my-6">
      {/* Heading with Download Button */}
      <div className="relative mb-6">
        <h1 className="text-3xl font-bold mb-6 text-center text-black-600 dark:text-blue-600">
          ADMET & Descriptors Report
        </h1>
        <button
          onClick={downloadReport}
          className="absolute right-0 top-1/2 transform -translate-y-1/2 bg-blue-500 hover:bg-blue-600 text-white px-4 py-2 rounded shadow-md ml-4 flex items-center gap-2"
        >
          <FiDownload
            size={20}
            className="transition-transform duration-200 group-hover:translate-x-1 text-white group-hover:text-yellow-300"
          />
          Download Report
        </button>
      </div>

      {/* Grid for 4x2 layout */}
      <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
        {/* Row 1: Molecule and Radar */}
        <div className="flex flex-col justify-center items-center border border-blue-200 rounded-xl shadow-md p-2 bg-white dark:bg-gray-800">
          <h3 className="text-2xl font-bold mb-2 text-blue-600 dark:text-blue-400 flex items-center justify-center gap-2">
            Molecular Structure
          </h3>
          <img src={molImage} alt="Molecule" className="max-h-80 rounded-lg mb-3" />
          <h1 className="text-xl font-semibold text-gray-800 dark:text-gray-300 break-words text-center">
            SMILES: {result?.smiles}</h1>
        </div>
        <div className="flex justify-center items-center border border-blue-200 rounded-xl shadow-md p-2 bg-white dark:bg-gray-800">
          <h3 className="text-2xl font-bold mb-2 text-blue-600 dark:text-blue-400 flex items-center justify-center gap-2">
            Radar Plot
          </h3>
          <img src={radarImage} alt="Radar" className="max-h-80 rounded-lg" />
        </div>

        {/* Row 2: Descriptors and Absorption */}
        <div className="border border-blue-200 rounded-xl shadow-md p-2 bg-white dark:bg-gray-800">
          <h3 className="text-2xl font-bold mb-2 text-blue-600 dark:text-blue-400 flex items-center justify-center gap-2">
            <img src={descriptorIcon} alt="Descriptors" className="w-20 h-20" />
            Molecular Descriptors
          </h3>
          {renderDescriptors()}
        </div>
        <div className="border border-blue-200 rounded-xl shadow-md p-2 bg-white dark:bg-gray-800">
          <h3 className="text-2xl font-bold mb-2 text-blue-600 dark:text-blue-400 flex items-center justify-center gap-2">
            <img src={absorptionIcon} alt="Absorption" className="w-20 h-20" />
            Absorption
          </h3>
          {renderPropertyTable("Absorption")}
        </div>

        {/* Row 3: Distribution and Metabolism */}
        <div className="border border-blue-200 rounded-xl shadow-md p-2 bg-white dark:bg-gray-800">
          <h3 className="text-2xl font-bold mb-2 text-blue-600 dark:text-blue-400 flex items-center justify-center gap-2">
            <img src={distributionIcon} alt="Distribution" className="w-20 h-20" />
            Distribution
          </h3>
          {renderPropertyTable("Distribution")}
        </div>
        <div className="border border-blue-200 rounded-xl shadow-md p-2 bg-white dark:bg-gray-800">
          <h3 className="text-2xl font-bold mb-2 text-blue-600 dark:text-blue-400 flex items-center justify-center gap-2">
            <img src={metabolismIcon} alt="Metabolism" className="w-20 h-20" />
            Metabolism
          </h3>
          {renderPropertyTable("Metabolism")}
        </div>

        {/* Row 4: Excretion and Toxicity */}
        <div className="border border-blue-200 rounded-xl shadow-md p-2 bg-white dark:bg-gray-800">
          <h3 className="text-2xl font-bold mb-2 text-blue-600 dark:text-blue-400 flex items-center justify-center gap-2">
            <img src={excretionIcon} alt="Excretion" className="w-20 h-20" />
            Excretion
          </h3>
          {renderPropertyTable("Excretion")}
        </div>
        <div className="border border-blue-200 rounded-xl shadow-md p-2 bg-white dark:bg-gray-800">
          <h3 className="text-2xl font-bold mb-2 text-blue-600 dark:text-blue-400 flex items-center justify-center gap-2">
            <img src={toxicityIcon} alt="Toxicity" className="w-20 h-20" />
            Toxicity
          </h3>
          {renderPropertyTable("Toxicity")}
        </div>
      </div>
    </div>
  );
};

export default PredictionOutput;