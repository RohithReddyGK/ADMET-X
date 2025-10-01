import { useState } from "react";
import DrawMolecule from "./DrawMolecule";
import { motion } from "framer-motion";
import Lottie from "lottie-react";
import predictAnimation from "../assets/Predict.json"; 
import loadingAnimation from "../assets/Loading.json"; 
import PredictionPanel from "./PredictionPanel";

// Relaxed SMILES validation to allow real-world SMILES
const isValidSmiles = (smiles) =>
  /^[A-Za-z0-9@+\-\[\]\(\)=#$%./\\]+$/.test(smiles);

const Predict = () => {
  const [selectedOption, setSelectedOption] = useState("text");
  const [smilesInput, setSmilesInput] = useState("");
  const [file, setFile] = useState(null);
  const [results, setResults] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState("");

  // File upload handler
  const handleFileChange = (e) => {
    const uploadedFile = e.target.files[0];
    if (!uploadedFile) return;

    if (!["text/plain", "text/csv"].includes(uploadedFile.type)) {
      alert("Please upload a .txt or .csv file");
      return;
    }

    setFile(uploadedFile);
    const reader = new FileReader();
    reader.onload = (evt) => {
      const text = evt.target.result || "";
      const parsed = text.split("\n").map((s) => s.trim()).filter(Boolean);
      setSmilesInput(parsed.join("\n"));
    };
    reader.readAsText(uploadedFile);
  };

  // Example SMILES
  const handleExample = () => {
    setSmilesInput(
      `CC(=O)NC1=CC=C(C=C1)O Aspirin
CC(C)CC1=CC=C(C=C1)C(C)C Ibuprofen
C1CCCCC1 Cyclohexane
CCN(CC)CCO Triethylamine
C1=CC=CC=C1 Benzene`
    );
  };

  // Prediction handler
  const handlePredict = async () => {
    setError("");
    setResults(null);

    const processedInput = smilesInput
      .split("\n")
      .map((line) => line.trim())
      .filter(Boolean)
      .map((line) => {
        const [smiles, ...nameParts] = line.split(" ");
        return { smiles, name: nameParts.join(" ") };
      });

    if (!processedInput.length) {
      setError("Please provide at least one SMILES.");
      return;
    }

    const invalidSmiles = processedInput.filter((item) => !isValidSmiles(item.smiles));
    if (invalidSmiles.length > 0) {
      setError(`Invalid SMILES detected: ${invalidSmiles.map((i) => i.smiles).join(", ")}`);
      return;
    }

    try {
      setLoading(true); // START fullscreen loading

      const resp = await fetch("http://127.0.0.1:5000/predict", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ smiles: processedInput.map((i) => i.smiles) }),
      });

      if (!resp.ok) {
        const errData = await resp.json();
        throw new Error(errData.error || "Prediction failed");
      }

      const data = await resp.json();
      const molecules =
        Array.isArray(data)
          ? data.map((mol, idx) => ({ ...mol, name: processedInput[idx]?.name || "" }))
          : data?.molecules?.map((mol, idx) => ({ ...mol, name: processedInput[idx]?.name || "" })) || [];

      setResults({ molecules });
    } catch (err) {
      console.error("Prediction error:", err);
      setError(err.message || "Prediction failed. Check backend logs.");
    } finally {
      setLoading(false); // STOP fullscreen loading
    }
  };

  const downloadAllCSV = () => {
    if (!results?.molecules?.length) return;

    const csvRows = ["SMILES,Drug Name,Category,Property,Prediction,Units,Status,Druglikeness"];
    results.molecules.forEach((mol) => {
      const smi = mol.smiles;
      const name = mol.name || "";
      if (!mol.ADMET) {
        csvRows.push(`${smi},${name},,,,,,`);
        return;
      }
      Object.keys(mol.ADMET).forEach((cat) => {
        const catData = mol.ADMET[cat];
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
    });

    const csvContent = "data:text/csv;charset=utf-8," + csvRows.join("\n");
    const encodedUri = encodeURI(csvContent);
    const link = document.createElement("a");
    link.setAttribute("download", "smiles_predictions.csv");
    link.setAttribute("href", encodedUri);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  };

  return (
    <>
      {/* Fullscreen loading overlay */}
      {loading && (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-white dark:bg-gray-900 bg-opacity-80 backdrop-blur-sm">
          <Lottie animationData={loadingAnimation} loop className="w-72 h-72 sm:w-96 sm:h-96" />
        </div>
      )}

      <motion.div
        className={`p-4 sm:p-6 max-w-6xl mx-auto bg-white dark:bg-gray-800 rounded-lg shadow-2xl transition-all duration-300 ${
          loading ? "pointer-events-none opacity-40" : ""
        }`}
      >
        <h1 className="text-2xl sm:text-3xl font-bold mb-4 text-center text-gray-900 dark:text-white">
          Prediction of ADMET Properties
        </h1>

        <div className="flex flex-col md:flex-row gap-6">
          <div className="w-full md:w-1/2">
            {/* Existing decorative animation */}
            <Lottie animationData={predictAnimation} loop className="w-full h-64 sm:h-80 md:h-[400px]" />
          </div>

          <div className="w-full md:w-1/2">
            <div className="flex flex-wrap gap-2 justify-center md:justify-start mb-4">
              {["text", "file", "draw", "example"].map((option) => (
                <button
                  key={option}
                  className={`px-4 py-2 rounded font-medium ${
                    selectedOption === option ? "bg-blue-600 text-white" : "bg-blue-100 text-blue-800"
                  }`}
                  onClick={() => {
                    setSelectedOption(option);
                    if (option === "example") handleExample();
                  }}
                >
                  {option === "text"
                    ? "Text Input"
                    : option === "file"
                    ? "Upload File"
                    : option === "draw"
                    ? "Draw Molecule"
                    : "Example"}
                </button>
              ))}
            </div>

            {/* Input Sections */}
            {selectedOption === "text" && (
              <textarea
                className="w-full h-40 p-2 border rounded mb-4 text-gray-900 dark:text-white bg-white dark:bg-gray-800"
                placeholder="Paste SMILES here (optionally with drug name, e.g., SMILES DrugName)"
                value={smilesInput}
                onChange={(e) => setSmilesInput(e.target.value)}
              />
            )}

            {selectedOption === "file" && (
              <div className="mb-4">
                <label className="font-medium text-gray-900 dark:text-white">Upload file (.txt or .csv)</label>
                <input
                  type="file"
                  accept=".txt,.csv"
                  onChange={handleFileChange}
                  className="block border p-2 rounded mb-2 w-full"
                />
                {smilesInput && (
                  <textarea
                    readOnly
                    className="w-full h-40 p-2 border rounded text-gray-900 dark:text-white bg-white dark:bg-gray-800"
                    value={smilesInput}
                  />
                )}
              </div>
            )}

            {selectedOption === "draw" && (
              <div className="mb-4 border p-2 rounded">
                <DrawMolecule
                  onSmilesGenerated={(smi) =>
                    setSmilesInput((prev) => (prev ? prev + "\n" + smi : smi))
                  }
                />
              </div>
            )}

            {selectedOption === "example" && (
              <textarea
                readOnly
                className="w-full h-40 p-2 border rounded mb-4 text-gray-900 dark:text-white bg-white dark:bg-gray-800"
                value={smilesInput}
              />
            )}

            {/* Buttons */}
            <div className="flex flex-wrap gap-2 justify-center md:justify-start mt-2">
              <button
                onClick={handlePredict}
                className={`bg-green-600 text-white px-6 py-2 rounded-lg ${
                  loading ? "opacity-50 cursor-not-allowed" : ""
                }`}
                disabled={loading}
              >
                {loading ? "Predicting..." : "Predict"}
              </button>

              {results?.molecules?.length > 0 && (
                <button
                  onClick={downloadAllCSV}
                  className="bg-blue-600 text-white px-6 py-2 rounded-lg"
                >
                  Download All CSV
                </button>
              )}
            </div>

            {error && <p className="mt-2 text-red-600 font-medium">{error}</p>}
          </div>
        </div>

        {results && (
          <PredictionPanel
            results={results}
            smilesInput={smilesInput}
            setSmilesInput={setSmilesInput}
          />
        )}
      </motion.div>
    </>
  );
};

export default Predict;