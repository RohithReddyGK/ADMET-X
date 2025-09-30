import { useState } from "react";
import DrawMolecule from "./DrawMolecule";
import { motion } from "framer-motion";
import Lottie from "lottie-react";
import predictAnimation from "../assets/Predict.json"; // Your new animation
import loadingAnimation from "../assets/Loading.json"; //Loading animation

const Predict = ({ setResult }) => {
  const [selectedOption, setSelectedOption] = useState("text");
  const [smilesInput, setSmilesInput] = useState("");
  const [file, setFile] = useState(null);
  const [loading, setLoading] = useState(false);

  const handleFileChange = (e) => {
    const uploadedFile = e.target.files[0];
    if (uploadedFile) {
      const fileType = uploadedFile.type;
      if (fileType === "text/plain" || fileType === "text/csv") {
        setFile(uploadedFile);
        const reader = new FileReader();
        reader.onload = (event) => {
          const text = event.target.result;
          if (fileType === "text/csv") {
            const parsedSMILES = text.split("\n").map(line => line.trim());
            setSmilesInput(parsedSMILES.join("\n"));
          } else setSmilesInput(text);
        };
        reader.readAsText(uploadedFile);
      } else alert("Please upload a valid .txt or .csv file.");
    }
  };

  const handleExample = () => setSmilesInput(`CC(=O)NC1=CC=C(C=C1)O`);

  const handlePredict = async () => {
    if (!smilesInput.trim()) return alert("Please provide at least one SMILE input.");

    setLoading(true); //Start Full-screen overlay.
    setResult(null);

    console.log("Predicting for SMILES:\n", smilesInput);

    try {
      // Send SMILES to Flask backend
      const response = await fetch("http://127.0.0.1:5000/predict", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({ smiles: smilesInput }),
      });

      if (!response.ok) throw new Error("Prediction request failed");

      const data = await response.json();
      setResult(data); // Send the backend JSON response to PredictionOutput.
      console.log("Prediction Result:", data);

    } catch (error) {
      console.error(error);
      alert("Error while predicting. Make sure the backend is running.");
    } finally {
      setLoading(false); //Stop Spinner.
    }
  };

  return (
    <motion.div
      className="p-6 max-w-6xl mx-auto bg-white dark:bg-gray-800 rounded-lg shadow-2xl animate-slide-up"
      initial={{ opacity: 0, y: 30 }}
      whileInView={{ opacity: 1, y: 0 }}
      viewport={{ once: true }}
      transition={{ duration: 0.6 }}
    >
      {/* 🔹 Fullscreen Loading Overlay */}
      {loading && (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-black bg-opacity-60">
          <Lottie animationData={loadingAnimation} loop={true} className="w-48 h-48" />
        </div>
      )}

      <h1 className="text-3xl font-bold mb-6 text-center">Prediction of ADMET Properties</h1>

      <div className="flex flex-col md:flex-row items-center md:items-start gap-8">
        {/* Left: Lottie Animation */}
        <motion.div
          className="w-full md:w-1/2"
          initial={{ opacity: 0, x: -50 }}
          whileInView={{ opacity: 1, x: 0 }}
          viewport={{ once: true }}
          transition={{ duration: 0.8 }}
        >
          <Lottie animationData={predictAnimation} loop={true} className="w-full h-80 md:h-[400px]" />
        </motion.div>

        {/* Right: Input Tabs + Content */}
        <motion.div
          className="w-full md:w-1/2"
          initial={{ opacity: 0, x: 50 }}
          whileInView={{ opacity: 1, x: 0 }}
          viewport={{ once: true }}
          transition={{ duration: 0.8 }}
        >
          {/* Tabs */}
          <div className="flex justify-center space-x-4 mb-4 md:justify-start">
            {["text", "file", "draw", "example"].map(option => (
              <button
                key={option}
                className={`px-4 py-2 rounded font-medium transition ${selectedOption === option ? "bg-blue-600 text-white dark:bg-blue-800" : "bg-blue-100 text-blue-800 dark:bg-gray-700 dark:text-white"}`}
                onClick={() => { setSelectedOption(option); if (option === "example") handleExample(); }}
              >
                {option === "text" ? "Text Input" : option === "file" ? "Upload File" : option === "draw" ? "Draw Molecule" : "Example"}
              </button>
            ))}
          </div>

          {/* Content */}
          {selectedOption === "text" && (
            <textarea
              className="w-full h-40 p-2 border rounded-lg focus:ring-2 focus:ring-blue-500 dark:bg-gray-700 dark:text-white dark:border-gray-600 mb-4"
              placeholder="Paste SMILES here (one per line)..."
              value={smilesInput}
              onChange={e => setSmilesInput(e.target.value)}
            />
          )}

          {selectedOption === "file" && (
            <div className="space-y-3 mb-4">
              <label className="font-medium">Upload a .txt or .csv file containing SMILES</label>
              <input
                type="file"
                accept=".txt,.csv"
                onChange={handleFileChange}
                className="block border p-2 rounded-lg dark:bg-gray-700 dark:text-white dark:border-gray-600"
              />
              {smilesInput && (
                <textarea
                  className="w-full h-40 p-2 border rounded-lg focus:ring-2 focus:ring-blue-500"
                  value={smilesInput}
                  readOnly
                />
              )}
            </div>
          )}

          {selectedOption === "draw" && (
            <div className="w-full p-2 border rounded-lg mb-4">
              <DrawMolecule onSmilesGenerated={(smiles) => setSmilesInput(smiles)} />
            </div>
          )}

          {selectedOption === "example" && (
            <textarea
              className="w-full h-40 p-2 border rounded-lg focus:ring-2 focus:ring-blue-500 dark:bg-gray-700 dark:text-white dark:border-gray-600 mb-4"
              value={smilesInput}
              readOnly
            />
          )}

          <div className="flex justify-center md:justify-start">
            <button
              onClick={handlePredict}
              disabled={loading}
              className="bg-green-600 text-white px-6 py-3 rounded-lg shadow-lg hover:bg-green-700 transition-all duration-300 dark:bg-green-500 dark:hover:bg-green-600"
            >
              {loading ? "Predicting..." : "Predict"}
            </button>
          </div>
        </motion.div>
      </div>
    </motion.div>
  );
};

export default Predict;
