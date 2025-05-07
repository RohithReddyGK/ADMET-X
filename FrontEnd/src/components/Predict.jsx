import { useState } from "react";
import DrawMolecule from "./DrawMolecule";

const Predict = () => {
    const [selectedOption, setSelectedOption] = useState("text");
    const [smilesInput, setSmilesInput] = useState("");
    const [file, setFile] = useState(null);

    const handleFileChange = (e) => {
        const uploadedFile = e.target.files[0];
        if (uploadedFile) {
            const fileType = uploadedFile.type;

            if (fileType === "text/plain" || fileType === "text/csv") {
                setFile(uploadedFile);

                const reader = new FileReader();
                reader.onload = (event) => {
                    const text = event.target.result;

                    // If it's a CSV, parse it line by line (simplified handling for CSV).
                    if (fileType === "text/csv") {
                        const parsedSMILES = text.split("\n").map(line => line.trim());
                        setSmilesInput(parsedSMILES.join("\n"));
                    } else {
                        setSmilesInput(text); // For plain text, directly set input.
                    }
                };
                reader.readAsText(uploadedFile);
            } else {
                alert("Please upload a valid .txt or .csv file.");
            }
        }
    };

    const handleExample = () => {
        const example = `CC(=O)NC1=CC=C(C=C1)O`;
        setSmilesInput(example);
    };

    const handlePredict = () => {
        if (!smilesInput.trim()) {
            alert("Please provide at least one SMILE input.");
            return;
        }
        // Send "smilesInput" to backend for prediction.
        console.log("Predicting for SMILES:\n", smilesInput);
    };

    return (
        <div className="p-6 max-w-3xl mx-auto bg-white dark:bg-gray-800 rounded-lg shadow-2xl">
            <div className="container mx-auto">
                <h1 className="text-3xl font-bold mb-4 text-center">Prediction of ADMET Properties</h1>

                {/* Tabs */}
                <div className="flex justify-center space-x-4 mb-6">
                    {["text", "file", "draw", "example"].map((option) => (
                        <button
                            key={option}
                            className={`px-6 py-2 rounded font-medium transition ${selectedOption === option
                                    ? "bg-blue-600 text-white dark:bg-blue-800 dark:text-white"
                                    : "bg-blue-100 text-blue-800 dark:bg-gray-700 dark:text-white"
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

                {/* Content */}
                <div className="mb-6">
                    {selectedOption === "text" && (
                        <textarea
                            className="w-full h-40 p-2 border rounded-lg focus:ring-2 focus:ring-blue-500 dark:bg-gray-700 dark:text-white dark:border-gray-600"
                            placeholder="Paste SMILES here (one per line)..."
                            value={smilesInput}
                            onChange={(e) => setSmilesInput(e.target.value)}
                        />
                    )}

                    {selectedOption === "file" && (
                        <div className="space-y-3">
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
                        <div className="w-full p-2 border rounded-lg">
                            <DrawMolecule onSmilesGenerated={(smiles) => setSmilesInput(smiles)} />
                        </div>
                    )}

                    {selectedOption === "example" && (
                        <textarea
                            className="w-full h-40 p-2 border rounded-lg focus:ring-2 focus:ring-blue-500 dark:bg-gray-700 dark:text-white dark:border-gray-600"
                            value={smilesInput}
                            readOnly
                        />
                    )}
                </div>

                <div className="flex justify-center">
                    <button
                        onClick={handlePredict}
                        className="bg-green-600 text-white px-6 py-3 rounded-lg shadow-lg hover:bg-green-700 transition-all duration-300 dark:bg-green-500 dark:hover:bg-green-600"
                    >
                        Predict
                    </button>
                </div>
            </div>
        </div>
    );
};

export default Predict;
