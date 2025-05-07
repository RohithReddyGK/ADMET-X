import React from "react";
import { FaBrain, FaBolt } from "react-icons/fa";
import { GiMicroscope } from "react-icons/gi";
import { FiDatabase } from "react-icons/fi";

const MainContent = () => {
    return (
        // <div className="bg-white px-4 pt-6 pb-12 flex flex-col items-center text-center">
        <div className="bg-cover bg-center px-4 pb-12 flex flex-col items-center text-center" style={{backgroundImage: "url(/DNA.jpg)"}}>
            {/* Introduction Block */}
            <section className="max-w-3xl mb-7 animate-fade-in">
                <h2 className="text-4xl font-bold text-gray-800 dark:text-gray-600 mb-4">What is ADMET?</h2>
                <p className="text-lg text-gray-600 dark:text-gary-300 leading-relaxed">
                    ADMET stands for <strong>Absorption, Distribution, Metabolism, Excretion, and Toxicity</strong>.
                    This AI-Powered model utilizes Graph Neural Networks (GNNs) and Cheminformatics to predict these critical drug properties, accelerating drug discovery.
                </p>
            </section>

            {/* Key Advantages Grid */}
            <section className="grid grid-cols-1 sm:grid-cols-2 gap-8 max-w-4xl w-full animate-fade-in">
                <FeatureCard
                    icon={<GiMicroscope size={36} className="text-blue-600 dark:text-blue-400 mb-3" />}
                    title="In-Silico Screening"
                    description="Virtually screen candidate molecules to predict ADMET properties before investing in costly lab experiments. This accelerates early-stage drug research significantly."
                />
                <FeatureCard
                    icon={<FaBrain size={36} className="text-blue-600 dark:text-blue-400 mb-3" />}
                    title="AI-Driven Engine"
                    description="Leverages Graph Neural Networks (GNNs) built using PyTorch Geometric and DeepChem to capture molecular structure and predict outcomes with high precision."
                />
                <FeatureCard
                    icon={<FiDatabase size={36} className="text-blue-600 dark:text-blue-400 mb-3" />}
                    title="Trusted Data"
                    description={
                        <>
                            Trained using curated open-source datasets like:&nbsp;
                            <a
                                href="https://tdcommons.ai/"
                                target="_blank"
                                rel="noopener noreferrer"
                                className="text-blue-700 dark:text-blue-300 underline"
                            >
                                TDC
                            </a>{" "}
                            (Therapeutics Data Commons) and&nbsp;

                            <a
                                href="https://www.ebi.ac.uk/chembl/"
                                target="_blank"
                                rel="noopener noreferrer"
                                className="text-blue-700 dark:text-blue-300 underline"
                            >
                                ChEMBL
                            </a>, ensuring broad coverage and robustness.
                        </>
                    }
                />
                <FeatureCard
                    icon={<FaBolt size={36} className="text-blue-600 dark:text-blue-400 mb-3" />}
                    title="Fast Predictions"
                    description="Generate predictions within seconds, whether for a single molecule or in bulk, optimized by backend APIs, powered by Flask and RDKit processing."
                />
            </section>
        </div>
    );
};

const FeatureCard = ({ icon, title, description }) => (
    <div className="bg-gray-50 dark:bg-gray-800 p-6 rounded-xl shadow-md hover:shadow-lg transition transform hover:scale-110 text-gray-800 dark:text-gray-100">
        <div className="flex flex-col items-center justify-center">
            {icon}
            <h3 className="text-xl font-semibold mb-2 text-gray-900 dark:text-white">{title}</h3>
            <p className="text-sm text-justify text-gray-600 dark:text-gray-300">{description}</p>
        </div>
    </div>
);

export default MainContent;
