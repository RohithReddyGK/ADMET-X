import React from "react";
import { FaBrain, FaBolt } from "react-icons/fa";
import { GiMicroscope } from "react-icons/gi";
import { FiDatabase } from "react-icons/fi";
import { motion } from "framer-motion";

const MainContent = () => {
  const cards = [
    {
      icon: <GiMicroscope size={36} className="text-blue-600 dark:text-blue-400 mb-3" />,
      title: "In-Silico Screening",
      description:
        "Virtually screen candidate molecules to predict ADMET properties before investing in costly lab experiments. This accelerates early-stage drug research significantly.",
    },
    {
      icon: <FaBrain size={36} className="text-blue-600 dark:text-blue-400 mb-3" />,
      title: "AI-Driven Engine",
      description:
        "Leverages Machine Learning algorithms like RandomForestClassifier, RandomForestRegressor, XGBoostClassifier and XGBoostRegressor for training models, for prediction of ADMET properties and to overcome the traditional methods.",
    },
    {
      icon: <FiDatabase size={36} className="text-blue-600 dark:text-blue-400 mb-3" />,
      title: "Trusted Data",
      description: (
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
            href="https://pubchem.ncbi.nlm.nih.gov/"
            target="_blank"
            rel="noopener noreferrer"
            className="text-blue-700 dark:text-blue-300 underline"
          >
            PubChem
          </a>
          , ensuring broad coverage and robustness.
        </>
      ),
    },
    {
      icon: <FaBolt size={36} className="text-blue-600 dark:text-blue-400 mb-3" />,
      title: "Fast Predictions",
      description:
        "Generate predictions within seconds, whether for a single molecule or in bulk, optimized by backend APIs, powered by Flask and RDKit processing.",
    },
  ];

  return (
    <div
      className="bg-cover bg-center px-4 pb-12 flex flex-col items-center text-center"
      style={{ backgroundImage: "url(/DNA.jpg)" }}
    >
      {/* Introduction Block */}
      <motion.section
        className="max-w-3xl mb-12"
        initial={{ opacity: 0, y: 30 }}
        whileInView={{ opacity: 1, y: 0 }}
        viewport={{ once: true }}
        transition={{ duration: 0.8 }}
      >
        <h2 className="text-4xl font-bold text-gray-800 dark:text-gray-800 mb-4">
          What is ADMET?
        </h2>
        <p className="text-lg text-gray-600 dark:text-gray-600 leading-relaxed">
          ADMET stands for{" "}
          <strong>Absorption, Distribution, Metabolism, Excretion, and Toxicity</strong>.
          This AI-Powered model utilizes Graph Neural Networks (GNNs) and Cheminformatics
          to predict these critical drug properties, accelerating drug discovery.
        </p>
      </motion.section>

      {/* Key Advantages Grid */}
      <div className="grid grid-cols-1 sm:grid-cols-2 gap-8 max-w-4xl w-full">
        {cards.map((card, index) => (
          <motion.div
            key={index}
            className="bg-gray-50 dark:bg-gray-800 p-6 rounded-xl shadow-md hover:shadow-lg transition transform hover:scale-110 text-gray-800 dark:text-gray-100"
            initial={{ opacity: 0, y: 30 }}
            whileInView={{ opacity: 1, y: 0 }}
            viewport={{ once: true }}
            transition={{ duration: 0.6, delay: index * 0.2 }}
          >
            <div className="flex flex-col items-center justify-center">
              {card.icon}
              <h3 className="text-xl font-semibold mb-2 text-gray-900 dark:text-white">{card.title}</h3>
              <p className="text-sm text-center text-gray-600 dark:text-gray-300">{card.description}</p>
            </div>
          </motion.div>
        ))}
      </div>
    </div>
  );
};

export default MainContent;
