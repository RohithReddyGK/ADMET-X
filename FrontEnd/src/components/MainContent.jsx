import React from "react";
import { Helmet } from "react-helmet-async";
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
    <>
      {/* SEO & Structured Data for Main Content */}
      <Helmet>
        <meta
          name="keywords"
          content="ADMET Prediction, Drug Discovery, Machine Learning, Cheminformatics, Toxicity, Deep Learning, AI in Pharmacology, In-Silico Drug Design"
        />
        <meta property="og:title" content="ADMET-X | AI-Driven Platform for In-Silico Drug Discovery" />
        <meta
          property="og:description"
          content="Predict ADMET properties of drug molecules using AI and cheminformatics. Accelerate drug discovery with in-silico prediction tools."
        />
        <meta property="og:type" content="website" />
        <meta property="og:url" content="https://admet-x.vercel.app/" />
        <meta property="og:image" content="https://admet-x.vercel.app/ADMET-X.png" />
        <meta name="twitter:card" content="summary_large_image" />
        <meta name="twitter:title" content="ADMET-X | AI-Driven Platform for Drug Discovery" />
        <meta
          name="twitter:description"
          content="AI-based platform to predict Absorption, Distribution, Metabolism, Excretion, and Toxicity (ADMET) properties."
        />
        <meta name="twitter:image" content="https://admet-x.vercel.app/ADMET-X.png" />

        {/* Schema.org structured data for better Google results */}
        <script type="application/ld+json">
          {JSON.stringify({
            "@context": "https://schema.org",
            "@type": "SoftwareApplication",
            name: "ADMET-X",
            operatingSystem: "Web",
            applicationCategory: "AI Drug Discovery Tool",
            description:
              "AI-driven platform for predicting ADMET properties (Absorption, Distribution, Metabolism, Excretion, and Toxicity) of drug molecules.",
            url: "https://admet-x.vercel.app/",
            author: {
              "@type": "Person",
              name: "Rohith Reddy G K",
            },
            image: "https://admet-x.vercel.app/ADMET-X.png",
          })}
        </script>
      </Helmet>

      {/* Actual Content */}
      <div
        className="bg-cover bg-center px-4 pb-12 flex flex-col items-center text-center"
        style={{ backgroundImage: "url(/DNA.jpg)" }}
      >
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
            <strong>Absorption, Distribution, Metabolism, Excretion and Toxicity</strong>.
            This AI-Powered model utilizes Machine Learning Models and Cheminformatics
            to predict these critical drug properties, accelerating drug discovery.
          </p>
        </motion.section>

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
                <h3 className="text-xl font-semibold mb-2 text-gray-900 dark:text-white">
                  {card.title}
                </h3>
                <p className="text-sm text-center text-gray-600 dark:text-gray-300">
                  {card.description}
                </p>
              </div>
            </motion.div>
          ))}
        </div>
      </div>
    </>
  );
};

export default MainContent;
