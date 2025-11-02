import { FaReact } from "react-icons/fa";
import { SiTailwindcss, SiFlask } from "react-icons/si";
import { MdScience } from "react-icons/md";
import { TbDatabaseCog } from "react-icons/tb";
import { motion } from "framer-motion";
import React, { useState, useEffect } from "react";

const Footer = () => {
  const [showCite, setShowCite] = useState(false);
  const [accessedDate, setAccessedDate] = useState("");

  useEffect(() => {
    const date = new Date();
    const month = date.toLocaleDateString("default", { month: 'long' });
    const year = date.getFullYear();
    setAccessedDate(`${month} ${year}`);
  }, []);

  const copyToClipboard = (text) => {
    navigator.clipboard.writeText(text);
    alert("Copied to Clipboard!");
  };

  const generalCitation = `Rohith Reddy.G.K, Sheik Arshad Ibrahim, et al. (2025). ADMET-X: AI-Driven Platform for In-Silico ADMET Prediction. [Online]. Available: https://admet-x.vercel.app`;

  const ieeeCitation = `[1] Sheik Arshad Ibrahim, Rohith Reddy.G.K, et al., “ADMET-X: AI-Driven Platform for In-Silico ADMET Prediction,” 2025. [Online]. Available: https://admet-x.vercel.app`;

  const bibtexCitation = `@software{admetx2025,
  author = {Sheik Arshad Ibrahim, Rohith Reddy.G.K and others},
  title = {ADMET-X: AI-Driven Platform for In-Silico ADMET Prediction},
  year = {2025},
  url = {https://admet-x.vercel.app},
  note = {Accessed: ${accessedDate}}
}`;

  return (
    <>
      <motion.footer
        className="bg-blue-100 dark:bg-gray-800 text-gray-800 dark:text-white py-1 shadow-md"
        initial={{ opacity: 0, y: 50 }}
        whileInView={{ opacity: 1, y: 0 }}
        viewport={{ once: true }}
        transition={{ duration: 1 }}
      >
        <div className="max-w-6xl mx-auto px-6 flex flex-col sm:flex-row justify-between gap-4 sm:gap-16">

          {/* ABOUT */}
          <div className="flex-1 flex flex-col items-center sm:items-start sm:ml-[-2rem]">
            <h3 className="font-semibold text-lg mb-3">ABOUT</h3>
            <p className="text-lgs dark:text-gray-300 mb-2 text-center sm:text-left">
              <a
                href="https://admet-x.vercel.app"
                target=" _blank"
                rel="noopener noreferrer"
                className="text-blue-600 dark:text-blue-400 hover:underline font-semibold text:justify"
              >
                ADMET-X
              </a>{" "}
              is an AI-powered platform for predicting ADMET - Absorption,
              Distribution, Metabolism, Excretion, and Toxicity properties of drug
              molecules, aiding faster drug discovery.
            </p>
            <a
              href="https://github.com/RohithReddyGK/ADMET-X"
              target="_blank"
              rel="noopener noreferrer"
              className="text-blue-600 dark:text-blue-400 hover:underline font-semibold"
            >
              View Project on GitHub
            </a>
          </div>

          {/* DEV TOOLS */}
          <div className="flex-1 flex flex-col items-center sm:items-start sm:ml-6">
            <h3 className="font-semibold text-lg mb-3">DEV TOOLS:</h3>
            <div className="flex gap-5 text-5xl sm:text-5xl">
              <FaReact className="text-cyan-500 dark:text-cyan-400 hover:scale-110 transition-transform duration-200" title="React" />
              <SiTailwindcss className="text-sky-400 dark:text-sky-400 hover:scale-110 transition-transform duration-200" title="Tailwind CSS" />
              <SiFlask className="text-black dark:text-white hover:scale-110 transition-transform duration-200" title="Flask" />
              <MdScience className="text-green-600 dark:text-green-400 hover:scale-110 transition-transform duration-200" title="RDKit" />
              <TbDatabaseCog className="text-purple-600 dark:text-purple-400 hover:scale-110 transition-transform duration-200" title="TDC" />
            </div>
          </div>

          {/* CITE ADMET-X */}
          <div className="flex-1 flex flex-col items-center sm:items-start sm:ml-12">
            <h3 className="font-semibold text-lg mb-3 text-center sm:text-left">CURIOUS TO CITE
              <a
                href="https://admet-x.vercel.app"
                target="_blank"
                rel="noopener noreferrer"
                className="text-blue-600 dark:text-blue-400 hover:underline font-semibold mx-1"
              >
                ADMET-X
              </a>{" "}:
            </h3>
            <div className="w-full flex justify-center sm:justify-start">
              <button
                onClick={() => setShowCite(true)}
                className="bg-blue-600 text-white px-4 py-2 rounded-lg hover:bg-blue-700 transition duration-200"
              >
                Cite ADMET-X
              </button>
            </div>
          </div>
        </div>

        <div className="mt-4 text-center text-lg font-semibold dark:text-gray-400">
          © {new Date().getFullYear()}
          <a
            href="https://admet-x.vercel.app"
            target="_blank"
            rel="noopener noreferrer"
            className="text-blue-600 dark:text-blue-400 hover:underline font-semibold mx-1"
          >
            ADMET-X
          </a>{" "}
          — All rights reserved.
        </div>
      </motion.footer>

      {/* Modal Popup for Cite ADMET-X */}
      {showCite && (
        <div className="fixed inset-0 bg-black bg-opacity-50 flex justify-center items-center z-50">
          <div className="bg-white dark:bg-gray-800 text-gray-800 dark:text-gray-200 rounded-2xl shadow-2xl w-11/12 sm:w-2/3 md:w-1/2 p-6 relative">
            <button
              onClick={() => setShowCite(false)}
              className="absolute top-3 right-3 text-gray-500 hover:text-gray-800 dark:hover:text-gray-300 text-2xl"
            >
              ✕
            </button>
            <h2 className="text-2xl font-bold text-center mb-6">🧬 How to Cite ADMET-X?</h2>

            {/* General Citation */}
            <div className="mb-5">
              <h4 className="font-semibold mb-1">General Citation Format</h4>
              <textarea
                readOnly
                value={generalCitation}
                className="text-sm bg-gray-100 dark:bg-white dark:text-black p-2 rounded-md w-full resize-none mb-2"
                rows="2"
              />
              <button
                onClick={() => copyToClipboard(generalCitation)}
                className="bg-blue-500 text-white text-sm px-3 py-1 rounded hover:bg-blue-600 float-right"
              >
                Copy
              </button>
            </div>

            {/* IEEE Format */}
            <div className="mb-5 clear-both">
              <h4 className="font-semibold mb-1">IEEE Format</h4>
              <textarea
                readOnly
                value={ieeeCitation}
                className="text-sm bg-gray-100 dark:bg-white dark:text-black p-2 rounded-md w-full resize-none mb-2"
                rows="2"
              />
              <button
                onClick={() => copyToClipboard(ieeeCitation)}
                className="bg-blue-500 text-white text-sm px-3 py-1 rounded hover:bg-blue-600 float-right"
              >
                Copy
              </button>
            </div>

            {/* BibTeX Format */}
            <div className="clear-both">
              <h4 className="font-semibold mb-1">BibTeX Format</h4>
              <textarea
                readOnly
                value={bibtexCitation}
                className="text-sm bg-gray-100 dark:bg-white dark:text-black p-2 rounded-md w-full resize-none mb-2"
                rows="4"
              />
              <button
                onClick={() => copyToClipboard(bibtexCitation)}
                className="bg-blue-500 text-white text-sm px-3 py-1 rounded hover:bg-blue-600 float-right"
              >
                Copy
              </button>
            </div>
          </div>
        </div>
      )}
    </>
  );
};

export default Footer;
