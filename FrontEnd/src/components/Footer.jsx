import { FaReact, FaGithub } from "react-icons/fa";
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
    const month = date.toLocaleDateString("default", { month: "long" });
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
        className="bg-blue-100 dark:bg-gray-800 text-gray-800 dark:text-white py-4 shadow-md border-t border-gray-300 dark:border-gray-700"
        initial={{ opacity: 0, y: 30 }}
        whileInView={{ opacity: 1, y: 0 }}
        viewport={{ once: true }}
        transition={{ duration: 0.6 }}
      >
        <div className="max-w-6xl mx-auto px-4 flex flex-col sm:flex-row justify-between items-start sm:items-center gap-6 sm:gap-4">

          {/* ABOUT SECTION */}
          <div className="flex-1 text-center sm:text-left space-y-2 pb-4 border-b border-gray-300 dark:border-gray-700 sm:border-none">
            <h3 className="font-semibold text-lg tracking-wide">ABOUT</h3>
            <p className="text-sm dark:text-gray-300 leading-relaxed">
              <a
                href="https://admet-x.vercel.app"
                target="_blank"
                rel="noopener noreferrer"
                className="text-blue-600 dark:text-blue-400 font-semibold hover:underline"
              >
                ADMET-X
              </a>{" "}
              is an AI-powered platform for predicting ADMET properties — Absorption,
              Distribution, Metabolism, Excretion and Toxicity — of drug molecules,
              supporting accelerated drug discovery.
            </p>
            <div className="flex flex-col sm:flex-row sm:items-center sm:gap-3 justify-center sm:justify-start text-sm">
              <a
                href="https://github.com/RohithReddyGK/ADMET-X"
                target="_blank"
                rel="noopener noreferrer"
                className="flex items-center gap-1 text-blue-600 dark:text-blue-400 hover:underline font-semibold"
              >
                <FaGithub className="text-black dark:text-white text-xl" />
                <span>Lead Developer's GitHub</span>
              </a>
              <a
                href="https://github.com/arshadibrahim882/ADMET-X"
                target="_blank"
                rel="noopener noreferrer"
                className="flex items-center gap-1 text-blue-600 dark:text-blue-400 hover:underline font-semibold"
              >
                <FaGithub className="text-black dark:text-white text-xl" />
                <span>Contributor's GitHub</span>
              </a>
            </div>
          </div>

          {/* DEV TOOLS */}
          <div className="flex-1 text-center sm:text-left pb-4 border-b border-gray-300 dark:border-gray-700 sm:border-none">
            <h3 className="font-semibold text-lg mb-2 tracking-wide">DEV TOOLS</h3>
            <div className="flex justify-center sm:justify-start gap-4 text-4xl">
              <FaReact className="text-cyan-500 dark:text-cyan-400 hover:scale-110 transition-transform" title="React" />
              <SiTailwindcss className="text-sky-400 dark:text-sky-400 hover:scale-110 transition-transform" title="Tailwind CSS" />
              <SiFlask className="text-black dark:text-white hover:scale-110 transition-transform" title="Flask" />
              <MdScience className="text-green-600 dark:text-green-400 hover:scale-110 transition-transform" title="RDKit" />
              <TbDatabaseCog className="text-purple-600 dark:text-purple-400 hover:scale-110 transition-transform" title="TDC" />
            </div>
          </div>

          {/* CITE SECTION */}
          <div className="flex-1 text-center sm:text-left pb-4 border-b border-gray-300 dark:border-gray-700 sm:border-none">
            <h3 className="font-semibold text-lg mb-2 tracking-wide">
              CURIOUS TO CITE ADMET-X?
            </h3>
            <button
              onClick={() => setShowCite(true)}
              className="bg-blue-600 text-white px-4 py-2 text-sm rounded-lg hover:bg-blue-700 transition"
            >
              Cite ADMET-X
            </button>
          </div>
        </div>

        <div className="mt-4 text-center text-sm font-bold dark:text-gray-400">
          © {new Date().getFullYear()} ADMET-X Team — All rights reserved.
        </div>
      </motion.footer>

      {/* Citation Modal */}
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
