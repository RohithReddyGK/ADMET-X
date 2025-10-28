import { FaReact, FaLinkedin, FaGithub } from "react-icons/fa";
import { SiTailwindcss, SiFlask } from "react-icons/si";
import { MdScience } from "react-icons/md";
import { TbDatabaseCog } from "react-icons/tb";
import { motion } from "framer-motion";

const Footer = () => {
  return (
    <motion.footer
      className="bg-blue-100 dark:bg-gray-800 text-gray-800 dark:text-white py-10 shadow-md"
      initial={{ opacity: 0, y: 50 }}
      whileInView={{ opacity: 1, y: 0 }}
      viewport={{ once: true }}
      transition={{ duration: 1 }}
    >
      <div className="max-w-6xl mx-auto px-6 flex flex-col sm:flex-row justify-between gap-10">

        {/* Column 1: ABOUT */}
        <div className="flex-1 flex flex-col items-center sm:items-start">
          <h3 className="font-semibold text-lg mb-3">ABOUT</h3>
          <p className="text-lgs dark:text-gray-300 mb-2 text-center sm:text-left">
            <a
              href="https://admet-x.vercel.app"
              target=" _blank"
              rel="noopener noreferrer"
              className="text-blue-600 dark:text-blue-400 hover:underline font-semibold"
            >
              ADMET-X
            </a>{" "}
            is an AI-powered platform for predicting ADMET - Absorption, Distribution, Metabolism, Excretion, and Toxicity properties of drug molecules, aiding faster drug discovery.
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

        {/* Column 2: DEV TOOLS */}
        <div className="flex-1 flex flex-col items-center sm:items-start">
          <h3 className="font-semibold text-lg mb-3">DEV TOOLS:</h3>
          <div className="flex gap-5 text-5xl sm:text-5xl">
            <FaReact className="text-cyan-500 dark:text-cyan-400 hover:scale-110 transition-transform duration-200" title="React" />
            <SiTailwindcss className="text-sky-400 dark:text-sky-400 hover:scale-110 transition-transform duration-200" title="Tailwind CSS" />
            <SiFlask className="text-black dark:text-white hover:scale-110 transition-transform duration-200" title="Flask" />
            <MdScience className="text-green-600 dark:text-green-400 hover:scale-110 transition-transform duration-200" title="RDKit" />
            <TbDatabaseCog className="text-purple-600 dark:text-purple-400 hover:scale-110 transition-transform duration-200" title="TDC" />
          </div>
        </div>

      </div>

      {/* Bottom Center: Rights */}
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
    </motion.footer >
  );
};

export default Footer;
