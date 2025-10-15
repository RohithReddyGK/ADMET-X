import logo from "../assets/ADMET-X.png";
import ThemeToggle from "./ThemeToggle";
import { motion } from "framer-motion";
import { Helmet } from "react-helmet-async";

function Header() {
  return (
    <>
      {/* SEO Metadata for Header */}
      <Helmet>
        <title>ADMET-X | AI-Driven Platform for In-Silico ADMET Prediction</title>
        <meta
          name="description"
          content="AI-powered platform for predicting Absorption, Distribution, Metabolism, Excretion, and Toxicity (ADMET) properties of drug molecules using advanced machine learning and cheminformatics."
        />
        <meta name="author" content="Sheik Arshad Ibrahim, Rohith Reddy G K, Sayed Jahangir Ali, Thirumurugan" />
        <meta name="robots" content="index, follow" />
        <link rel="canonical" href="https://admet-x.vercel.app/" />
        <link rel="icon" href="/ADMET-X.png" type="image/png" />

      </Helmet>

      {/* Header UI */}
      <motion.header
        className="bg-blue-100 dark:bg-gray-800 py-5 shadow-md"
        initial={{ opacity: 0, y: -50 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 1 }}
      >
        <div className="container mx-auto flex flex-col md:flex-row items-center justify-between gap-4 px-4">
          <img
            src={logo}
            alt="ADMET-X Logo"
            className="w-20 h-20 rounded-full mb-2 md:mb-0"
          />
          <h1 className="text-2xl md:text-4xl text-center md:text-left font-bold text-gray-800 dark:text-white">
            AI-Driven Platform for In-Silico ADMET Prediction
          </h1>
          <ThemeToggle />
        </div>
      </motion.header>
    </>
  );
}

export default Header;
