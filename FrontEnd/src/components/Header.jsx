import logo from "../assets/ADMET-AI Logo.png";
import ThemeToggle from "./ThemeToggle";
import { motion } from "framer-motion";

function Header() {
  return (
    <motion.header
      className="bg-blue-100 dark:bg-gray-800 py-5 shadow-md"
      initial={{ opacity: 0, y: -50 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 1 }}
    >
      <div className="container mx-auto flex items-center justify-between gap-4 px-4">
        <img
          src={logo}
          alt="ADMET-AI Logo"
          className="w-20 h-20 rounded-full"
        />
        <h1 className="whitespace-nowrap text-3xl md:text-5xl font-bold text-gray-800 dark:text-white text-center">
          AI-Driven Platform for In-Silico ADMET Prediction
        </h1>
        <ThemeToggle />
      </div>
    </motion.header>
  );
}

export default Header;
