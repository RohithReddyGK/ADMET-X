import './index.css';
import { useEffect, useState } from "react";
import { motion, AnimatePresence } from "framer-motion";
import Header from "./components/Header";
import MainContent from "./components/MainContent";
import Predict from "./components/Predict";
import Footer from "./components/Footer";
import introLogo from "./assets/ADMET-AI Into_Logo.png";

function App() {
  const [showIntro, setShowIntro] = useState(true);

  useEffect(() => {
    const timer = setTimeout(() => {
      setShowIntro(false);
    }, 3000); // 3 seconds

    return () => clearTimeout(timer);
  }, []);

  return (
    <div className="w-full min-h-screen bg-white text-black dark:bg-gray-900 dark:text-white transition-all duration-500">
      {showIntro ? (
        <div className="flex flex-col justify-center items-center h-screen bg-gradient-to-r from-blue-500 to-cyan-500 transition-all duration-1000">
          <img
            src={introLogo}
            alt="Intro Logo"
            className="w-40 h-40 animate-bounce drop-shadow-2xl"
          />
          <div className="mt-6 text-white text-lg animate-pulse">
            Loading ADMET-AI.....
          </div>
        </div>
      ) : (
        <div className="bg-white transition-all duration-1000">
          <Header />
          <MainContent />
          <Predict />
          <div className="h-8" />
          <Footer />
        </div>
      )}
    </div>
  );
};

export default App;
