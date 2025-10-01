import './index.css';
import { useEffect, useState } from "react";
import { motion, AnimatePresence } from "framer-motion";
import Lottie from "lottie-react";
import Header from "./components/Header";
import MainContent from "./components/MainContent";
import Predict from "./components/Predict";
import Footer from "./components/Footer";
import introAnimation from "./assets/Doctor.json";

function App() {
  const [showIntro, setShowIntro] = useState(true);
  const [result, setResult] = useState(null);

  useEffect(() => {
    const timer = setTimeout(() => setShowIntro(false), 3000);
    return () => clearTimeout(timer);
  }, []);

  return (
    <div className="w-full min-h-screen bg-white text-black dark:bg-gray-900 dark:text-white transition-all duration-500">
      <AnimatePresence>
        {showIntro ? (
          <motion.div
            key="intro"
            initial={{ opacity: 0 }}
            animate={{ opacity: 2 }}
            exit={{ opacity: 0 }}
            className="flex justify-center items-center h-screen bg-white"
          >
            <Lottie
              animationData={introAnimation}
              loop={true}
              className="w-96 h-96 md:w-[800px] md:h-[800px]"
            />
          </motion.div>
        ) : (
          <motion.div
            key="main"
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            exit={{ opacity: 0 }}
            transition={{ duration: 1 }}
            className="bg-white transition-all duration-1000"
          >
            <Header />
            <MainContent />
            <Predict />
            <div className="h-8" />
            <Footer />
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}

export default App;
