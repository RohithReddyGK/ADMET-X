import React, { useState, useEffect } from "react";
import { FaSun, FaMoon } from "react-icons/fa";
import { motion } from "framer-motion";

const ThemeToggle = () => {
  const [darkMode, setDarkMode] = useState(false);

  useEffect(() => {
    const storedTheme = localStorage.getItem("theme");
    if (storedTheme === "dark") {
      setDarkMode(true);
      document.documentElement.classList.add("dark");
    } else {
      setDarkMode(false);
      document.documentElement.classList.remove("dark");
    }
  }, []);

  const toggleDarkMode = () => {
    setDarkMode((prev) => {
      const newMode = !prev;
      localStorage.setItem("theme", newMode ? "dark" : "light");
      if (newMode) document.documentElement.classList.add("dark");
      else document.documentElement.classList.remove("dark");
      return newMode;
    });
  };

  return (
    <motion.button
      onClick={toggleDarkMode}
      className={`p-3 rounded-full shadow-lg focus:outline-none transition-all duration-300
        ${darkMode 
          ? "bg-gradient-to-r from-purple-700 to-pink-500 hover:from-purple-600 hover:to-pink-400"
          : "bg-gradient-to-r from-yellow-400 to-orange-300 hover:from-yellow-300 hover:to-orange-200"
        }`}
      whileTap={{ scale: 0.9 }}
    >
      <motion.div
        key={darkMode ? "sun" : "moon"} // re-mounts icon for toggle animation
        initial={{ rotate: -180, opacity: 0 }}
        animate={{ rotate: 0, opacity: 1 }}
        transition={{ type: "spring", stiffness: 300, damping: 20 }}
        whileHover={{
          rotate: [0, 360], // rotate full circle
          transition: { repeat: Infinity, duration: 2, ease: "linear" },
        }}
      >
        {darkMode ? (
          <FaSun className="text-white text-xl sm:text-2xl" />
        ) : (
          <FaMoon className="text-white text-xl sm:text-2xl" />
        )}
      </motion.div>
    </motion.button>
  );
};

export default ThemeToggle;
