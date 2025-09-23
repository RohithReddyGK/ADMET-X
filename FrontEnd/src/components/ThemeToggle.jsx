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
    setDarkMode(prev => {
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
      className="text-2xl p-2 focus:outline-none"
      whileHover={{ scale: 1.2 }}
      whileTap={{ scale: 0.9 }}
    >
      {darkMode ? <FaSun className="text-yellow-500" /> : <FaMoon className="text-blue-500" />}
    </motion.button>
  );
};

export default ThemeToggle;
