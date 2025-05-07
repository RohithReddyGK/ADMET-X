import React, {useState,useEffect} from "react";
import {FaSun, FaMoon} from "react-icons/fa";

const ThemeToggle=()=>{
    const [darkMode,setdarkMode]=useState(false);

    useEffect(()=>{
        const storedTheme=localStorage.getItem("theme");
        if (storedTheme==="dark") {
            setdarkMode(true);
            document.documentElement.classList.add("dark");
        }
        else {
            setdarkMode(false);
            document.documentElement.classList.remove("dark");
        }
    },[]);

    const toggleDarkMode=()=>{
        setdarkMode((prevMode)=>{
            const newMode=!prevMode;
            localStorage.setItem("theme",newMode ? "dark":"light");

            if (newMode) {
                document.documentElement.classList.add("dark");
            }
            else {
                document.documentElement.classList.remove("dark");
            }

            return newMode;
        });
    };

    return (
        <button
        onClick={toggleDarkMode}
        className="text-2xl p-2 focus:outline-none"
        >
            {darkMode ? (
                <FaSun className="text-yellow-500" />
            ):(
                <FaMoon className="text-blue-500" />
            )}
        </button>
    );
};

export default ThemeToggle;
