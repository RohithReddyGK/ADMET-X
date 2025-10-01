import { useEffect, useRef, useState } from "react";
import { motion } from "framer-motion";

const DrawMolecule = ({ onSmilesGenerated }) => {
  const jsmeRef = useRef(null);
  const [drawnSmiles, setDrawnSmiles] = useState([]);
  const [jsmeApplet, setJsmeApplet] = useState(null);

  useEffect(() => {
    window.jsmeOnLoad = function () {
      if (window.JSApplet && window.JSApplet.JSME) {
        if (jsmeRef.current && !jsmeRef.current.hasChildNodes()) {
          const applet = new window.JSApplet.JSME(
            jsmeRef.current.id,
            "100%", // full width for mobile
            "300px" // smaller height on mobile
          );
          setJsmeApplet(applet);
        }
      } else console.error("JSApplet.JSME is not available");
    };

    if (!document.getElementById("jsme_script")) {
      const script = document.createElement("script");
      script.id = "jsme_script";
      script.src = "/jsme/jsme.nocache.js";
      script.type = "text/javascript";
      script.async = true;
      document.body.appendChild(script);
    } else {
      if (typeof window.jsmeOnLoad === "function") window.jsmeOnLoad();
    }

    return () => {
      delete window.jsmeOnLoad;
    };
  }, []);

  const handleAddMolecule = () => {
    if (!jsmeApplet) return;
    const smi = jsmeApplet.smiles();
    if (smi && !drawnSmiles.includes(smi)) {
      const newList = [...drawnSmiles, smi];
      setDrawnSmiles(newList);
      onSmilesGenerated(newList.join("\n"));
    }
  };

  const handleClear = () => {
    setDrawnSmiles([]);
    onSmilesGenerated("");
    if (jsmeApplet) jsmeApplet.reset();
  };

  return (
    <motion.div
      className="p-4 animate-slide-up"
      initial={{ opacity: 0, y: 30 }}
      whileInView={{ opacity: 1, y: 0 }}
      viewport={{ once: true }}
      transition={{ duration: 0.6 }}
    >
      <h2 className="text-xl mb-2 font-semibold text-gray-900 dark:text-white">
        Draw Molecule
      </h2>

      <div
        id="jsme_container"
        ref={jsmeRef}
        className="w-full border rounded overflow-hidden mb-2"
      ></div>

      <div className="flex flex-col sm:flex-row gap-2 mt-2">
        <button
          onClick={handleAddMolecule}
          className="flex-1 px-4 py-2 bg-blue-600 text-white rounded hover:bg-blue-700 transition"
        >
          Add Molecule
        </button>
        <button
          onClick={handleClear}
          className="flex-1 px-4 py-2 bg-red-600 text-white rounded hover:bg-red-700 transition"
        >
          Clear All
        </button>
      </div>

      {drawnSmiles.length > 0 && (
        <div className="mt-2">
          <strong className="text-gray-900 dark:text-white">
            Drawn Molecules:
          </strong>
          <ul className="list-disc pl-6 text-gray-800 dark:text-gray-200 break-words">
            {drawnSmiles.map((smi, i) => (
              <li key={i}>{smi}</li>
            ))}
          </ul>
        </div>
      )}
    </motion.div>
  );
};

export default DrawMolecule;
