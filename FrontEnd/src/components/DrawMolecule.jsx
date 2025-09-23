import { useEffect, useRef } from 'react';
import { motion } from "framer-motion";

const DrawMolecule = ({ onSmilesGenerated }) => {
  const jsmeRef = useRef(null);

  useEffect(() => {
    window.jsmeOnLoad = function () {
      if (window.JSApplet && window.JSApplet.JSME) {
        if (jsmeRef.current && !jsmeRef.current.hasChildNodes()) {
          new window.JSApplet.JSME(jsmeRef.current.id, '500px', '400px');
        }
      } else console.error('JSApplet.JSME is not available');
    };

    if (!document.getElementById('jsme_script')) {
      const script = document.createElement('script');
      script.id = 'jsme_script';
      script.src = '/jsme/jsme.nocache.js';
      script.type = 'text/javascript';
      script.async = true;
      document.body.appendChild(script);
    } else {
      if (typeof window.jsmeOnLoad === 'function') window.jsmeOnLoad();
    }

    return () => { delete window.jsmeOnLoad; };
  }, []);

  return (
    <motion.div
      className="p-4 animate-slide-up"
      initial={{ opacity: 0, y: 30 }}
      whileInView={{ opacity: 1, y: 0 }}
      viewport={{ once: true }}
      transition={{ duration: 0.6 }}
    >
      <h2 className="text-xl mb-2 font-semibold">Draw Molecule</h2>
      <div id="jsme_container" ref={jsmeRef}></div>
    </motion.div>
  );
};

export default DrawMolecule;
