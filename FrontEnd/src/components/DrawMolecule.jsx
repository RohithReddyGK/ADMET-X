import { useEffect, useRef } from 'react';

const DrawMolecule = () => {
  const jsmeRef = useRef(null);

  useEffect(() => {
    // Define global callback before script loads.
    window.jsmeOnLoad = function () {
      if (window.JSApplet && window.JSApplet.JSME) {
        // Only initialize if the container is empty.
        if (jsmeRef.current && !jsmeRef.current.hasChildNodes()) {
          new window.JSApplet.JSME(jsmeRef.current.id, '500px', '400px');
        }
      } else {
        console.error('JSApplet.JSME is not available');
      }
    };

    // Prevent duplicate script injection.
    if (!document.getElementById('jsme_script')) {
      const script = document.createElement('script');
      script.id = 'jsme_script';
      script.src = '/jsme/jsme.nocache.js';
      script.type = 'text/javascript';
      script.async = true;
      document.body.appendChild(script);
    } else {
      // If already loaded, call jsmeOnLoad manually.
      if (typeof window.jsmeOnLoad === 'function') {
        window.jsmeOnLoad();
      }
    }

    return () => {
      delete window.jsmeOnLoad;
    };
  }, []);

  return (
    <div className="p-4">
      <h2 className="text-xl mb-2 font-semibold">Draw Molecule</h2>
      <div id="jsme_container" ref={jsmeRef}></div>
    </div>
  );
};

export default DrawMolecule;
