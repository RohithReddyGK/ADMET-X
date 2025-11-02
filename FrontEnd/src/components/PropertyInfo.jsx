import React, { useRef, useEffect, useState } from "react";

const PropertyInfo = ({ info, visible, anchorRef }) => {
  const tooltipRef = useRef();
  const [coords, setCoords] = useState({ left: 0, top: 0, opacity: 0 });

  useEffect(() => {
    if (!visible) {
      setCoords((p) => ({ ...p, opacity: 0 }));
      return;
    }
    const trigger = anchorRef.current;
    const tooltip = tooltipRef.current;
    if (trigger && tooltip) {
      const triggerRect = trigger.getBoundingClientRect();
      const tooltipRect = tooltip.getBoundingClientRect();

      // Center tooltip horizontally above trigger
      let left = triggerRect.left + (triggerRect.width - tooltipRect.width) / 2;
      let top = triggerRect.top - tooltipRect.height - 8;

      // Prevent horizontal overflow
      if (left < 8) left = 8;
      if (left + tooltipRect.width > window.innerWidth - 8)
        left = window.innerWidth - tooltipRect.width - 8;

      // If not enough space above, place below
      if (top < 8) top = triggerRect.bottom + 8;
      setCoords({ left, top, opacity: 1 });
    }
  }, [visible, anchorRef, info]);

  if (!info || !visible) return null;

  return (
    <div
      ref={tooltipRef}
      style={{
        position: "fixed",
        left: coords.left,
        top: coords.top,
        opacity: coords.opacity,
        pointerEvents: "auto",
        zIndex: 9999,
        width: "320px",
        transition: "opacity 0.2s ease",
      }}
      className="bg-blue-50 text-gray-800 text-sm p-3 rounded-lg shadow-lg dark:bg-gray-900 dark:text-white"
    >
      <p className="mb-2">{info.description}</p>
      {info.link && (
        <a
          href={info.link}
          target="_blank"
          rel="noopener noreferrer"
          className="text-blue-600 dark:text-blue-400 font-medium hover:text-blue-800 dark:hover:text-blue-800"
        >
          🔗 Click Here
        </a>
      )}
      <div className="absolute left-1/2 top-full -translate-x-1/2 border-8 border-transparent dark:border-t-gray-900 border-t-white" />
    </div>
  );
};

export default PropertyInfo;