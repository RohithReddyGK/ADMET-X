import { useRef, useState } from "react";
import PropertyInfo from "./PropertyInfo";
import { propertyInfo } from "../data/PropertyInfo";

const ToolTipCell = ({ name }) => {
  const displayName = name?.replace(/_/g, " ") || "Unknown";
  const cellRef = useRef(null);
  const [showTooltip, setShowTooltip] = useState(false);

  return (
    <div
      ref={cellRef}
      className="relative inline-block cursor-pointer text-back-700 font-medium underline dark:text-black-300 hover:text-blue-900 dark:hover:text-blue-400"
      onMouseEnter={() => setShowTooltip(true)}
      onMouseLeave={() => setShowTooltip(false)}
    >
      {displayName}
      <PropertyInfo
        name={name}
        info={propertyInfo[name]}
        visible={showTooltip}
        anchorRef={cellRef}
      />
    </div>
  );
};

export default ToolTipCell;