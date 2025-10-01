// Conclusion.jsx
import React from "react";

// Tag component for Good/Moderate/Poor
const StatusTag = ({ status }) => {
  let bgColor = "";
  let textColor = "text-white";
  let label = status;

  switch (status) {
    case "Good":
      bgColor = "bg-green-500";
      label = "Good";
      break;
    case "Moderate":
      bgColor = "bg-yellow-400";
      textColor = "text-black";
      label = "Moderate";
      break;
    case "Poor":
      bgColor = "bg-red-500";
      label = "Poor";
      break;
    default:
      bgColor = "bg-gray-400";
      label = "N/A";
  }

  return (
    <span
      className={`px-3 py-1 rounded-full font-semibold ${bgColor} ${textColor}`}
    >
      {label}
    </span>
  );
};

const Conclusion = ({ admet }) => {
  if (!admet || Object.keys(admet).length === 0) {
    return (
      <div className="mt-6 p-4 rounded bg-blue-50 dark:bg-blue-900 text-gray-900 dark:text-gray-100">
        <h4 className="font-bold text-lg mb-3">Conclusion</h4>
        <p>ℹ️ No ADMET data available.</p>
      </div>
    );
  }

  // Compute overall status for each ADMET category
  const categoryStatus = Object.keys(admet).map((cat) => {
    const catRows = admet[cat];
    const statusCount = { green: 0, yellow: 0, red: 0 };

    catRows.forEach((row) => {
      if (row.status === "green") statusCount.green++;
      else if (row.status === "yellow") statusCount.yellow++;
      else if (row.status === "red") statusCount.red++;
    });

    let overallStatus = "Good"; // default
    if (statusCount.red > statusCount.green + statusCount.yellow)
      overallStatus = "Poor";
    else if (
      statusCount.yellow > statusCount.green &&
      statusCount.yellow >= statusCount.red
    )
      overallStatus = "Moderate";

    return { category: cat, overallStatus };
  });

  // --- Market Release Suitability Logic ---
  let suitability = { message: "", bgColor: "", textColor: "" };

  const toxicity = categoryStatus.find((c) => c.category === "Toxicity");
  const absorption = categoryStatus.find((c) => c.category === "Absorption");
  const metabolism = categoryStatus.find((c) => c.category === "Metabolism");

  const anyModerate = categoryStatus.some((c) => c.overallStatus === "Moderate");
  const allGood = categoryStatus.every((c) => c.overallStatus === "Good");

  if (
    (toxicity && toxicity.overallStatus === "Poor") ||
    (absorption && absorption.overallStatus === "Poor") ||
    (metabolism && metabolism.overallStatus === "Poor")
  ) {
    suitability = {
      message:
        "This drug shows critical ADMET issues (e.g., Poor Toxicity/Absorption/Metabolism) and is unlikely suitable for market release.",
      bgColor: "bg-red-500",
      textColor: "text-white",
    };
  } else if (anyModerate) {
    suitability = {
      message:
        "This drug shows a moderate ADMET profile. Further optimization may be required before market evaluation.",
      bgColor: "bg-yellow-400",
      textColor: "text-black",
    };
  } else if (allGood) {
    suitability = {
      message:
        "This drug demonstrates favorable ADMET properties and is likely suitable for market evaluation.",
      bgColor: "bg-green-500",
      textColor: "text-white",
    };
  }

  return (
    <div className="mt-6 p-4 rounded-lg shadow bg-blue-50 dark:bg-gray-800 text-gray-900 dark:text-gray-100">
      <h4 className="font-bold text-lg mb-3">Overall Drug Conclusion</h4>

      {/* ADMET Property Table */}
      <div className="overflow-x-auto">
        <table className="min-w-full border text-sm bg-white dark:bg-gray-900 rounded">
          <thead>
            <tr className="bg-gray-200 dark:bg-gray-700 text-center">
              <th className="px-2 py-1 border">ADMET Category</th>
              <th className="px-2 py-1 border">Overall Status</th>
            </tr>
          </thead>
          <tbody>
            {categoryStatus.map((c, idx) => (
              <tr key={idx}>
                <td className="px-2 py-1 border">{c.category}</td>
                <td className="px-2 py-1 border text-center">
                  <StatusTag status={c.overallStatus} />
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      {/* Market Suitability Message */}
      <div
        className={`mt-4 p-4 rounded-lg text-center font-semibold text-lg ${suitability.bgColor} ${suitability.textColor}`}
      >
        {suitability.message}
      </div>
    </div>
  );
};

export default Conclusion;
