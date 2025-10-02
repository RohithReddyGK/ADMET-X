// Conclusion.jsx
import React from "react";
import {
  Radar, RadarChart, PolarGrid, PolarAngleAxis, PolarRadiusAxis, ResponsiveContainer
} from "recharts";

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

// Utility to shorten long SMILES for display
const shortenSmiles = (smi, maxLen = 20) => {
  if (!smi) return "";
  if (smi.length <= maxLen) return smi;
  return smi.slice(0, maxLen / 2) + "..." + smi.slice(-maxLen / 2);
};

const Conclusion = ({ admet, smiles, drugName }) => {
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

  // Convert status to numeric value for radar chart
  const radarData = categoryStatus.map((c) => {
    let value = 3;
    if (c.overallStatus === "Good") value = 3;
    else if (c.overallStatus === "Moderate") value = 2;
    else if (c.overallStatus === "Poor") value = 1;

    return { category: c.category, value, status: c.overallStatus };
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

  // Make a stable-ish unique gradient id per mounted component
  const gradientId = React.useMemo(
    () => "admetGradient-" + Math.random().toString(36).slice(2, 9),
    []
  );

  // drug label: prefer explicit prop, else admet.smiles, else empty
  const drugLabel = (smiles && String(smiles)) || admet?.smiles || "";

  return (
    <div className="mt-6 p-4 rounded-lg shadow bg-blue-50 dark:bg-gray-800 text-gray-900 dark:text-gray-100">
      <h4 className="font-bold text-2xl mb-3">Overall Drug Conclusion</h4>

      {/* ADMET Property Table with Blue Theme */}
      <div className="overflow-x-auto">
        <table className="min-w-full border text-sm text-gray-900 dark:text-gray-100 bg-blue-50 dark:bg-gray-900 rounded-lg overflow-hidden">
          <thead>
            <tr className="bg-blue-200 dark:bg-gray-800 text-gray-900 dark:text-gray-100 font-semibold">
              <th className="px-2 py-1 border border-blue-300 dark:border-gray-700">
                ADMET Category
              </th>
              <th className="px-2 py-1 border border-blue-300 dark:border-blue-700">
                Overall Status
              </th>
            </tr>
          </thead>
          <tbody>
            {categoryStatus.map((c, idx) => (
              <tr
                key={idx}
                className="odd:bg-blue-50 even:bg-blue-100 dark:odd:bg-gray-900 dark:even:bg-gray-800 hover:bg-gray-200 dark:hover:bg-gray-700"
              >
                <td className="px-2 py-1 border border-blue-300 dark:border-gray-700">
                  {c.category}
                </td>
                <td className="px-2 py-1 border border-blue-300 dark:border-gray-700 text-center">
                  <StatusTag status={c.overallStatus} />
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      {/* Spider (Radar) Plot */}
      <div className="mt-6 w-full h-72 p-4 rounded-lg bg-white shadow-lg relative">
        {/* ResponsiveContainer will size the chart; the parent is 'relative' so we can overlay text inside it */}
        <ResponsiveContainer>
          <RadarChart cx="50%" cy="50%" outerRadius="80%" data={radarData}>
            {/* Subtle dashed grid */}
            <PolarGrid stroke="#CBD5E1" strokeDasharray="3 3" />

            {/* Category axis labels */}
            <PolarAngleAxis
              dataKey="category"
              stroke="#1F2937"
              tick={{ fontSize: 12, fontWeight: 500 }}
            />
            <PolarRadiusAxis domain={[0, 3]} tick={false} />

            {/* Gradient fill for radar polygon to reflect Good/Moderate/Poor */}
            <defs>
              <linearGradient id={gradientId} x1="0" y1="0" x2="1" y2="1">
                <stop
                  offset="0%"
                  stopColor={
                    radarData[0]?.status === "Good"
                      ? "#22C55E"
                      : radarData[0]?.status === "Moderate"
                      ? "#FACC15"
                      : "#EF4444"
                  }
                  stopOpacity={0.4}
                />
                <stop
                  offset="50%"
                  stopColor={
                    radarData[Math.floor(radarData.length / 2)]?.status === "Good"
                      ? "#22C55E"
                      : radarData[Math.floor(radarData.length / 2)]?.status === "Moderate"
                      ? "#FACC15"
                      : "#EF4444"
                  }
                  stopOpacity={0.4}
                />
                <stop
                  offset="100%"
                  stopColor={
                    radarData[radarData.length - 1]?.status === "Good"
                      ? "#22C55E"
                      : radarData[radarData.length - 1]?.status === "Moderate"
                      ? "#FACC15"
                      : "#EF4444"
                  }
                  stopOpacity={0.4}
                />
              </linearGradient>
            </defs>

            {/* Radar polygon with gradient fill and soft shadow */}
            <Radar
              name="ADMET Status"
              dataKey="value"
              stroke="#2563EB"
              fill={`url(#${gradientId})`}
              fillOpacity={0.6}
              style={{ filter: "drop-shadow(0px 2px 4px rgba(0,0,0,0.2))" }}
            />
             {/* 👇 Center labels (Drug name + shortened SMILES, both with tooltips) */}
            {drugName && (
              <text
                x="50%"
                y="47%"  // shift slightly up
                textAnchor="middle"
                dominantBaseline="middle"
                style={{ fontSize: "13px", fontWeight: "700", fill: "#111827" }}
              >
                {drugName}
                <title>{drugName}</title>
              </text>
            )}
            <text
              x="50%"
              y={drugName ? "55%" : "50%"} // shift down if name exists
              textAnchor="middle"
              dominantBaseline="middle"
              style={{ fontSize: "11px", fontWeight: "600", fill: "#374151" }}
            >
              {shortenSmiles(smiles)}
              <title>{smiles}</title>
            </text>
          </RadarChart>
        </ResponsiveContainer>
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
