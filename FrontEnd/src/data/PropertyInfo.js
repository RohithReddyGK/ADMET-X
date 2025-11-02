export const propertyInfo = {
  // Absorption 
  Bioavailability_Ma: {
    description: "Oral bioavailability is defined as “the rate and extent to which the active ingredient or active moiety is absorbed from a drug product and becomes available at the site of action”.",
    link: "https://tdcommons.ai/single_pred_tasks/adme/#bioavailability-ma-et-al"
  },
  Caco2_Wang: {
    description: "The human colon epithelial cancer cell, Caco-2 - is used as an in vitro model to simulate the human intestinal tissue. The experimental result on the rate of drug passing through the Caco-2 cells can approximate the rate at which the drug permeates through the human intestinal tissue.",
    link: "https://tdcommons.ai/single_pred_tasks/adme/#caco-2-cell-effective-permeability-wang-et-al"
  },
  HIA_Hou: {
    description: "When a drug is orally administered, it needs to be absorbed from the human gastrointestinal system into the bloodstream of the human body. This ability of absorption is called human intestinal absorption (HIA) and it is crucial for a drug to be delivered to the target.",
    link: "https://tdcommons.ai/single_pred_tasks/adme/#hia-human-intestinal-absorption-hou-et-al"
  },
  HydrationFreeEnergy_FreeSolv: {
    description: "The Free Solvation Database, FreeSolv(SAMPL), provides experimental and calculated hydration free energy of small molecules in water. The calculated values are derived from alchemical free energy calculations using molecular dynamics simulations.",
    link: "https://tdcommons.ai/single_pred_tasks/adme/#hydration-free-energy-freesolv"
  },
  Lipophilicity_AstraZeneca: {
    description: "Lipophilicity measures the ability of a drug to dissolve in a lipid (e.g. fats, oils) environment. High lipophilicity often leads to high rate of metabolism, poor solubility, high turn-over, and low absorption.",
    link: "https://tdcommons.ai/single_pred_tasks/adme/#lipophilicity-astrazeneca"
  },  
  PAMPA_NCATS: {
    description: "PAMPA (parallel artificial membrane permeability assay) is a commonly employed assay to evaluate drug permeability across the cellular membrane. Although PAMPA does not model active and efflux transporters, it still provides permeability values that are useful for absorption prediction because the majority of drugs are absorbed by passive diffusion through the membrane.",
    link: "https://tdcommons.ai/single_pred_tasks/adme/#pampa-permeability-ncats"    
  },
  Pgp_Broccatelli: {
    description: "P-glycoprotein (Pgp) is an ABC transporter protein involved in intestinal absorption, drug metabolism, and brain penetration, and its inhibition can seriously alter a drug's bioavailability and safety. In addition, inhibitors of Pgp can be used to overcome multidrug resistance.",
    link: "https://tdcommons.ai/single_pred_tasks/adme/#pgp-p-glycoprotein-inhibition-broccatelli-et-al"
  },
  Solubility_AqSolDB: {
    description: "Aqeuous solubility measures a drug's ability to dissolve in water. Poor water solubility could lead to slow drug absorptions, inadequate bioavailablity and even induce toxicity. More than 40% of new chemical entities are not soluble.",
    link: "https://tdcommons.ai/single_pred_tasks/adme/#solubility-aqsoldb"
  },

  // Distribution
  BBB_Martins: {
    description: "As a membrane separating circulating blood and brain extracellular fluid, the blood-brain barrier (BBB) is the protection layer that blocks most foreign drugs. Thus the ability of a drug to penetrate the barrier to deliver to the site of action forms a crucial challenge in development of drugs for central nervous system.",
    link: "https://tdcommons.ai/single_pred_tasks/adme/#bbb-blood-brain-barrier-martins-et-al"
  },
  PPBR_AZ: {
    description: "The human plasma protein binding rate (PPBR) is expressed as the percentage of a drug bound to plasma proteins in the blood. This rate strongly affect a drug's efficiency of delivery. The less bound a drug is, the more efficiently it can traverse and diffuse to the site of actions.",
    link:"https://tdcommons.ai/single_pred_tasks/adme/#ppbr-plasma-protein-binding-rate-astrazeneca"
  },
  VDss_Lombardo: {
    description: "The volume of distribution at steady state (VDss) measures the degree of a drug's concentration in body tissue compared to concentration in blood. Higher VD indicates a higher distribution in the tissue and usually indicates the drug with high lipid solubility, low plasma protein binidng rate.",
    link: "https://tdcommons.ai/single_pred_tasks/adme/#vdss-volumn-of-distribution-at-steady-state-lombardo-et-al"
  },

  // Excretion 
  Clearance_Hepatocyte_AZ: {
    description: "Drug clearance is defined as the volume of plasma cleared of a drug over a specified time period and it measures the rate at which the active drug is removed from the body.It contains clearance measures from two experiments types, hepatocyte and microsomes.",
    link: "https://tdcommons.ai/single_pred_tasks/adme/#clearance-astrazeneca"
  },
  Half_Life_Obach: {
    description: "Half life of a drug is the duration for the concentration of the drug in the body to be reduced by half. It measures the duration of actions of a drug.",
    link: "https://tdcommons.ai/single_pred_tasks/adme/#half-life-obach-et-al"
  },

  // Metabolism
  CYP1A2_Veith: {
    description: "The CYP P450 genes are involved in the formation and breakdown (metabolism) of various molecules and chemicals within cells. Specifically, CYP1A2 localizes to the endoplasmic reticulum and its expression is induced by some polycyclic aromatic hydrocarbons (PAHs), some of which are found in cigarette smoke.",
    link: "https://tdcommons.ai/single_pred_tasks/adme/#cyp-p450-1a2-inhibition-veith-et-al"
  },
  CYP2C19_Veith: {
    description: "The CYP P450 genes are essential in the breakdown (metabolism) of various molecules and chemicals within cells. A drug that can inhibit these enzymes would mean poor metabolism to this drug and other drugs, which could lead to drug-drug interactions and adverse effects.",
    link: "https://tdcommons.ai/single_pred_tasks/adme/#cyp-p450-2c19-inhibition-veith-et-al"
  },
  CYP2C9_Substrate_CarbonMangels: {
    description: "CYP P450 2C9 plays a major role in the oxidation of both xenobiotic and endogenous compounds. Substrates are drugs that are metabolized by the enzyme.",
    link: "https://tdcommons.ai/single_pred_tasks/adme/#cyp2c9-substrate-carbon-mangels-et-al"
  },
  CYP2C9_Veith: {
    description: "The CYP P450 genes are involved in the formation and breakdown (metabolism) of various molecules and chemicals within cells. Specifically, the CYP P450 2C9 plays a major role in the oxidation of both xenobiotic and endogenous compounds.",
    link: "https://tdcommons.ai/single_pred_tasks/adme/#cyp-p450-2c9-inhibition-veith-et-al"
  },
  CYP2D6_Substrate_CarbonMangels: {
    description: "CYP2D6 is primarily expressed in the liver. It is also highly expressed in areas of the central nervous system, including the substantia nigra.",
    link: "https://tdcommons.ai/single_pred_tasks/adme/#cyp2d6-substrate-carbon-mangels-et-al"
  },
  CYP2D6_Veith: {
    description: "The CYP P450 genes are involved in the formation and breakdown (metabolism) of various molecules and chemicals within cells. Specifically, CYP2D6 is primarily expressed in the liver. It is also highly expressed in areas of the central nervous system, including the substantia nigra.",
    link:"https://tdcommons.ai/single_pred_tasks/adme/#cyp-p450-2d6-inhibition-veith-et-al"
  },
  CYP3A4_Substrate_CarbonMangels: {
    description: "CYP3A4 is an important enzyme in the body, mainly found in the liver and in the intestine. It oxidizes small foreign organic molecules (xenobiotics), such as toxins or drugs, so that they can be removed from the body.",
    link: "https://tdcommons.ai/single_pred_tasks/adme/#cyp3a4-substrate-carbon-mangels-et-al"
  },
  CYP3A4_Veith: {
    description: "The CYP P450 genes are involved in the formation and breakdown (metabolism) of various molecules and chemicals within cells. Specifically, CYP3A4 is an important enzyme in the body, mainly found in the liver and in the intestine. It oxidizes small foreign organic molecules (xenobiotics), such as toxins or drugs, so that they can be removed from the body.",
    link: "https://tdcommons.ai/single_pred_tasks/adme/#cyp-p450-3a4-inhibition-veith-et-al"
  },

  // Toxicity
  AMES: {
    description: "Mutagenicity means the ability of a drug to induce genetic alterations. Drugs that can cause damage to the DNA can result in cell death or other severe adverse effects. The Ames test is a short-term bacterial reverse mutation assay detecting a large number of compounds which can induce genetic damage and frameshift mutations.",
    link: "https://tdcommons.ai/single_pred_tasks/tox/#ames-mutagenicity"
  },
  Carcinogens_Lagunin: {
    description: "A carcinogen is any substance, radionuclide, or radiation that promotes carcinogenesis, the formation of cancer. This may be due to the ability to damage the genome or to the disruption of cellular metabolic processes.",
    link: "https://tdcommons.ai/single_pred_tasks/tox/#carcinogens"
  },
  ClinTox: {
    description: "Clinical toxicity likelihood indicator.",
    link: "https://tdcommons.ai/single_pred_tasks/tox/#clintox"
  },
  DILI: {
    description: "Drug-induced liver injury (DILI) is fatal liver disease caused by drugs and it has been the single most frequent cause of safety-related drug marketing withdrawals (e.g. iproniazid, ticrynafen, benoxaprofen). ",
    link: "https://tdcommons.ai/single_pred_tasks/tox/#dili-drug-induced-liver-injury"
  },
  hERG: {
    description: "Human ether-à-go-go related gene (hERG) is crucial for the coordination of the heart's beating. Thus, if a drug blocks the hERG, it could lead to severe adverse effects. Therefore, reliable prediction of hERG liability in the early stages of drug design is quite important to reduce the risk of cardiotoxicity-related attritions in the later development stages.",
    link: "https://tdcommons.ai/single_pred_tasks/tox/#herg-blockers"
  },
  herg_central: {
    description: "Human ether-à-go-go related gene (hERG) is crucial for the coordination of the heart's beating. Thus, if a drug blocks the hERG, it could lead to severe adverse effects. Therefore, reliable prediction of hERG liability in the early stages of drug design is quite important to reduce the risk of cardiotoxicity-related attritions in the later development stages.",
    link: "https://tdcommons.ai/single_pred_tasks/tox/#herg-central"
  },
  hERG_Karim: {
    description: "A integrated Ether-a-go-go-related gene (hERG) dataset consisting of molecular structures labelled as hERG (<10uM) and non-hERG (>=10uM) blockers in the form of SMILES strings.",
    link: "https://tdcommons.ai/single_pred_tasks/tox/#herg-karim-et-al"
  },
  LD50_Zhu: {
    description: "Acute toxicity LD50 measures the most conservative dose that can lead to lethal adverse effects. The higher the dose, the more lethal of a drug.",
    link: "https://tdcommons.ai/single_pred_tasks/tox/#acute-toxicity-ld50"
  },
  Skin_Reaction: {
    description: "Repetitive exposure to a chemical agent can induce an immune reaction in inherently susceptible individuals that leads to skin sensitization.",
    link: "https://tdcommons.ai/single_pred_tasks/tox/#skin-reaction"
  },
  Tox21: {
    description: "Tox21 is a data challenge which contains qualitative toxicity measurements for 7,831 compounds on 12 different targets, such as nuclear receptors and stree response pathways.",
    link: "https://tdcommons.ai/single_pred_tasks/tox/#tox21"
  },
  Toxcast: {
    description: "ToxCast includes qualitative results of over 600 experiments on 8k compounds.",
    link: "https://tdcommons.ai/single_pred_tasks/tox/#toxcast"
  },
}