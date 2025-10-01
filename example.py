import joblib
m = joblib.load("BackEnd/Models/Absorption/Bioavailability_Ma.joblib")
print(type(m))
print(m.keys() if isinstance(m, dict) else m)
