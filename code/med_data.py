from pandas import read_csv, unique

raw_data_dir = "/data/penn_new/vitals_labs_anon"
Vitals = read_csv("%s/Vitals.csv"%(raw_data_dir))
Diagnosis = read_csv("%s/Dx.csv"%(raw_data_dir))
Encounters = read_csv("%s/Encounter.csv"%(raw_data_dir))
ENCOUNTER_COLUMNS = Encounters.columns.tolist()
PANS = unique(Vitals['PAN'])
