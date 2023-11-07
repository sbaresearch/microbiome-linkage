import os, sys
import prettytable
import numpy as np
from scipy.spatial import distance

from nearest_neighbour import DATA
from nearest_neighbour import DIST

########################################################################################################
# SETUP
file_path = os.path.dirname(os.path.abspath(__file__))
print("Working directory path:" + file_path)

outputprefix = os.path.join(file_path + '/output/')
outputpathTables = os.path.join(outputprefix + 'tables/nearest_neighbour/')
outputpathLists = os.path.join(outputprefix + 'lists/nearest_neighbour/')

# PARAMETERS
k = 1

# Abundance thresholds
taxon_cla = [-np.inf, 0.00005, 0.0005, 0.005, 0.05, 0.5, np.inf]
gene_cla = [-np.inf, 0.005, 0.05, 0.5, 5, 50, np.inf]

# Distance metrics
thresh='Strict'
bt_meas=distance.correlation
ballt_meas='correlation'
ts_meas=distance.correlation
thresh_meas='correlation_standard'

# Computes precision, recall, F1
def STATS(tp, fp, n):
    prec = tp / (tp + fp)
    rec = tp / n
    fscore = 2 * (prec * rec) / (prec + rec)

    return prec, rec, fscore


# The feature sets to iterate over
features = ["Species", "Markers", "KBWindows", "OTUs"]

# the sites to iterate over, plus their pretty names to be used for tables
sites_generic = [["anterior_nares", "buccal_mucosa", "tongue_dorsum", "supragingival_plaque", "stool", "posterior_fornix"]]
sites_generic_prettyNames = [['Ant. Nares', 'Buccal Mucosa', 'Tongue Dorsum', 'Supra. Plaque', 'Stool', 'Post. Fornix']]

# for OTUs, we have more sites; we group them into 3 lists, as we put them into three tables
sites_otus = [["l_antecubital_fossa", "r_antecubital_fossa", "l_retroauricular_crease", "r_retroauricular_crease", "anterior_nares", "buccal_mucosa"],
              ["tongue_dorsum", "hard_palate", "saliva", "throat", "palatine_tonsils", "supragingival_plaque"],
              ["subgingival_plaque", "keratinized_gingiva", "stool", "vaginal_introitus", "mid_vagina", "posterior_fornix"]]
sites_otus_prettyNames = [['L Ante. Fossa', 'R Ante. Fossa', 'L Retro. Crease', 'R Retro. Crease', 'Ant. Nares', 'Buccal Mucosa'],
                          ['Tongue Dorsum', 'Hard Palate', 'Saliva', 'Throat', 'Palatine Tonsils', 'Sprag. Plaque'],
                          ['Subg. Plaque', 'Kerat. Gingiva', 'Stool', 'Vaginal Introitus', 'Mid Vagina', 'Post. Fornix']]

########################################################################################################
# Iterate over sites, compute matches & stats

# Initialise dictionary
stats = {}
for subsite in sites_otus:
    for site in subsite:
        stats[site] = {}

# Process all features
for feature in features:
    # for OTUs, we use different body sites --> check which ones to use, and set as the current working version
    if (feature == "OTUs"):
        current_sites_group = sites_otus
        current_sites_group_prettyNames = sites_otus_prettyNames
    else:
        current_sites_group = sites_generic
        current_sites_group_prettyNames = sites_generic_prettyNames

    print("\n********* Processing feature '" + feature + "'")
    # We need different feature abundance limits for taxon vs. gene leve features --> select which ones
    if (feature in ["Markers", "KBWindows"]):
        cla = gene_cla
    else:
        cla = taxon_cla

    # Process each of the available body sites of the current feature type
    for index in range(len(current_sites_group)):
        results = {}

        subsite = current_sites_group[index]
        subsite_prettyNames = current_sites_group_prettyNames[index]
        print("Processing subsite '" + str(subsite) + "' of " + str(current_sites_group))
        siteName = feature + (str(index + 1) if len(current_sites_group) > 1 else "")

        for site in subsite:
            print("\tProcessing site '" + site + "'")

            # load data
            firstVisit = os.path.join(file_path, 'data/' + feature + "/" + feature.lower() + "-" + site + "-visit1.pcl")
            secondVisit = os.path.join(file_path, 'data/' + feature + "/" + feature.lower() + "-" + site + "-visit2.pcl")
            dataFirst, dataSecond = DATA(firstVisit, secondVisit, cla)

            # compute distances / matches
            results[site] = DIST(dataFirst, dataSecond, k, bt_meas, ts_meas, thresh)

            # write list of matches to file
            listFile = os.path.join(outputpathLists, feature + '_' + site + "_" + str(ballt_meas) + '_' + str(thresh_meas) + '_' + str(thresh) + '.txt')
            with open(listFile, 'w') as w:
                w.write(str(results[site][5]))
                print("Wrote matches to: " + listFile)

            # compute overall stats from this site
            stats[site][feature] = STATS (results[site][0], results[site][1], results[site][0] + results[site][1] + results[site][2] + results[site][3])


        # table with stats for this group of sites
        x = prettytable.PrettyTable([siteName] + subsite_prettyNames)
        x.add_row(["TP"] + [results[key][0] for key in results.keys()])
        x.add_row(["FP"] + [results[key][1] for key in results.keys()])
        x.add_row(["TP+FP"] + [results[key][2] for key in results.keys()])
        x.add_row(["FN"] + [results[key][3] for key in results.keys()])
        print(x)

        # write table to file
        table = os.path.join(outputpathTables, siteName + '_' + str(ballt_meas) + '_' + str(thresh_meas) + '_' + str(thresh) + '.txt')
        with open(table, 'w') as w:
            w.write(str(x))

########################################################################################################
# Compute overall stats

metrics = ["MeanPrec", "MeanRec", "MeanFScore"]
x = prettytable.PrettyTable(["Summary", 'Ant. Nares', 'Buccal Mucosa', 'Tongue Dorsum', 'Supra. Plaque', 'Stool', 'Post. Fornix'])

current_sites_group = sites_generic[0]

# compute all defined metrics as averages from the individual metrics
for i in range(len(metrics)):
    metric = metrics[i]
    row = [metric,
         np.round(np.mean([stats[current_sites_group[0]][feature][i] for feature in features]), 3),
         np.round(np.mean([stats[current_sites_group[1]][feature][i] for feature in features]), 3),
         np.round(np.mean([stats[current_sites_group[2]][feature][i] for feature in features]), 3),
         np.round(np.mean([stats[current_sites_group[3]][feature][i] for feature in features]), 3),
         np.round(np.mean([stats[current_sites_group[4]][feature][i] for feature in features]), 3),
         np.round(np.mean([stats[current_sites_group[5]][feature][i] for feature in features]), 3),
         ]
    x.add_row(row)
print(x)

# write table to file
table = os.path.join(outputpathTables, 'Summary_' + str(ballt_meas) + '_' + str(thresh_meas) + '_' + str(thresh) + '.txt')
with open(table, 'w') as w:
    w.write(str(x))
