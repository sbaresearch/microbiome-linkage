#! /usr/bin/env python3

"""
Runs idability.py on all available feature types and body-sites
"""

import os
import subprocess

## SETUP paths and constants for calling of the python script
file_path = os.path.dirname(os.path.abspath(__file__))
print("Working directory path:" + file_path)

mainScript = ["/usr/bin/python3", file_path + "/idability.py"]
print("Script static part: " + ' '.join(mainScript))

# The feature types
features = ["KBWindows", "Markers", "Species", "OTUs"]

# Body sites are the same for everything, except OTUs, which have some in addition
sites_generic = ["anterior_nares", "buccal_mucosa", "posterior_fornix", "stool", "supragingival_plaque", "tongue_dorsum"]
sites_otus = ["anterior_nares", "buccal_mucosa", "hard_palate", "keratinized_gingiva", "l_antecubital_fossa", "l_retroauricular_crease", "mid_vagina", "palatine_tonsils", "posterior_fornix",
              "r_antecubital_fossa", "r_retroauricular_crease", "saliva", "stool", "subgingival_plaque", "supragingival_plaque", "throat", "tongue_dorsum", "vaginal_introitus"]

# Iterate over all feature types
for feature in features:
    # check which body sites to iterate over
    if feature == "OTUs":
        sites_all = sites_generic + sites_otus
    else:
        sites_all = sites_generic

    # KBWindows and Markers features are in "reads per kilobase per million reads (RPKM)", while Species and OTUs are with "relative abundance ((relab)
    if feature == "KBWindows" or feature == "Markers":
        metamode = "rpkm"
    else:
        metamode = "relab"

    print("\n********* Processing feature '" + feature + "' with meta_mode '" + metamode + "' on sites: " + str(sites_all))

    # Iterate over all body sites
    for site_base in sites_all:
        site = feature.lower() + "-" + site_base

        # construct codes from first visit
        input_visit1 = "data/" + feature + "/" + site + "-visit1.pcl"
        output_visit1 = "output/idability-codes/" + feature + "-" + site_base + "-visit1.codes.txt"
        visit1Args = ["--meta_mode", metamode, input_visit1, "--output", output_visit1]
        print("Running: " + ' '.join(mainScript + visit1Args))
        subprocess.run(mainScript + visit1Args)

        # now construct codes and compare to second visit
        input_visit2 = "data/" + feature + "/" + site + "-visit2.pcl"
        output_matches = "output/lists/idability/" + feature + "-" + site_base + "-visit2." + feature + "-" + site_base + "-visit1.hits.txt"
        visit2Args = [input_visit2, "--codes", output_visit1, "--meta_mode", metamode, "--output", output_matches]
        print("Running: " + ' '.join(mainScript + visit2Args))
        subprocess.run(mainScript + visit2Args)
