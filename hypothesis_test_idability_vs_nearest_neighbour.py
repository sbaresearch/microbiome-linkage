import os
from collections import OrderedDict

from statsmodels.stats.contingency_tables import mcnemar
import shlex

# SETUP paths and constants for paths
file_path = os.path.dirname(os.path.abspath(__file__))
print("Working directory path:" + file_path)

thresh = 'Strict'
ballt_meas = 'correlation'
thresh_meas = 'correlation_standard'

# The feature types
features = ["Species", "OTUs", "Markers", "KBWindows"]

# Body sites are the same for everything, except OTUs, which have some in addition
sites_generic = ["anterior_nares", "buccal_mucosa", "tongue_dorsum", "supragingival_plaque", "stool", "posterior_fornix"]
sites_otus = ["anterior_nares", "buccal_mucosa", "hard_palate", "keratinized_gingiva", "l_antecubital_fossa", "l_retroauricular_crease", "mid_vagina", "palatine_tonsils", "posterior_fornix",
              "r_antecubital_fossa", "r_retroauricular_crease", "saliva", "stool", "subgingival_plaque", "supragingival_plaque", "throat", "tongue_dorsum", "vaginal_introitus"]
sites_otus = ["l_antecubital_fossa", "r_antecubital_fossa", "l_retroauricular_crease", "r_retroauricular_crease", "anterior_nares", "buccal_mucosa",
              "tongue_dorsum", "hard_palate", "saliva", "throat", "palatine_tonsils", "supragingival_plaque",
              "subgingival_plaque", "keratinized_gingiva", "stool", "vaginal_introitus", "mid_vagina", "posterior_fornix"]

debugInfo = True

# Iterate over all feature types
for feature in features:
    # check which body sites to iterate over (OTU has more)
    if feature == "OTUs":
        sites_all = sites_otus
    else:
        sites_all = sites_generic

    print("\n********* Processing feature '" + feature + "' on sites: " + str(sites_all))

    # Iterate over all body sites
    for site_base in sites_all:
        print("--- Processing site '" + site_base + "'")

        # read files for idability & nearest neighbour lists (need to first run nearest_neighbour-all_experiments and idability-all_experiments!)
        site = feature + "-" + site_base
        site_path = feature + "-" + site_base
        input_idability = "output/lists/idability/" + site_path + "-visit2." + feature + "-" + site_base + "-visit1.hits.txt"
        input_nn = "output/lists/nearest_neighbour/" + feature + '_' + site_base + "_" + str(ballt_meas) + '_' + str(thresh_meas) + '_' + str(thresh) + '.txt'

        # variables for scores / results
        scores_idability_optimistic = {}
        scores_idability_pessimistic = {}
        scores_nn = {}
        tp = fp = fn = tn = tpfp = fpfn = na = 0

        # Read IDAbility scores
        # We account for two cases of how to solve TP+FP: optimistic (every TP+FP ==> TP) and pessimisitc (every TP+FP => FP)
        # FIXME: does not account for TN yet
        with open(input_idability) as f:
            for line in f:
                fields = shlex.split(line, comments=True)
                # ignore empty lines
                if not fields:
                    continue

                targetRecord = fields[0] # the record to match up
                status = fields[1] # the matching status
                numMatches = len(fields) - 2 # the target record(s)
                if status == "matches": # if it was a match
                    targets = fields[2:]
                    if targetRecord in targets: # check if there was a match to the target record to match
                        if numMatches <= 1: # TP if there is only ONE hit
                            tp += 1
                            scores_idability_pessimistic[targetRecord] = True # TP
                            scores_idability_optimistic[targetRecord] = True # TP
                        else: # TP + FP if there are more than one hits
                            tpfp += 1
                            scores_idability_pessimistic[targetRecord] = False # TP+FP
                            scores_idability_optimistic[targetRecord] = True # TP+FP
                    else: # False Positives if we had matches but not to the target
                        scores_idability_optimistic[targetRecord] = False
                        scores_idability_pessimistic[targetRecord] = False
                        if numMatches <= 1: # FP if there is only one wrong hit
                            fp += 1
                        else:
                            fpfn += 1 # FN + FP otherwise
                else: # FN, NAs, ...
                    scores_idability_optimistic[targetRecord] = False
                    scores_idability_pessimistic[targetRecord] = False
                    if status == "no_code": # distinguish the two cases of no matches
                        na += 1
                    elif status == "no_matches":
                        fn += 1

        # Now the same for the nearest neighbour method
        # FIXME: does not account for TN yet
        with open(input_nn) as f:
            for line in f:
                fields = shlex.split(line, comments=True)
                if not fields:
                    continue
                targetRecord = fields[0]
                status = fields[1]
                numMatches = len(fields) - 2 # FIXME: does not work for the k>1 case yet!
                scores_nn[targetRecord] = True if status == "matches" else False #it is an easier setting here, as we are more explicit in the result file already

        # Sort the dictionaries for nicer readability in print output
        scores_nn = dict(sorted(scores_nn.items()))
        scores_idability_optimistic = dict(sorted(scores_idability_optimistic.items()))
        scores_idability_pessimistic = dict(sorted(scores_idability_pessimistic.items()))

        # Print some status
        if (debugInfo):
            print("IDability status: TP=" + str(tp) + " | FP=" + str(fp) + " | TPFP=" + str(tpfp) + " | FN=" + str(fn) + " | TN=" + str(tn) + " | FPFN=" + str(fpfn) + " | NA=" + str(na))
            print("IDability scores pessimistic, True+False=Total: " + str(sum(x == True for x in scores_idability_pessimistic.values())) + " + "  + str(sum(x == False for x in scores_idability_pessimistic.values())) + " = " + str(len(scores_idability_pessimistic) )+ "\t\t" + str(scores_idability_pessimistic))
            print("Nearest Neighbor scores, True+False=Total: " + str(sum(x == True for x in scores_nn.values())) + " + " + str(sum(x == False for x in scores_nn.values())) + " = " + str(len(scores_nn)) + "\t\t\t" + str(scores_nn))
            print("IDability scores optimistic, True+False=Total: " + str(sum(x == True for x in scores_idability_optimistic.values())) + " + "  + str(sum(x == False for x in scores_idability_optimistic.values())) + " = " + str(len(scores_idability_optimistic)) +  "\t\t" + str(scores_idability_optimistic))

        # Compute significance to both modes (optimistic and pessimistic)
        for mode in ["pessimistic", "optimistic"]:
            if mode == "optimistic":
                scores_to_compare = scores_idability_optimistic
            elif mode == "pessimistic":
                scores_to_compare = scores_idability_pessimistic
            else:
                print("Unsupported mode '" + mode + "'")
                continue

            # compute the contingency matrix for McNemar's test (see e.g. https://en.wikipedia.org/wiki/McNemar%27s_test)
            n00 = n01 = n10 = n11 = 0
            for key in scores_to_compare.keys():
                result_idability = scores_to_compare[key]
                result_nn = scores_nn[key]
                if result_idability and result_nn:
                    n11 += 1
                elif not result_idability and result_nn:
                    n01 += 1
                elif result_idability and not result_nn:
                    n10 += 1
                else:
                    n00 += 1

            # McNemar's test has two ways to compute the statistic, depending on the number of samples.
            exactTest = True if n01 < 25 and n10 < 25 else False
            print("\tComparing to " + mode + " idability scores, using exact mode: " + str(exactTest))

            contingency_table = [[n00, n01], [n10, n11]]
            result = mcnemar(contingency_table, exact=exactTest, correction=True)

            test_result = '\t\t' + str(contingency_table) + '==> statistic=%.3f, p-value=%.3f' % (result.statistic, result.pvalue)

            alpha = 0.05
            if result.pvalue > alpha:
                print(test_result + ' ==> same proportions of errors (fail to reject H0)')
            else:
                print(test_result + ' ==> different proportions of errors (reject H0)')
