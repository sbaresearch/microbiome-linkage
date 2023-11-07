import numpy as np
import pandas as pd
from sklearn.neighbors import BallTree


def DATA(maindata, data, cla):
    print("Comparing databases " + maindata + " and " + data)

    main_df = pd.read_csv(maindata, sep='\t', skipinitialspace=True)
    main_df = main_df.drop(["RANDSID"], axis=1)
    main_df = main_df.transpose()
    main_index = list(main_df.index.values)
    main_columns = list(main_df.columns.values)

    mult_df = pd.read_csv(data, sep='\t', skipinitialspace=True)
    mult_df = mult_df.drop(["RANDSID"], axis=1)
    mult_df = mult_df.transpose()
    mult_index = list(mult_df.index.values)
    mult_columns = list(mult_df.columns.values)

    main_df = pd.DataFrame(np.digitize(main_df, cla),
                           index=list(main_index),
                           columns=list(main_columns)).subtract(1)
    mult_df = pd.DataFrame(np.digitize(mult_df, cla),
                           index=list(mult_index),
                           columns=list(mult_columns)).subtract(1)

    selected_feat = np.nonzero(np.any(main_df != 0, axis=0).tolist())[0].tolist()
    print('Non-zero Features: ' + str(len(selected_feat)) + '/'
          + str(len(main_columns)))
    main_df = main_df.iloc[:, selected_feat]
    mult_df = mult_df.iloc[:, selected_feat]

    return main_df, mult_df


def DIST(main_df, mult_df, kn, balltreemeasure, thresholdmeasure, thresh):
    header_main = list(main_df.index.values)

    # Build Tree on Original
    bat = BallTree(main_df, metric=balltreemeasure)

    mult_neighborlist = []
    dynamic_thresh = []

    for i in range(len(mult_df)):
        pt = mult_df.values[i].reshape(1, -1)
        ind = bat.query(pt, k=kn)[1][0].tolist()
        dist = bat.query(pt, k=kn)[0][0].tolist()

        ids = [header_main[i] for i in ind]
        for d in range(len(dist)):
            d = round(len(main_df.columns.values) * d)
        mult_neighborlist.append([mult_df.index.values[i],
                                  ids, dist])

        thresholds = []
        measure_scores = []
        for n in range(kn):
            sample_scores = []
            for no in header_main:
                sample_scores.append(thresholdmeasure(
                    main_df.loc[mult_neighborlist[i][1][n]].tolist(),
                    main_df.loc[no].values.tolist())
                )

            sample_scores = np.sort(sample_scores)
            if thresh == 'Relaxed':
                thresholds.append(sample_scores[2])
            elif thresh == 'Strict':
                thresholds.append(sample_scores[1])
            elif thresh == 'None':
                thresholds.append(float("inf"))

            measure_scores.append(thresholdmeasure(
                main_df.loc[mult_neighborlist[i][1][n]].tolist(),
                mult_df.loc[mult_neighborlist[i][0]].tolist())
            )

        mult_neighborlist[i].append(measure_scores)
        dynamic_thresh.append(thresholds)

    tpcount = 0
    fpcount = 0
    tpfpcount = 0
    fncount = 0

    neighbourlist_idability_style = ""

    for i in range(len(mult_df)):
        for n in range(kn):
            index = kn - 1 - n
            if mult_neighborlist[i][-1][index] >= dynamic_thresh[i][index]:
                mult_neighborlist[i][1].remove(mult_neighborlist[i][1][index])

        if mult_neighborlist[i][0] not in mult_neighborlist[i][1]:
            if len(mult_neighborlist[i][1]) > 0:
                fpcount += 1
                neighbourlist_idability_style += mult_neighborlist[i][0] + "\tmatches_fp\t" + "\t".join(mult_neighborlist[i][1]) + "\n"
            else:
                fncount += 1
                neighbourlist_idability_style += mult_neighborlist[i][0] + "\tno_matches\n"
        else:
            if len(mult_neighborlist[i][1]) > 1:
                mult_neighborlist[i].append('Match')
                tpfpcount += 1
                neighbourlist_idability_style += mult_neighborlist[i][0] + "\tmatches_tpfp\t" + "\t".join(mult_neighborlist[i][1]) + "\n"
            else:
                mult_neighborlist[i].append('Match')
                tpcount += 1
                neighbourlist_idability_style += mult_neighborlist[i][0] + "\tmatches\t" + "\t".join(mult_neighborlist[i][1]) + "\n"

    neighbourlist_idability_style = "# 1|TP: " + str(tpcount) + "\n" + \
                                    "# 2|TP+FP: " + str(tpfpcount) + "\n" + \
                                    "# 3|FP: " + str(fpcount) + "\n" + \
                                    "# 4|FN: " + str(fncount) + "\n" + \
                                    "# 5|NA: 0\n" + neighbourlist_idability_style

    return tpcount, fpcount, tpfpcount, fncount, mult_neighborlist, neighbourlist_idability_style
