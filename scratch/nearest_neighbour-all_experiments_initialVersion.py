import os,sys
import prettytable
import numpy as np
from scipy.spatial import distance

from nearest_neighbour import DATA
from nearest_neighbour import DIST

########################################################################################################
# SETUP
file_path = os.path.dirname(os.path.abspath(__file__))
print ("Working directory path:" + file_path)

outputprefix = os.path.join(file_path + '/output/')
outputpathTables = os.path.join(outputprefix + 'tables/nearest_neighbour/')
outputpathLists = os.path.join(outputprefix + 'lists/nearest_neighbour/')


#PARAMETERS
k=1

taxon_cla=[-np.inf, 0.00005,0.0005,0.005,0.05,0.5, np.inf]
gene_cla=[-np.inf, 0.005,0.05,0.5,5,50, np.inf]

#def h_minkowski(x,y):
#    return distance.minkowski(x,y,1.5)

thresh='Strict'
bt_meas=distance.correlation
ballt_meas='correlation'
ts_meas=distance.correlation
thresh_meas='correlation_standard'

stop='no'

def STATS(tp,fp,n):
    
    prec=tp/(tp+fp)
    rec=tp/n
    fscore=2*(prec*rec)/(prec+rec)

    return prec, rec, fscore



########################################################################################################
########################################################################################################


# SPECIES

print("\n********* Processing feature '" + "Species" + "'")

an_main = os.path.join(file_path, 'data/Species/species-anterior_nares-visit1.pcl')
an_mult = os.path.join(file_path, 'data/Species/species-anterior_nares-visit2.pcl')
bm_main = os.path.join(file_path, 'data/Species/species-buccal_mucosa-visit1.pcl')
bm_mult = os.path.join(file_path, 'data/Species/species-buccal_mucosa-visit2.pcl')
pf_main = os.path.join(file_path, 'data/Species/species-posterior_fornix-visit1.pcl')
pf_mult = os.path.join(file_path, 'data/Species/species-posterior_fornix-visit2.pcl')
st_main = os.path.join(file_path, 'data/Species/species-stool-visit1.pcl')
st_mult = os.path.join(file_path, 'data/Species/species-stool-visit2.pcl')
sp_main = os.path.join(file_path, 'data/Species/species-supragingival_plaque-visit1.pcl')
sp_mult = os.path.join(file_path, 'data/Species/species-supragingival_plaque-visit2.pcl')
td_main = os.path.join(file_path, 'data/Species/species-tongue_dorsum-visit1.pcl')
td_mult = os.path.join(file_path, 'data/Species/species-tongue_dorsum-visit2.pcl')

an_main, an_mult = DATA(an_main, an_mult, taxon_cla)
bm_main, bm_mult = DATA(bm_main, bm_mult, taxon_cla)
td_main, td_mult = DATA(td_main, td_mult, taxon_cla)
sp_main, sp_mult = DATA(sp_main, sp_mult, taxon_cla)
st_main, st_mult = DATA(st_main, st_mult, taxon_cla)
pf_main, pf_mult = DATA(pf_main, pf_mult, taxon_cla)

an_tp, an_fp, an_tpfp, an_fn, an_list = DIST(an_main,an_mult,k,bt_meas,ts_meas,thresh)
bm_tp, bm_fp, bm_tpfp, bm_fn, bm_list = DIST(bm_main,bm_mult,k,bt_meas,ts_meas,thresh)
td_tp, td_fp, td_tpfp, td_fn, td_list = DIST(td_main,td_mult,k,bt_meas,ts_meas,thresh)
sp_tp, sp_fp, sp_tpfp, sp_fn, sp_list = DIST(sp_main,sp_mult,k,bt_meas,ts_meas,thresh)
st_tp, st_fp, st_tpfp, st_fn, st_list = DIST(st_main,st_mult,k,bt_meas,ts_meas,thresh)
pf_tp, pf_fp, pf_tpfp, pf_fn, pf_list = DIST(pf_main,pf_mult,k,bt_meas,ts_meas,thresh)

row1=["TP", an_tp, bm_tp, td_tp, sp_tp, st_tp, pf_tp]
row2=["FP", an_fp, bm_fp, td_fp, sp_fp, st_fp, pf_fp]
row3=["TP+FP", an_tpfp, bm_tpfp, td_tpfp, sp_tpfp, st_tpfp, pf_tpfp]
row4=["FN", an_fn, bm_fn, td_fn, sp_fn, st_fn, pf_fn]

x=prettytable.PrettyTable(["Species", 'Ant. Nares', 'Buccal Mucosa', 'Tongue Dorsum', 'Supra. Plaque', 'Stool', 'Post. Fornix'])
x.add_row(row1)
x.add_row(row2)
x.add_row(row3)
x.add_row(row4)

print(x)

table = os.path.join(outputpathTables, 'Species_' + str(ballt_meas) + '_' + str(thresh_meas) + '_' + str(thresh) + '.txt')
with open(table, 'w') as w:
    w.write(str(x))

sites_generic = ["anterior_nares", "buccal_mucosa", "posterior_fornix", "stool", "supragingival_plaque", "tongue_dorsum"]
lists = [an_list, bm_list, pf_list, st_list, sp_list, td_list]

for site in sites_generic:
    listFile = os.path.join(outputpathLists, 'Species_' + site + "_" + str(ballt_meas) + '_' + str(thresh_meas) + '_' + str(thresh) + '.txt')
    with open(listFile, 'w') as w:
        w.write(str(an_list))
        print("wrote matches to: " + listFile)

if stop == 'yes':
    sys.exit()
    
spec_an_prec, spec_an_rec, spec_an_fscore = STATS(an_tp, an_fp, an_tp+an_tpfp+an_fp+an_fn)
spec_bm_prec, spec_bm_rec, spec_bm_fscore = STATS(bm_tp, bm_fp, bm_tp+bm_tpfp+bm_fp+bm_fn)
spec_pf_prec, spec_pf_rec, spec_pf_fscore = STATS(pf_tp, pf_fp, pf_tp+pf_tpfp+pf_fp+pf_fn)
spec_st_prec, spec_st_rec, spec_st_fscore = STATS(st_tp, st_fp, st_tp+st_tpfp+st_fp+st_fn)
spec_sp_prec, spec_sp_rec, spec_sp_fscore = STATS(sp_tp, sp_fp, sp_tp+sp_tpfp+sp_fp+sp_fn)
spec_td_prec, spec_td_rec, spec_td_fscore = STATS(td_tp, td_fp, td_tp+td_tpfp+td_fp+td_fn)

########################################################################################################
########################################################################################################

# MARKERS

print("\n********* Processing feature '" + "Markers" + "'")

an_main = os.path.join(file_path, 'data/Markers/markers-anterior_nares-visit1.pcl')
an_mult = os.path.join(file_path, 'data/Markers/markers-anterior_nares-visit2.pcl')
bm_main = os.path.join(file_path, 'data/Markers/markers-buccal_mucosa-visit1.pcl')
bm_mult = os.path.join(file_path, 'data/Markers/markers-buccal_mucosa-visit2.pcl')
pf_main = os.path.join(file_path, 'data/Markers/markers-posterior_fornix-visit1.pcl')
pf_mult = os.path.join(file_path, 'data/Markers/markers-posterior_fornix-visit2.pcl')
st_main = os.path.join(file_path, 'data/Markers/markers-stool-visit1.pcl')
st_mult = os.path.join(file_path, 'data/Markers/markers-stool-visit2.pcl')
sp_main = os.path.join(file_path, 'data/Markers/markers-supragingival_plaque-visit1.pcl')
sp_mult = os.path.join(file_path, 'data/Markers/markers-supragingival_plaque-visit2.pcl')
td_main = os.path.join(file_path, 'data/Markers/markers-tongue_dorsum-visit1.pcl')
td_mult = os.path.join(file_path, 'data/Markers/markers-tongue_dorsum-visit2.pcl')

an_main, an_mult = DATA(an_main, an_mult, gene_cla)
bm_main, bm_mult = DATA(bm_main, bm_mult, gene_cla)
td_main, td_mult = DATA(td_main, td_mult, gene_cla)
sp_main, sp_mult = DATA(sp_main, sp_mult, gene_cla)
st_main, st_mult = DATA(st_main, st_mult, gene_cla)
pf_main, pf_mult = DATA(pf_main, pf_mult, gene_cla)

an_tp, an_fp, an_tpfp, an_fn, an_list = DIST(an_main,an_mult,k,bt_meas,ts_meas,thresh)
bm_tp, bm_fp, bm_tpfp, bm_fn, bm_list = DIST(bm_main,bm_mult,k,bt_meas,ts_meas,thresh)
td_tp, td_fp, td_tpfp, td_fn, td_list = DIST(td_main,td_mult,k,bt_meas,ts_meas,thresh)
sp_tp, sp_fp, sp_tpfp, sp_fn, sp_list = DIST(sp_main,sp_mult,k,bt_meas,ts_meas,thresh)
st_tp, st_fp, st_tpfp, st_fn, st_list = DIST(st_main,st_mult,k,bt_meas,ts_meas,thresh)
pf_tp, pf_fp, pf_tpfp, pf_fn, pf_list = DIST(pf_main,pf_mult,k,bt_meas,ts_meas,thresh)

row1=["TP", an_tp, bm_tp, td_tp, sp_tp, st_tp, pf_tp]
row2=["FP", an_fp, bm_fp, td_fp, sp_fp, st_fp, pf_fp]
row3=["TP+FP", an_tpfp, bm_tpfp, td_tpfp, sp_tpfp, st_tpfp, pf_tpfp]
row4=["FN", an_fn, bm_fn, td_fn, sp_fn, st_fn, pf_fn]

x=prettytable.PrettyTable(["Markers", 'Ant. Nares', 'Buccal Mucosa', 'Tongue Dorsum', 'Supra. Plaque', 'Stool', 'Post. Fornix'])
x.add_row(row1)
x.add_row(row2)
x.add_row(row3)
x.add_row(row4)

print(x)

table = os.path.join(outputpathTables, 'Markers_' + str(ballt_meas) + '_' + str(thresh_meas) + '_' + str(thresh) + '.txt')
with open(table, 'w') as w:
    w.write(str(x))


mark_an_prec, mark_an_rec, mark_an_fscore = STATS(an_tp, an_fp, an_tp+an_tpfp+an_fp+an_fn)
mark_bm_prec, mark_bm_rec, mark_bm_fscore = STATS(bm_tp, bm_fp, bm_tp+bm_tpfp+bm_fp+bm_fn)
mark_pf_prec, mark_pf_rec, mark_pf_fscore = STATS(pf_tp, pf_fp, pf_tp+pf_tpfp+pf_fp+pf_fn)
mark_st_prec, mark_st_rec, mark_st_fscore = STATS(st_tp, st_fp, st_tp+st_tpfp+st_fp+st_fn)
mark_sp_prec, mark_sp_rec, mark_sp_fscore = STATS(sp_tp, sp_fp, sp_tp+sp_tpfp+sp_fp+sp_fn)
mark_td_prec, mark_td_rec, mark_td_fscore = STATS(td_tp, td_fp, td_tp+td_tpfp+td_fp+td_fn)


########################################################################################################
########################################################################################################



# KBWINDOWS

print("\n********* Processing feature '" + "KBWindows" + "'")

an_main = os.path.join(file_path, 'data/KBWindows/kbwindows-anterior_nares-visit1.pcl')
an_mult = os.path.join(file_path, 'data/KBWindows/kbwindows-anterior_nares-visit2.pcl')
bm_main = os.path.join(file_path, 'data/KBWindows/kbwindows-buccal_mucosa-visit1.pcl')
bm_mult = os.path.join(file_path, 'data/KBWindows/kbwindows-buccal_mucosa-visit2.pcl')
pf_main = os.path.join(file_path, 'data/KBWindows/kbwindows-posterior_fornix-visit1.pcl')
pf_mult = os.path.join(file_path, 'data/KBWindows/kbwindows-posterior_fornix-visit2.pcl')
st_main = os.path.join(file_path, 'data/KBWindows/kbwindows-stool-visit1.pcl')
st_mult = os.path.join(file_path, 'data/KBWindows/kbwindows-stool-visit2.pcl')
sp_main = os.path.join(file_path, 'data/KBWindows/kbwindows-supragingival_plaque-visit1.pcl')
sp_mult = os.path.join(file_path, 'data/KBWindows/kbwindows-supragingival_plaque-visit2.pcl')
td_main = os.path.join(file_path, 'data/KBWindows/kbwindows-tongue_dorsum-visit1.pcl')
td_mult = os.path.join(file_path, 'data/KBWindows/kbwindows-tongue_dorsum-visit2.pcl')

an_main, an_mult = DATA(an_main, an_mult, gene_cla)
bm_main, bm_mult = DATA(bm_main, bm_mult, gene_cla)
td_main, td_mult = DATA(td_main, td_mult, gene_cla)
sp_main, sp_mult = DATA(sp_main, sp_mult, gene_cla)
st_main, st_mult = DATA(st_main, st_mult, gene_cla)
pf_main, pf_mult = DATA(pf_main, pf_mult, gene_cla)

an_tp, an_fp, an_tpfp, an_fn, an_list = DIST(an_main,an_mult,k,bt_meas,ts_meas,thresh)
bm_tp, bm_fp, bm_tpfp, bm_fn, bm_list = DIST(bm_main,bm_mult,k,bt_meas,ts_meas,thresh)
td_tp, td_fp, td_tpfp, td_fn, td_list = DIST(td_main,td_mult,k,bt_meas,ts_meas,thresh)
sp_tp, sp_fp, sp_tpfp, sp_fn, sp_list = DIST(sp_main,sp_mult,k,bt_meas,ts_meas,thresh)
st_tp, st_fp, st_tpfp, st_fn, st_list = DIST(st_main,st_mult,k,bt_meas,ts_meas,thresh)
pf_tp, pf_fp, pf_tpfp, pf_fn, pf_list = DIST(pf_main,pf_mult,k,bt_meas,ts_meas,thresh)

row1=["TP", an_tp, bm_tp, td_tp, sp_tp, st_tp, pf_tp]
row2=["FP", an_fp, bm_fp, td_fp, sp_fp, st_fp, pf_fp]
row3=["TP+FP", an_tpfp, bm_tpfp, td_tpfp, sp_tpfp, st_tpfp, pf_tpfp]
row4=["FN", an_fn, bm_fn, td_fn, sp_fn, st_fn, pf_fn]

x=prettytable.PrettyTable(["KBWindows", 'Ant. Nares', 'Buccal Mucosa', 'Tongue Dorsum', 'Supra. Plaque', 'Stool', 'Post. Fornix'])
x.add_row(row1)
x.add_row(row2)
x.add_row(row3)
x.add_row(row4)

print(x)

table = os.path.join(outputpathTables, 'KBWindows_' + str(ballt_meas) + '_' + str(thresh_meas) + '_' + str(thresh) + '.txt')
with open(table, 'w') as w:
    w.write(str(x))
    
kbw_an_prec, kbw_an_rec, kbw_an_fscore = STATS(an_tp, an_fp, an_tp+an_tpfp+an_fp+an_fn)
kbw_bm_prec, kbw_bm_rec, kbw_bm_fscore = STATS(bm_tp, bm_fp, bm_tp+bm_tpfp+bm_fp+bm_fn)
kbw_pf_prec, kbw_pf_rec, kbw_pf_fscore = STATS(pf_tp, pf_fp, pf_tp+pf_tpfp+pf_fp+pf_fn)
kbw_st_prec, kbw_st_rec, kbw_st_fscore = STATS(st_tp, st_fp, st_tp+st_tpfp+st_fp+st_fn)
kbw_sp_prec, kbw_sp_rec, kbw_sp_fscore = STATS(sp_tp, sp_fp, sp_tp+sp_tpfp+sp_fp+sp_fn)
kbw_td_prec, kbw_td_rec, kbw_td_fscore = STATS(td_tp, td_fp, td_tp+td_tpfp+td_fp+td_fn)

########################################################################################################
########################################################################################################


# OTU

print("\n********* Processing feature '" + "OTUs" + "'")

laf_main = os.path.join(file_path, 'data/OTUs/otus-l_antecubital_fossa-visit1.pcl')
laf_mult = os.path.join(file_path, 'data/OTUs/otus-l_antecubital_fossa-visit2.pcl')
raf_main = os.path.join(file_path, 'data/OTUs/otus-r_antecubital_fossa-visit1.pcl')
raf_mult = os.path.join(file_path, 'data/OTUs/otus-r_antecubital_fossa-visit2.pcl')
lrc_main = os.path.join(file_path, 'data/OTUs/otus-l_retroauricular_crease-visit1.pcl')
lrc_mult = os.path.join(file_path, 'data/OTUs/otus-l_retroauricular_crease-visit2.pcl')
rrc_main = os.path.join(file_path, 'data/OTUs/otus-r_retroauricular_crease-visit1.pcl')
rrc_mult = os.path.join(file_path, 'data/OTUs/otus-r_retroauricular_crease-visit2.pcl')
an_main = os.path.join(file_path, 'data/OTUs/otus-anterior_nares-visit1.pcl')
an_mult = os.path.join(file_path, 'data/OTUs/otus-anterior_nares-visit2.pcl')
bm_main = os.path.join(file_path, 'data/OTUs/otus-buccal_mucosa-visit1.pcl')
bm_mult = os.path.join(file_path, 'data/OTUs/otus-buccal_mucosa-visit2.pcl')

laf_main, laf_mult = DATA(laf_main, laf_mult, taxon_cla)
raf_main, raf_mult = DATA(raf_main, raf_mult, taxon_cla)
lrc_main, lrc_mult = DATA(lrc_main, lrc_mult, taxon_cla)
rrc_main, rrc_mult = DATA(rrc_main, rrc_mult, taxon_cla)
an_main, an_mult = DATA(an_main, an_mult, taxon_cla)
bm_main, bm_mult = DATA(bm_main, bm_mult, taxon_cla)

laf_tp, laf_fp, laf_tpfp, laf_fn, laf_list = DIST(laf_main, laf_mult, k, bt_meas, ts_meas, thresh)
raf_tp, raf_fp, raf_tpfp, raf_fn, raf_list = DIST(raf_main, raf_mult, k, bt_meas, ts_meas, thresh)
lrc_tp, lrc_fp, lrc_tpfp, lrc_fn, lrc_list = DIST(lrc_main, lrc_mult, k, bt_meas, ts_meas, thresh)
rrc_tp, rrc_fp, rrc_tpfp, rrc_fn, rrc_list = DIST(rrc_main, rrc_mult, k, bt_meas, ts_meas, thresh)
an_tp, an_fp, an_tpfp, an_fn, an_list = DIST(an_main, an_mult, k, bt_meas, ts_meas, thresh)
bm_tp, bm_fp, bm_tpfp, bm_fn, bm_list = DIST(bm_main, bm_mult, k, bt_meas, ts_meas, thresh)

row1 = ["TP", laf_tp, raf_tp, lrc_tp, rrc_tp, an_tp, bm_tp]
row2 = ["FP", laf_fp, raf_fp, lrc_fp, rrc_fp, an_fp, bm_fp]
row3 = ["TP+FP", laf_tpfp, raf_tpfp, lrc_tpfp, rrc_tpfp, an_tpfp, bm_tpfp]
row4 = ["FN", laf_fn, raf_fn, lrc_fn, rrc_fn, an_fn, bm_fn]

x = prettytable.PrettyTable(
    ["OTU1", 'L Ante. Fossa', 'R Ante. Fossa', 'L Retro. Crease', 'R Retro. Crease', 'Ant. Nares', 'Buccal Mucosa'])
x.add_row(row1)
x.add_row(row2)
x.add_row(row3)
x.add_row(row4)

print(x)

table = os.path.join(outputpathTables, 'OTU1_' + str(ballt_meas) + '_' + str(thresh_meas) + '_' + str(thresh) + '.txt')
with open(table, 'w') as w:
    w.write(str(x))

########################################################################################################

td_main = os.path.join(file_path, 'data/OTUs/otus-tongue_dorsum-visit1.pcl')
td_mult = os.path.join(file_path, 'data/OTUs/otus-tongue_dorsum-visit2.pcl')
hp_main = os.path.join(file_path, 'data/OTUs/otus-hard_palate-visit1.pcl')
hp_mult = os.path.join(file_path, 'data/OTUs/otus-hard_palate-visit2.pcl')
sv_main = os.path.join(file_path, 'data/OTUs/otus-saliva-visit1.pcl')
sv_mult = os.path.join(file_path, 'data/OTUs/otus-saliva-visit2.pcl')
th_main = os.path.join(file_path, 'data/OTUs/otus-throat-visit1.pcl')
th_mult = os.path.join(file_path, 'data/OTUs/otus-throat-visit2.pcl')
pt_main = os.path.join(file_path, 'data/OTUs/otus-palatine_tonsils-visit1.pcl')
pt_mult = os.path.join(file_path, 'data/OTUs/otus-palatine_tonsils-visit2.pcl')
sp_main = os.path.join(file_path, 'data/OTUs/otus-supragingival_plaque-visit1.pcl')
sp_mult = os.path.join(file_path, 'data/OTUs/otus-supragingival_plaque-visit2.pcl')

td_main, td_mult = DATA(td_main, td_mult, taxon_cla)
hp_main, hp_mult = DATA(hp_main, hp_mult, taxon_cla)
sv_main, sv_mult = DATA(sv_main, sv_mult, taxon_cla)
th_main, th_mult = DATA(th_main, th_mult, taxon_cla)
pt_main, pt_mult = DATA(pt_main, pt_mult, taxon_cla)
sp_main, sp_mult = DATA(sp_main, sp_mult, taxon_cla)

td_tp, td_fp, td_tpfp, td_fn, td_list = DIST(td_main, td_mult, k, bt_meas, ts_meas, thresh)
hp_tp, hp_fp, hp_tpfp, hp_fn, hp_list = DIST(hp_main, hp_mult, k, bt_meas, ts_meas, thresh)
sv_tp, sv_fp, sv_tpfp, sv_fn, sv_list = DIST(sv_main, sv_mult, k, bt_meas, ts_meas, thresh)
th_tp, th_fp, th_tpfp, th_fn, th_list = DIST(th_main, th_mult, k, bt_meas, ts_meas, thresh)
pt_tp, pt_fp, pt_tpfp, pt_fn, pt_list = DIST(pt_main, pt_mult, k, bt_meas, ts_meas, thresh)
sp_tp, sp_fp, sp_tpfp, sp_fn, sp_list = DIST(sp_main, sp_mult, k, bt_meas, ts_meas, thresh)

row1 = ["TP", td_tp, hp_tp, sv_tp, th_tp, pt_tp, sp_tp]
row2 = ["FP", td_fp, hp_fp, sv_fp, th_fp, pt_fp, sp_fp]
row3 = ["TP+FP", td_tpfp, hp_tpfp, sv_tpfp, th_tpfp, pt_tpfp, sp_tpfp]
row4 = ["FN", td_fn, hp_fn, sv_fn, th_fn, pt_fn, sp_fn]

x = prettytable.PrettyTable(
    ["OTU2", 'Tongue Dorsum', 'Hard Palate', 'Saliva', 'Throat', 'Palatine Tonsils', 'Sprag. Plaque'])
x.add_row(row1)
x.add_row(row2)
x.add_row(row3)
x.add_row(row4)

print(x)

table = os.path.join(outputpathTables, 'OTU2_' + str(ballt_meas) + '_' + str(thresh_meas) + '_' + str(thresh) + '.txt')
with open(table, 'w') as w:
    w.write(str(x))

########################################################################################################

sb_main = os.path.join(file_path, 'data/OTUs/otus-subgingival_plaque-visit1.pcl')
sb_mult = os.path.join(file_path, 'data/OTUs/otus-subgingival_plaque-visit2.pcl')
kt_main = os.path.join(file_path, 'data/OTUs/otus-keratinized_gingiva-visit1.pcl')
kt_mult = os.path.join(file_path, 'data/OTUs/otus-keratinized_gingiva-visit2.pcl')
st_main = os.path.join(file_path, 'data/OTUs/otus-stool-visit1.pcl')
st_mult = os.path.join(file_path, 'data/OTUs/otus-stool-visit2.pcl')
vi_main = os.path.join(file_path, 'data/OTUs/otus-vaginal_introitus-visit1.pcl')
vi_mult = os.path.join(file_path, 'data/OTUs/otus-vaginal_introitus-visit2.pcl')
mv_main = os.path.join(file_path, 'data/OTUs/otus-mid_vagina-visit1.pcl')
mv_mult = os.path.join(file_path, 'data/OTUs/otus-mid_vagina-visit2.pcl')
pf_main = os.path.join(file_path, 'data/OTUs/otus-posterior_fornix-visit1.pcl')
pf_mult = os.path.join(file_path, 'data/OTUs/otus-posterior_fornix-visit2.pcl')

sb_main, sb_mult = DATA(sb_main, sb_mult, taxon_cla)
kt_main, kt_mult = DATA(kt_main, kt_mult, taxon_cla)
st_main, st_mult = DATA(st_main, st_mult, taxon_cla)
vi_main, vi_mult = DATA(vi_main, vi_mult, taxon_cla)
mv_main, mv_mult = DATA(mv_main, mv_mult, taxon_cla)
pf_main, pf_mult = DATA(pf_main, pf_mult, taxon_cla)

sb_tp, sb_fp, sb_tpfp, sb_fn, sb_list = DIST(sb_main, sb_mult, k, bt_meas, ts_meas, thresh)
kt_tp, kt_fp, kt_tpfp, kt_fn, kt_list = DIST(kt_main, kt_mult, k, bt_meas, ts_meas, thresh)
st_tp, st_fp, st_tpfp, st_fn, st_list = DIST(st_main, st_mult, k, bt_meas, ts_meas, thresh)
vi_tp, vi_fp, vi_tpfp, vi_fn, vi_list = DIST(vi_main, vi_mult, k, bt_meas, ts_meas, thresh)
mv_tp, mv_fp, mv_tpfp, mv_fn, mv_list = DIST(mv_main, mv_mult, k, bt_meas, ts_meas, thresh)
pf_tp, pf_fp, pf_tpfp, pf_fn, pf_list = DIST(pf_main, pf_mult, k, bt_meas, ts_meas, thresh)

row1 = ["TP", sb_tp, kt_tp, st_tp, vi_tp, mv_tp, pf_tp]
row2 = ["FP", sb_fp, kt_fp, st_fp, vi_fp, mv_fp, pf_fp]
row3 = ["TP+FP", sb_tpfp, kt_tpfp, st_tpfp, vi_tpfp, mv_tpfp, pf_tpfp]
row4 = ["FN", sb_fn, kt_fn, st_fn, vi_fn, mv_fn, pf_fn]

x = prettytable.PrettyTable(
    ["OTU3", 'Subg. Plaque', 'Kerat. Gingiva', 'Stool', 'Vaginal Introitus', 'Mid Vagina', 'Post. Fornix'])
x.add_row(row1)
x.add_row(row2)
x.add_row(row3)
x.add_row(row4)

print(x)

table = os.path.join(outputpathTables, 'OTU3_' + str(ballt_meas) + '_' + str(thresh_meas) + '_' + str(thresh) + '.txt')
with open(table, 'w') as w:
    w.write(str(x))

otu_an_prec, otu_an_rec, otu_an_fscore = STATS(an_tp, an_fp, an_tp + an_tpfp + an_fp + an_fn)
otu_bm_prec, otu_bm_rec, otu_bm_fscore = STATS(bm_tp, bm_fp, bm_tp + bm_tpfp + bm_fp + bm_fn)
otu_pf_prec, otu_pf_rec, otu_pf_fscore = STATS(pf_tp, pf_fp, pf_tp + pf_tpfp + pf_fp + pf_fn)
otu_st_prec, otu_st_rec, otu_st_fscore = STATS(st_tp, st_fp, st_tp + st_tpfp + st_fp + st_fn)
otu_sp_prec, otu_sp_rec, otu_sp_fscore = STATS(sp_tp, sp_fp, sp_tp + sp_tpfp + sp_fp + sp_fn)
otu_td_prec, otu_td_rec, otu_td_fscore = STATS(td_tp, td_fp, td_tp + td_tpfp + td_fp + td_fn)

########################################################################################################
########################################################################################################


row1=["MeanPrec", np.round(np.mean([spec_an_prec,otu_an_prec,mark_an_prec,kbw_an_prec]),3),
      np.round(np.mean([spec_bm_prec,otu_bm_prec,mark_bm_prec,kbw_bm_prec]),3),
      np.round(np.mean([spec_td_prec,otu_td_prec,mark_td_prec,kbw_td_prec]),3),
      np.round(np.mean([spec_sp_prec,otu_sp_prec,mark_sp_prec,kbw_sp_prec]),3),
      np.round(np.mean([spec_st_prec,otu_st_prec,mark_st_prec,kbw_st_prec]),3),
      np.round(np.mean([spec_pf_prec,otu_pf_prec,mark_pf_prec,kbw_pf_prec]),3)]
row2=["MeanRec", np.round(np.mean([spec_an_rec,otu_an_rec,mark_an_rec,kbw_an_rec]),3),
      np.round(np.mean([spec_bm_rec,otu_bm_rec,mark_bm_rec,kbw_bm_rec]),3),
      np.round(np.mean([spec_td_rec,otu_td_rec,mark_td_rec,kbw_td_rec]),3),
      np.round(np.mean([spec_sp_rec,otu_sp_rec,mark_sp_rec,kbw_sp_rec]),3),
      np.round(np.mean([spec_st_rec,otu_st_rec,mark_st_rec,kbw_st_rec]),3),
      np.round(np.mean([spec_pf_rec,otu_pf_rec,mark_pf_rec,kbw_pf_rec]),3)]
row3=["MeanFScore", np.round(np.mean([spec_an_fscore,otu_an_fscore,mark_an_fscore,kbw_an_fscore]),3),
      np.round(np.mean([spec_bm_fscore,otu_bm_fscore,mark_bm_fscore,kbw_bm_fscore]),3),
      np.round(np.mean([spec_td_fscore,otu_td_fscore,mark_td_fscore,kbw_td_fscore]),3),
      np.round(np.mean([spec_sp_fscore,otu_sp_fscore,mark_sp_fscore,kbw_sp_fscore]),3),
      np.round(np.mean([spec_st_fscore,otu_st_fscore,mark_st_fscore,kbw_st_fscore]),3),
      np.round(np.mean([spec_pf_fscore,otu_pf_fscore,mark_pf_fscore,kbw_pf_fscore]),3)]

x=prettytable.PrettyTable(["Summary", 'Ant. Nares', 'Buccal Mucosa', 'Tongue Dorsum', 'Supra. Plaque', 'Stool', 'Post. Fornix'])
x.add_row(row1)
x.add_row(row2)
x.add_row(row3)

print(x)

table = os.path.join(outputpathTables, 'Summary_' + str(ballt_meas) + '_' + str(thresh_meas) + '_' + str(thresh) + '.txt')
with open(table, 'w') as w:
    w.write(str(x))




