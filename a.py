import sys
from base64 import b64decode
from array import array
from model import msms
import pandas as pd
import numpy as np
from model1 import ppm_error, gaussian_similarity, cal_theory_masses, getd3, get_peaks_from_xml,plotMsMs, addcolor, draw_plot

import matplotlib.pyplot as plt

sigma = 0.4
filename = sys.argv[1] # 17mix_test2.mzxml
num = int(sys.argv[2])  # 1298
aaseq = sys.argv[3]  # TYDSYLGDDYVR

peakselt = msms()
peakselt.text = get_peaks_from_xml("17mix_test2.mzxml", num).text
l = len(aaseq)
b = []
y = []
for i in range(l-1):
    b.append(aaseq[:i+1])
    y.append(aaseq[i+1:])
d3 = getd3()
    

massb = cal_theory_masses(b,d3,1)
massy = cal_theory_masses(y,d3,19)
# print(massb, len(massb),l,massy, len(massy),len(b),len(y))

d2 = {}
d2["b"] = b
d2["b_index"] = [x[1] for x in massb]
d2["mwb"] = [x[0] for x in massb]
d2["y"] = y
d2["y_index"] = [x[1] for x in massy]
d2["mwy"] = [x[0] for x in massy]
df2 = pd.DataFrame(data=d2)
print(df2)

peaks = array('f',b64decode(peakselt.text))
if sys.byteorder != 'big':
    peaks.byteswap()
    
mzs = peaks[::2]
ints = peaks[1::2]
relative_abundance = [ints_i/max(ints)*100 for ints_i in ints]

d={}
d["mzs"] = mzs
d["ints"] = ints
d["relative_abundance"] = relative_abundance
d["gs"] = [""]*len(mzs)
d["gaussian_sim_match"] = [""]*len(mzs)
d["color_gs"] = [""]*len(mzs)

d["ppm"] = [""]*len(mzs)
d["ppm_match"] = [""]*len(mzs)
d["color"] = [""]*len(mzs)

df = pd.DataFrame(data=d)

# calculate gaussian similarity
for index, mzs_i in enumerate(mzs):
    sim = []
    if df.loc[index, "relative_abundance"]<5:
            continue
    for mw_i in massy:    
        sim.append((gaussian_similarity(mw_i[0], mzs_i, sigma=0.4), mw_i[1]))
    index2 = max(sim, key=lambda x: x[0])
    if max(sim, key=lambda x: x[0])[0] < 0.999:
        sim2=[]
        for mw_i in massb:
            sim2.append((gaussian_similarity(mw_i[0], mzs_i, sigma=0.4),mw_i[1]))
        if max(sim2, key=lambda x: x[0])[0] < 0.999:
            continue
        else:
            index2 = max(sim2, key=lambda x: x[0])
    df.loc[index, "gaussian_sim_match"] = index2[1]
    df.loc[index, "gs"] = index2[0]
    # print(f"Found matching b-ion: {mwb_i} ~ {mzs_i}")
            
            

# calculate ppm
for index, mzs_i in enumerate(mzs):
    ppm = ppm_error(mzs_i, massb, massy)
    if df.loc[index, "relative_abundance"]<5:
            continue
    if ppm[0] < 50:
        df.loc[index, "ppm"] = ppm[0]
        df.loc[index, "ppm_match"] = ppm[1]
        # prefer y than b
    ## 补丁ppm = NONE

# add color
addcolor(df, "ppm_match", "color")
addcolor(df, "gaussian_sim_match", "color_gs")

pd.set_option("display.max_rows", 200)
print(df)


fig, ((ax1, ax3), (ax2, ax4)) = plt.subplots(
    nrows=2,
    ncols=2,
    figsize=(16, 9.6),# (8,4.8)
    sharex=True)

fig.suptitle(
    f"Annotated MS/MS Spectrum of\n{num} {aaseq})",
    fontsize=16,
    fontweight="bold",
    y=0.98
)

draw_plot(ax1, 
          df,  
          "relative_abundance", 
          "color",
          "ppm_match", 
          "Relative Abundance (%)", 
          True, 
          title="PPM")
draw_plot(ax2, 
          df,  
          "ints", 
          "color",
          "ppm_match", 
          "Ints", 
          False)
draw_plot(ax3, 
          df,  
          "relative_abundance", 
          "color_gs",
          "gaussian_sim_match", 
          "Relative Abundance (%)", 
          True,
          title="Gaussian Similarity")
draw_plot(ax4, 
          df,  
          "ints", 
          "color_gs",
          "gaussian_sim_match", 
          "Ints", 
          False)



plt.tight_layout()
fig_name = f"/home/student/sf_sharedFolders/project_v3/msms2_{num}_{aaseq[:8]}.png"
plt.savefig(fig_name, dpi = 1000)

# html
# sequence mutiplt
# modules
# work flow
# reason choose this msms
# pylab create msms picture
# aa mass in it own table remove from main function
# amino acid masses
#  data
# tricks and bugs