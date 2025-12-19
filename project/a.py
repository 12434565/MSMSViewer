import sys
from base64 import b64decode
from array import array
from model import msms
import pandas as pd
import numpy as np
from model1 import ppm_error, gaussian_similarity, cal_theory_masses, getd3, get_peaks_from_xml

import matplotlib.pyplot as plt

sigma = 0.4
filename = sys.argv[1] # 17mix_test2.mzxml
num = int(sys.argv[2])  # 1298
aaseq = sys.argv[3]  # TYDSYLGDDYVR

peakselt = msms()
peakselt.text = get_peaks_from_xml("17mix_test2.mzxml", num).text
l = len(aaseq)
b = [None]*(l-1)
y = [None]*(l-1)

for i in range(l-1):
    b[i] = aaseq[:i+1]
    y[i] = aaseq[i+1:]

d3 = getd3()
    
massb = [0]*(l-1)
massy = [0]*(l-1)
massb = cal_theory_masses(b,massb,d3,1)
massy = cal_theory_masses(y,massy,d3,19)


d2 = {}
d2["b"] = b
d2["b_index"] = (f"b{m}" for m in range(1, len(b)+1))
d2["mwb"] = massb
d2["y"] = y
d2["y_index"] = (f"y{m}" for m in reversed(range(1, len(b)+1)))
d2["mwy"] = massy
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
d["gaussian_sim_match"] = [""]*len(mzs)
d["ppm"] = [""]*len(mzs)
d["ppm_match"] = [""]*len(mzs)
d["color"] = [""]*len(mzs)

df = pd.DataFrame(data=d)


# calculate gaussian similarity
for mw_i in massb+massy:
    for mzs_i in mzs:
        sim = gaussian_similarity(mw_i, mzs_i, sigma=0.4)
        if df.loc[df["mzs"] == mzs_i, "relative_abundance"].values[0]<5:
            continue
        if sim > 0.9:
            # print(f"Found matching b-ion: {mwb_i} ~ {mzs_i}")
            index2 = df2.loc[df2["mwb"] == mw_i, "b_index"]
            if index2.empty:
                index2 = df2.loc[df2["mwy"] == mw_i, "y_index"]
            df.loc[df["mzs"] == mzs_i, "gaussian_sim_match"] = index2.values[0]

# calculate ppm
for mzs_i in mzs:
    ppm = ppm_error(mzs_i, massb, massy)
    if ppm[0] < 100:
        df.loc[df["mzs"] == mzs_i, "ppm"] = ppm[0]
        index2 = df2.loc[df2["mwb"] == ppm[1], "b_index"]
        if index2.empty:
            index2 = df2.loc[df2["mwy"] == ppm[1], "y_index"]
        if df.loc[df["mzs"] == mzs_i, "relative_abundance"].values[0]<5:
            continue
        df.loc[df["mzs"] == mzs_i, "ppm_match"] = index2.values[0]
    
        # prefer b over y

# add color
for i in df.ppm_match.index:
    if df.loc[i, "ppm_match"] != "":
        if df.loc[i, "ppm_match"].startswith("b"):
            df.loc[i, "color"] = "blue"
        else:
            df.loc[i, "color"] = "red"
    else:
        df.loc[i, "color"] = "black"
pd.set_option("display.max_rows", 200)
print(df)


# plt.style.use('_mpl-gallery')
# make data

# # plot
# fig, ax = plt.subplots()
# markerline, stemlines, baseline = ax.stem(
#     x, y,
#     markerfmt=' ',   
#     basefmt=' ')      
# ax.set(ylim = (0,100),xlabel='m/z', ylabel='Relative Abundance (%)')
# stemlines.set_linewidth(1)
# plt.tight_layout()
# # plt.show()
# plt.savefig("msms.png", dpi = 1000)


x = df.mzs
y = df.relative_abundance
fig, ax = plt.subplots()
ax.vlines(x, 0, y,
       colors= df.color, 
       linestyles='solid',
       label='',
       lw = 1)
ax.set(ylim = (0,100),xlabel='m/z', ylabel='Relative Abundance (%)')
for _, row in df.dropna(subset=["ppm_match"]).iterrows():
    ax.text(
        row["mzs"],                      
        row["relative_abundance"] + 1,   
        row["ppm_match"],                
        ha="center",
        va="bottom",
        fontsize=8,
        color=row["color"]
    )
ax.set_title(
    f"Annotated MS/MS Spectrum of \n{num} {aaseq}",
    fontsize=14,
    fontweight="bold",
    pad=10
)
plt.tight_layout()
fig_name = f"/home/student/sf_sharedFolders/project/msms2_{num}_{aaseq[:8]}.png"
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