import numpy as np
import xml.etree.ElementTree as ET
import urllib.request
from model import msms

def get_peaks_from_xml(filename,num):
    peak = msms()
    # msRun -> instrument/dataProcessing/scan -> peak.text
    document = ET.parse(filename)
    root = document.getroot()
    # print(root.tag)
    ns = "{http://sashimi.sourceforge.net/schema/}"
    for ele in root.findall(ns+'scan'):
        if int(ele.attrib.get("num", 0)) == num:
            peak.text = ele.find(ns+'peaks').text
    if peak.text == "":
        print("the input number is out of range")
    return peak
                        
def getd3():
    d3 = {"G":57.02146,"A":71.03711,
      "S":87.03203,"P":97.05276,
      "V":99.06841,"T":101.04768,
      "C":103.00919,"L":113.08406,
      "I":113.08406,"N":114.04293,
      "D":115.02694,"Q":128.05858,
      "K":128.09496,"E":129.04259,
      "M":131.04049,"H":137.05891,
      "F":147.06841,"R":156.10111,
      "Y":163.06333,"W":186.07931}   
    return d3

def cal_theory_masses(b,d3,x=1):
    mass = []
    for i,seq in enumerate(b):
        mz = 0
        for j in seq:
            if j in d3:
                mz = float(mz) + float(d3[j])
        mz = mz + x
        if x == 1:
            mass.append((mz, f"b{i+1}"))
        elif x == 19:
            mass.append((mz,f"y{len(b)-i}"))
    return mass
                
def ppm_error(observe, massb, massy, bound):
    ppm_y = []
    ppm_b = []
    for theory, label_y in massy:
        ppm_y.append(((abs(theory - observe) / theory * 1e6), label_y))
    for theory, label_b in massb:
        ppm_b.append(((abs(theory - observe) / theory * 1e6), label_b))
    best_y = min(ppm_y, key=lambda x: x[0])
    best_b = min(ppm_b, key=lambda x: x[0])
    if (best_y[0] - bound) > best_b[0]:
        return best_b
    return best_y

def gaussian_similarity(observe, theory, sigma=0.4):
    sim = np.exp(-(theory - observe)**2 / (2 * sigma**2))
    return sim

def addcolor(df, methods, colorColumnName):
    for i in df.index:
        if df.loc[i, methods] == "":
            df.loc[i, colorColumnName] = "black"
        else:
            if df.loc[i, methods].startswith("b"):
                df.loc[i, colorColumnName] = "blue"
            else:
                df.loc[i, colorColumnName] = "red"
    return df

def plotMsMs(ax, df, colors,percent, colorColumnName):
    df_color = df[df[colorColumnName] == colors] 
    if colors == "black":
        zorders = 1
    else:
        zorders = 2
    if percent:
        ax.vlines(df_color.mzs, 
                0, 
                df_color.relative_abundance, 
                colors=colors, 
                zorder=zorders,
                linestyles='solid',
                lw = 1)
        return ax
    else:
        ax.vlines(df_color.mzs, 
                0, 
                df_color.ints, 
                colors=colors, 
                zorder=zorders,
                linestyles='solid',
                lw = 1)
        return ax
    
def draw_plot(ax, df,  
              ycolumn_name, # "relative_abundance" or "ints"
              colorColumnName, # "color" or "color_gs"
              method_score_column_name, # "ppm_match" or "gaussian_sim_match"
              ylabel_a, # "Relative Abundance (%)" or "init"
              percent, # True or False
              title=None,
              xlabel_a="m/z"):
    max = int(df[ycolumn_name].max()) + int(df[ycolumn_name].max()*0.05)
    plotMsMs(ax, df, "black", percent, colorColumnName)
    plotMsMs(ax, df, "blue", percent, colorColumnName)
    plotMsMs(ax, df, "red", percent, colorColumnName)
    ax.set(ylim = (0,max), ylabel= ylabel_a)
    for _, row in df.dropna(subset=[colorColumnName]).iterrows():
        ax.text(
            row["mzs"],                      
            row[ycolumn_name] + 1,   
            row[method_score_column_name],                
            ha="center",
            va="bottom",
            fontsize=8,
            color=row[colorColumnName]
        )
    if title != None:
        ax.set_title(
            f"{title}",
            fontsize=12,
            fontweight="bold",
            pad=10
        )
    return ax