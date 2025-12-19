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

def cal_theory_masses(b,massb,d3,x=1):
    for i,seq in enumerate(b):
        for j in seq:
            if j in d3:
                massb[i] = float(massb[i]) + float(d3[j])
        massb[i] = massb[i] + x
    return massb
                
def ppm_error(oberve, mwb, mwy):
    ppm = []
    for theory in mwb+mwy:
        ppm.append(((abs(theory - oberve) / theory * 1e6), theory))
    return min(ppm)

def gaussian_similarity(observe, theory, sigma=0.4):
    sim = np.exp(-(theory - observe)**2 / (2 * sigma**2))
    return sim

def plotMsMs(ax, df, colors,percent):
    df_color = df[df["color"] == colors] 
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