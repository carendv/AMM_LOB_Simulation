#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 13:34:42 2022

@author: caren
"""
from Simulation import simulation 
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

results, points = simulation(NAMM=1, NLOB=1, shocks=1, days=3, seed=100)

names = ["Price", "Unit spread", "1000 spread", "Available to buy", "Available to sell", "quote Buy", "quote Sell"]
data1 = results[0].exchange.statistics
data2 = results[1].exchange.statistics
bigSpread1 = [min(i, 10) for i in data1.bigSpread]
bigSpread2 = [min(i, 10) for i in data2.bigSpread]
dataS1 = [data1.prices, data1.spread, bigSpread1, data1.availableBuy, data1.availableSell, data1.unitBuy, data1.unitSell]
dataS2 = [data2.prices, data2.spread, bigSpread2, data2.availableBuy, data2.availableSell, data2.unitBuy, data2.unitSell]
name1 = "AMM"
name2 = "LOB"

num = len(names)
rows = 2
columns = num/rows

#fig = plt.figure()
for i in range(num):
#    plt.subplot(rows, int(np.ceil(columns)), i+1)
    plt.figure()
    plt.plot(data1.times, dataS1[i], alpha = 0.8, label=name1)
    plt.plot(data2.times, dataS2[i], alpha = 0.8, label=name2)
    plt.xlabel("Time (s)")
    plt.ylabel(names[i])
    plt.legend()
    plt.grid(True)    
    for i in range(round(len(points[0])/3)-1):
        plt.axvspan(data1.times[points[0][i*3+1]], data1.times[points[0][i*3+2]], alpha = 0.5, color='g')
        plt.axvspan(data1.times[points[0][i*3+2]], data1.times[points[0][i*3+3]], alpha = 0.5, color='r')
#plt.tight_layout() 

fig = plt.figure()
orders1 = pd.DataFrame(results[0].orders, columns=["Name", "Spot", "SpotA", "expP", "belP", "Amount", \
                                    "Time", "ActualTime", "Assets", \
                                    "Money", "CompletionPer", "Filled"])
orders2 = pd.DataFrame(results[1].orders, columns=["Name", "Spot", "SpotA", "expP", "belP", "Amount", \
                                    "Time", "ActualTime", "Assets", \
                                    "Money", "CompletionPer", "Filled"])
plt.figure()
plt.hist(orders1.CompletionPer, bins=50, label=name1, alpha = 0.5, color='b')
plt.axvline(orders1.CompletionPer.mean(), color='b', linestyle='dashed', linewidth=1)
plt.hist(orders2.CompletionPer, bins=50, label=name2, alpha = 0.5, color='g')
plt.axvline(orders2.CompletionPer.mean(), color='g', linestyle='dashed', linewidth=1)
locs, _ = plt.yticks() 
plt.yticks(locs,np.round(locs/len(orders1.CompletionPer),3))
plt.xlabel("Completion percentage")
plt.ylabel("Frequency")
plt.legend()
plt.show()