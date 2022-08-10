#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 16:07:57 2022

@author: caren
"""

import simpy
import random
import math
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from statsmodels.tsa.stattools import adfuller
import statsmodels.api as sm

from Settings import Settings
from Output import Output
from AMM import AMM
from LOB import LOB
from lobTrader import trade as lobTrade
from ammTrader import trade as ammTrade
import warnings

warnings.filterwarnings("ignore")

# Setup simulation, given the environment and settings
def setup(env, s, o, kindExchange, i):
    random.seed(101)
    
    if kindExchange == "AMM":
        o.exchange = AMM(env, 'AMM_' + str(i), s)
        trade = ammTrade
    elif kindExchange == "LOB":
        o.exchange = LOB(env, 'LOB_' + str(i), s)
        trade = lobTrade

    shockIn = 1
    incr = math.inf
    decr = math.inf
    start = -1
    trades = 0
        
    while True:
        yield env.timeout(random.randint(s.liqMin,s.liqMax))
        # Call process for a new client
        if (shockIn-1)*s.shockTime + s.shockWait < env.now:
            s.trueP += s.shockStep
            shockIn += 1
            incr = env.now + s.shockIncrTime
            start = env.now
        elif incr > env.now:
            s.infP = (s.maxInfP-s.minInfP)/(incr-start)*(env.now-start)+s.minInfP
        elif incr <= env.now:
            s.infP = s.maxInfP
            decr = env.now + s.shockDecrTime
            start = env.now
            incr = math.inf
        elif decr > env.now:  
            s.infP = s.maxInfP - (s.maxInfP-s.minInfP)/(decr-start)*(env.now-start)
        elif decr <= env.now:
            s.infP = s.minInfP
            decr = math.inf

        env.process(trade(env, trades, o.exchange, s, o, random))
        trades += 1

def visualizeResults(results):
    statistics = results.exchange.statistics
    statistics.buyVol.remove(0)
    statistics.sellVol.remove(0)
    
    print("#############")
    print(f"### {results.exchange.name} ###")
    print("#############")
    
    ##########################################################
    ### Find plots where number of informed trades changes ###
    ##########################################################
    points = []
    state = "c"
    infP = statistics.informedProb
    for (i, x) in enumerate(infP):
        if i == 0:
            continue
        elif infP[i]-infP[i-1] > 0 and state == "c":
            points.append(i)
            state = "i"
        elif infP[i]-infP[i-1] < 0 and state == "i":
            points.append(i)
            state = "d"
        elif infP[i] == results.exchange.s.minInfP and state == "d":
            points.append(i-1)
            state = "c"  
    
    ####################
    ### Create plots ###
    ####################
    plotWithInformed(statistics, points, results.exchange.name)
    
    #######################
    ### Price discovery ###
    #######################
    for i in range(round(len(points)/3)):
        a = results.exchange.s.orP + i*results.exchange.s.shockStep
        b = results.exchange.s.orP + (i+1)*results.exchange.s.shockStep
        index = next(j for j,v in enumerate(statistics.prices) if v >= b and j >= points[i*3])
        buyVol = round(statistics.buyVol[index]-statistics.buyVol[points[i*3]])
        sellVol = round(statistics.sellVol[index]-statistics.sellVol[points[i*3]])
        time = statistics.times[index]-statistics.times[points[i*3]]
        trades = index-points[i*3]
        print(f"Shock {i+1} took the price from {a} to {b}.")
        print(f"{time} seconds went by and there where {trades} trades.")
        print(f"The buy volume was {buyVol}, while the sell volume was {sellVol}.")
        print()
    
    ####################
    ### Random walks ###
    ####################  
    # There should be a random walk when the informed trades are on the lower threshold
    points.insert(0, 0)
    points.append(len(statistics.prices)-1)
    pvaluesAF = []
    pvaluesK = []
    for i in range(int(np.ceil(len(points)/3))):
        statPart = statistics.prices[points[i*3]:points[i*3+1]]
        
        # Do an adfuller test
        af = adfuller(statPart)
        pvaluesAF.append(round(af[1], 4))
        
        # Do an KPSS test
        k = sm.tsa.stattools.kpss(statPart, regression='ct')
        pvaluesK.append(round(k[1], 4))
    
    print()
    print(f"Found p-values of adfuller test stationary parts: {pvaluesAF}")
    print(f"Found p-values of kpss test stationary parts: {pvaluesK}")
    if all(i > 0.05 for i in pvaluesAF) and all(i < 0.05 for i in pvaluesK):
        print("The market has perfect information.")
    else:
        print("The market prices are influenced by the past.")
    print()
    
        
    #############################
    ### Completion Percentage ###
    ############################# 
    #(name, spot, belP, amount, time, guessedTime, env.now-now, assets, money, completionPer, filled)
    orders = pd.DataFrame(results.orders, columns=["Name", "Spot", "SpotA", "expP", "belP", "Amount", \
                                        "Time", "ActualTime", "Assets", \
                                        "Money", "CompletionPer", "Filled"])
    plt.figure()
    plt.hist(orders.CompletionPer, bins=50)
    plt.title(results.exchange.name)
    plt.axvline(orders.CompletionPer.mean(), color='k', linestyle='dashed', linewidth=1)
    locs, _ = plt.yticks() 
    plt.yticks(locs,np.round(locs/len(orders.CompletionPer),3))
    plt.xlabel("Completion percentage")
    plt.ylabel("Frequency")
    plt.show()
    
    print(f"Completion percentage mean: {orders.CompletionPer.mean()}")
    print(f"Completion percentage liquidity providers: {orders[orders.Filled == 'L'].CompletionPer.mean()}")
    print(f"Completion percentage market providers: {orders[orders.Filled == 'M'].CompletionPer.mean()}")
    
    return points
    
    
def plotWithInformed(data, points, name):
    names = ["Price", "Unit spread", "1000 spread", "Available to buy"]
    bigSpread = [min(i, 10) for i in data.bigSpread]
    dataS = [data.prices, data.spread, bigSpread, data.availableBuy]
    num = len(names)
    rows = 2
    columns = num/rows
    
    fig = plt.figure()
    fig.suptitle(name)
    for i in range(num):
        plt.subplot(rows, int(np.ceil(columns)), i+1)
        plt.plot(data.times, dataS[i])
        plt.xlabel("Time (s)")
        plt.ylabel(names[i])
        plt.grid(True)    
        for i in range(round(len(points)/3)):
            plt.axvspan(data.times[points[i*3]], data.times[points[i*3+1]], alpha = 0.5, color='g')
            plt.axvspan(data.times[points[i*3+1]], data.times[points[i*3+2]], alpha = 0.5, color='r')
    plt.tight_layout()  
    
def simulation(NAMM = 1, NLOB = 1, shocks=0, days=3, seed=100): 
    outputs = []
    points = []
    for i in range(NAMM):
        settings = Settings(NAMM=1, NLOB=0, shocks=shocks, days=days, seed=seed)
        output = Output(settings)
        random.seed(settings.seed)
        env = simpy.Environment()
        env.process(setup(env, settings, output, "AMM", i))
        env.run(until=settings.totTime)
        point = visualizeResults(output)
        outputs.append(output)
        points.append(point)
        
    for i in range(NLOB):
        settings = Settings(NAMM=0, NLOB=1, shocks=shocks, days=days, seed=seed)
        output = Output(settings)
        random.seed(settings.seed)
        env = simpy.Environment()
        env.process(setup(env, settings, output, "LOB", i))
        env.run(until=settings.totTime)
        point = visualizeResults(output)
        outputs.append(output)
        points.append(point)

    return outputs, points

#results, points = simulation(NAMM=1, NLOB=1, shocks=1, days=5, seed=100)

