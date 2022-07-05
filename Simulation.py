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

from Settings import Settings
from Output import Output
from AMM import AMM
from LOB import LOB
from Trader import trade

# Setup simulation, given the environment and settings
def setup(env, s, o, kindExchange, i):
    random.seed(101)
    
    if kindExchange == "AMM":
        o.exchange = AMM(env, 'AMM_' + str(i), s)
    elif kindExchange == "LOB":
        o.exchange = LOB(env, 'LOB_' + str(i), s)

    shockIn = 1
    incr = -math.inf
    decr = -math.inf
    start = -1
        
    while True:
        yield env.timeout(random.randint(s.liqMin,s.liqMax))
        # Call process for a new client
        if (shockIn-1)*s.shockTime + s.shockWait < env.now:
            s.trueP += s.shockStep
            shockIn += 1
            incr = env.now + s.shockIncrTime
            start = env.now
            o.points.append(o.numTrades)
        elif incr > env.now:
            s.infP = (s.maxInfP-s.minInfP)/(incr-start)*(env.now-start)+s.minInfP
        elif incr + 10 > env.now:
            s.infP = s.maxInfP
            o.points.append(o.numTrades)
            decr = env.now + s.shockDecrTime
            start = env.now
            incr = 0
        elif decr > env.now:  
            s.infP = s.maxInfP - (s.maxInfP-s.minInfP)/(decr-start)*(env.now-start)
        elif decr + 10 > env.now:
            s.infP = s.minInfP
            o.points.append(o.numTrades)
            decr = 0

        o.times.append(env.now)
        env.process(trade(env, o.numTrades, o.exchange, s, o, random))
        o.numTrades += 1

def visualizeLOBResults(results):
    orders = pd.DataFrame(results.orders, columns=['Name', 'Spot', 'belP', 'AmountWanted', 'Time', 'GuessedTime', 'PassedTime', 'Assets', 'Money', 'CompletionPer', 'Kind'])
    
    plt.figure()
    plt.plot(results.prices)    
    for i in range(round(len(results.points)/3)):
        plt.axvspan(results.points[i*3], results.points[i*3+1], alpha = 0.5, color='g')
        plt.axvspan(results.points[i*3+1], results.points[i*3+2], alpha = 0.5, color='r')
    plt.show()
    
    plt.figure()
    plt.plot([x for (_, x, _) in results.exchange.allPrices]) #Buy
    plt.figure()
    plt.plot([x for (_, _, x) in results.exchange.allPrices]) # Sell
    
    # Bid-Ask spread 
    plt.figure()
    plt.plot(results.spread)

    mo = orders[orders.Kind == "M"]
    moBuy = mo[mo.AmountWanted < 0]
    moSell = mo[mo.AmountWanted > 0]
    lo = orders[orders.Kind != "M"]
    loBF = lo[(lo.AmountWanted < 0) & (lo.Kind)]
    loBU = lo[(lo.AmountWanted < 0) & (lo.Kind == False)]
    loSF = lo[(lo.AmountWanted > 0) & (lo.Kind)]
    loSU = lo[(lo.AmountWanted > 0) & (lo.Kind == False)]
    
    print("Market buy orders:", len(moBuy))
    print("Market sell orders:", len(moSell))
    print("Filled limit buy orders:", len(loBF))
    print("Unfilled limit buy orders:", len(loBU))
    print("Filled limit sell orders:", len(loSF))
    print("Unfilled limit sell orders:", len(loSU))
    
    lo = orders[orders.Time>0]
    plt.figure()
    plt.hist(lo.GuessedTime-lo.PassedTime, bins = np.arange(-5000, 5000, 10))

def visualizeAMMResults(results):
    orders = pd.DataFrame(results.orders, columns=['Name', 'Spot', 'belP', 'AmountWanted', 'Time', 'GuessedTime', 'PassedTime', 'Assets', 'Money', 'CompletionPer', 'Kind'])
    plt.figure()
    
    plt.plot(results.prices)    
    for i in range(round(len(results.points)/3)):
        plt.axvspan(results.points[i*3], results.points[i*3+1], alpha = 0.5, color='g')
        plt.axvspan(results.points[i*3+1], results.points[i*3+2], alpha = 0.5, color='r')
    
    lo = orders[orders.Time>0]
    plt.figure()
    plt.hist(lo.GuessedTime-lo.PassedTime, bins = np.arange(-5000, 5000, 10))
        
    # Bid-Ask spread 
    plt.figure()
    plt.plot([x if x < 10 else None for x in results.spread])
    # TODO: Look at big values
    
def simulation(NAMM = 1, NLOB = 1, shocks=0, days=3, seed=100): 
    outputs = []
    for i in range(NAMM):
        settings = Settings(NAMM=1, NLOB=0, shocks=shocks, days=days, seed=seed)
        output = Output(settings)
        random.seed(settings.seed)
        env = simpy.Environment()
        env.process(setup(env, settings, output, "AMM", i))
        env.run(until=settings.totTime)
        visualizeAMMResults(output)
        outputs.append(output)
        
    for i in range(NLOB):
        settings = Settings(NAMM=0, NLOB=1, shocks=shocks, days=days, seed=seed)
        output = Output(settings)
        random.seed(settings.seed)
        env = simpy.Environment()
        env.process(setup(env, settings, output, "LOB", i))
        env.run(until=settings.totTime)
        visualizeLOBResults(output)
        outputs.append(output)

    return outputs

results = simulation(NAMM=1, NLOB=1, shocks=1, days=2, seed=100)
results[0].timeDiscovery()

