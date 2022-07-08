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
        elif incr + 10 > env.now:
            s.infP = s.maxInfP
            decr = env.now + s.shockDecrTime
            start = env.now
            incr = 0
        elif decr > env.now:  
            s.infP = s.maxInfP - (s.maxInfP-s.minInfP)/(decr-start)*(env.now-start)
        elif decr + 10 > env.now:
            s.infP = s.minInfP
            decr = 0

        env.process(trade(env, trades, o.exchange, s, o, random))
        trades += 1

def visualizeResults(results):
    statistics = results.exchange.statistics
    statistics.buyVol.remove(0)
    statistics.sellVol.remove(0)
    
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
    plotWithInformed(statistics, points)
    # plotWithInformed(statistics.times, statistics.prices, points, "Price")
    # plotWithInformed(statistics.times, statistics.spread, points, "Unit spread")
    # plotWithInformed(statistics.times, statistics.sellVol, points, "Cumulative sell volume")
    # plotWithInformed(statistics.times, statistics.buyVol, points, "Cumulative buy volume")

def plotWithInformed(data, points):
    names = ["Price", "Unit spread", "Cumulative sell volume", "cumulative buy volume"]
    dataS = [data.prices, data.spread, data.sellVol, data.buyVol]
    num = len(names)
    rows = 2
    columns = num/rows
    
    plt.figure()
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
    for i in range(NAMM):
        settings = Settings(NAMM=1, NLOB=0, shocks=shocks, days=days, seed=seed)
        output = Output(settings)
        random.seed(settings.seed)
        env = simpy.Environment()
        env.process(setup(env, settings, output, "AMM", i))
        env.run(until=settings.totTime)
        visualizeResults(output)
        outputs.append(output)
        
    for i in range(NLOB):
        settings = Settings(NAMM=0, NLOB=1, shocks=shocks, days=days, seed=seed)
        output = Output(settings)
        random.seed(settings.seed)
        env = simpy.Environment()
        env.process(setup(env, settings, output, "LOB", i))
        env.run(until=settings.totTime)
        visualizeResults(output)
        outputs.append(output)

    return outputs

results = simulation(NAMM=1, NLOB=1, shocks=1, days=2, seed=100)

