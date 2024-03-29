#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 16:10:26 2022

@author: caren
"""
import math
import numpy as np

minPriceRange = 950

class Settings(object):
    def __init__(self, NAMM=1, NLOB=0, shocks=1, days=3, seed=100):
        # Variables that can be set
        self.NAMM = NAMM
        self.NLOB = NLOB
        self.shocks = shocks
        
        self.seed = seed
        self.r = np.random
        self.r.seed(100)
                
        # Variables with respect to transaction
        self.orP = 1000
        self.trueP = self.orP
        self.fee = 0.0005
        self.transSize = [1000 , 5000]
        
        # Variables w.r.t. liquidity of market
        self.initialBidAsk = 2
        self.agressiveness = 0.5
        self.minLiquidity = (self.transSize[0]+self.transSize[1]) /2 * 10
        self.liqP = 1
        self.liqMin = int(4+(1-self.liqP)*6)
        self.liqMax = int(8+(1-self.liqP)*12)
        self.lookNTransactionsBack = 100
        self.AMMmax = self.trueP+2 # The maximal price in AMM
        self.minPriceRange = minPriceRange
        self.maxPriceRange = 1100
        
        # Variables w.r.t. duration of simulation
        self.days = days
        self.totTime = self.days*60*60*8
        
        # Variables w.r.t. informed traders
        self.infP = 0.1
        self.minInfP = self.infP
        self.maxInfP = 0.75
        
        # Variables w.r.t. shocks, when they happen and how long
        self.shockStep = 10
        self.shockTime = self.totTime/self.shocks if not self.shocks==0 else math.inf
        self.shockWait = self.shockTime*0.3
        self.shockIncrTime = self.shockTime*0.1
        self.shockDecrTime = self.shockTime*0.1
        self.shockAfterTime = self.shockTime-self.shockWait-self.shockIncrTime-self.shockDecrTime
        
    def getAmountTime(self, timeFunc, buy, kind = False,):
        time = round(max(timeFunc(),0))
        amount = self.r.uniform(self.transSize[0],self.transSize[1])
        
        if kind == "M":
            time = 0
        
        while kind == "L" and time == 0:
            time = round(max(timeFunc(),0))
        
        if time > 0:
            amount = amount*1.3662
        
        amount = round(amount*(1-2*buy))
    
        return (time, amount)
    
    def getBuy(self, kind=False):
        return self.getAmountTime(lambda : self.r.normal(46.92*60,72.31*60), True, kind)
    
    def getSell(self, kind=False):
        return self.getAmountTime(lambda : self.r.normal(34.15*60,53.94*60), False, kind)


        
        


