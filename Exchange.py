#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 16:07:07 2022

@author: caren
"""
import simpy
from Output import ExchangeStatistics

# Abstract version of the Exchange
class Exchange(object):
    def __init__(self, env, name, settings):
        self.env = env
        self.name = name
        self.s = settings
        self.capacity = simpy.Resource(env,1) # Only one person can trade against the LOB at a time
        self.marketSellTransactions = Queue(env, self.s)
        self.marketBuyTransactions = Queue(env, self.s)
        self.allPrices = []
        self.statistics = ExchangeStatistics()
    
    def getBuySellQuantity(self, time):
        return (self.marketBuyTransactions.getLiqSec()*time, self.marketSellTransactions.getLiqSec()*time)
    
    def addStatistics(self, sellVol = 0, buyVol = 0):
        (UnitBuy, UnitSell) = self.__quoteSize__()
        self.statistics.add(price = self.spot(), \
                            time = self.env.now, \
                            sellVol = sellVol, \
                            buyVol = buyVol, \
                            spread = self.bestAskPrice() - self.bestBidPrice(), \
                            bigSpread = self.__bigSpread__(),\
                            informedProb = self.s.infP,\
                            avSell = self.getLiquidityDown(),\
                            avBuy = self.getLiquidityUp(), \
                            unitBuy = UnitBuy, \
                            unitSell = UnitSell)

class Queue():
    def __init__(self, env, s):
        # self.liqMax = 8+(1-self.liqP)*12
        # self.initTimePerVol = ((self.liqMin+self.liqMax)/2) / (0.5 * 0.5 * (self.transSize[0]+self.transSize[1])/2)
        self.env = env
        self.queue = []
        avgTimeStep = (s.liqMin+s.liqMax)/2/4 #Divide by 4, since it can be sell/buy and market/liquidity order.
        avgSizeOrder = (s.transSize[0]+s.transSize[1])/2
        initTime = -(s.lookNTransactionsBack-1)*avgTimeStep
        for i in range(s.lookNTransactionsBack):
            self.queue.append((initTime+i*avgTimeStep, avgSizeOrder))
    
    def getLiqSec(self):
        totVol = sum([x[1] for x in self.queue])
        time = self.env.now-self.queue[0][0]
        return totVol/time
    
    def push(self, vol):
        self.queue.pop(0)
        self.queue.append((self.env.now, vol))
        
        

