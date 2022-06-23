#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 16:07:07 2022

@author: caren
"""
import simpy

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
    
    def getBuySellQuantity(self, time):
        return (self.marketBuyTransactions.getLiqSec()*time, self.marketSellTransactions.getLiqSec()*time)
    

class Queue():
    def __init__(self, env, s):
        self.env = env
        self.queue = []
        initTime = -(s.lookNTransactionsBack-1)*s.initTimePerVol
        for i in range(s.lookNTransactionsBack):
            self.queue.append((initTime+i*s.initTimePerVol, 1))
    
    def getLiqSec(self):
        totVol = sum([x[1] for x in self.queue])
        time = self.env.now-self.queue[0][0]
        return totVol/time
    
    def push(self, vol):
        self.queue.pop(0)
        self.queue.append((self.env.now, vol))
        
        

