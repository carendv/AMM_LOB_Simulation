#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 16:12:25 2022

@author: caren
"""

class Output(object):
    def __init__(self, settings):
        self.exchange = None
        self.settings = settings
        self.orders = []

class ExchangeStatistics(object):
    def __init__(self):
        self.prices = []
        self.times = []
        self.sellVol = [0]
        self.buyVol = [0]
        self.spread = []
        self.informedProb = []
    
    def add(self, price = None, time = None, sellVol = 0, buyVol=0, spread = None, informedProb = None):
        self.prices.append(price)
        self.times.append(time)
        self.sellVol.append(self.sellVol[-1]+sellVol)
        self.buyVol.append(self.buyVol[-1]+buyVol)
        self.spread.append(spread)
        self.informedProb.append(informedProb)
        
        