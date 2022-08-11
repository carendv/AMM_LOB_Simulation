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
        self.bigSpread = []
        self.informedProb = []
        self.availableBuy = []
        self.availableSell = []
        self.unitBuy = []
        self.unitSell = []
    
    def add(self, price, time, sellVol, buyVol, spread, bigSpread, informedProb, avBuy, avSell, unitBuy, unitSell):
        self.prices.append(price)
        self.times.append(time)
        self.sellVol.append(self.sellVol[-1]+sellVol)
        self.buyVol.append(self.buyVol[-1]+buyVol)
        self.spread.append(spread)
        self.bigSpread.append(bigSpread)
        self.informedProb.append(informedProb)
        self.availableBuy.append(avBuy)
        self.availableSell.append(avSell)
        self.unitBuy.append(unitBuy)
        self.unitSell.append(unitSell)
        
        