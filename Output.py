#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 16:12:25 2022

@author: caren
"""
from statsmodels.tsa.stattools import adfuller

class Output(object):
    def __init__(self):
        self.exchange = None
        self.numTrades = 0
        self.points = []
        self.prices = []
        self.orders = []
        self.spread = []
    
    def randomWalk(self):
        results = adfuller(self.prices)
        print(f"ADF Statistic: {results[0]}")
        print(f"p-value: {results[1]}")
        print("Critical Values:")
        for key, value in results[4].items():
            print("\t%s: %.3f" % (key, value))