#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 16:12:25 2022

@author: caren
"""
# Output statistics (to be collected)
#exchanges = []


# 'Name', 'Start', 'Time', 'Exchange', 'Buy?', 'Asset', 'Money', 'waitCosts', 'Bid', 'Ask'
#orders = []
#liquidityOrders = []

#numTrades = [0]
#undercut = []
#inqueue = []

class Output(object):
    def __init__(self):
        self.exchange = None
        self.numTrades = 0
        self.points = []
        self.prices = []
        self.orders = []
        self.spread = []