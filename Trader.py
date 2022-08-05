#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 16:07:47 2022

@author: caren
"""
import numpy as np
    
def retrieveFunds(amount, contract, name, belP, assets=0, money=0):
    if amount < 0:
        (assets2, money, filled) = contract.retrieve(money, name)
        assets += assets2
        completionPer = ((-amount*belP)/assets)/belP
    else:
        (assets, money2, filled) = contract.retrieve(assets, name)
        money += money2
        completionPer = money/amount/belP
    return (assets, money, completionPer, filled)

def probOrder(belP, marketPrice):
    return 1/(1+np.e**(-(marketPrice-belP)/1))

def updateBuySellQuantity(buyQuan, sellQuan, belP, spot):
    totQuan = buyQuan + sellQuan
    p = belP/spot
    weight = (1+abs(1-p))/2
    buyQuan = buyQuan*(1-weight) + totQuan*p/(1 + p)*weight
    sellQuan = totQuan - buyQuan
    return (buyQuan, sellQuan)
    