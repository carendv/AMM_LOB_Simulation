#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 14:38:49 2022

@author: caren
"""
import random
import numpy as np
from matplotlib import pyplot as plt
from Trader import updateBuySellQuantity, retrieveFunds, probOrder

random.seed(10)
    
def trade(env, name, lob, s, o, r):    
    #print(name)
    now = env.now
    forced = False
    
    ###################
    # Liquidity check #
    ################### 
    lst_condition_result = [(lob.getLiquidityDown() < s.minLiquidity, True), (lob.getLiquidityUp() < s.minLiquidity, False)]
    random.shuffle(lst_condition_result)
    
    for condition, param in lst_condition_result:
        if condition:
            forced = True
            buy = param
    
    ##################
    # Initialization #
    ##################
    spot = lob.spot()
    bid = lob.bestBidPrice()
    ask = lob.bestAskPrice()
    informed = r.random() < s.infP
    belP = round(s.trueP+(min(max(np.random.normal(0, 1), -3), 3))) if informed and not forced else spot
    buy = 1/(1+np.e**(-s.agressiveness*(belP-spot))) > r.random() if not forced else buy
    kind = "L" if forced else None
    
    ######################
    # Decision variables #
    ######################
    # If not forced, see whether a limit order or market order will be placed
    if not forced:
        if informed and ((buy and belP < ask) or (not buy and belP > bid)):
            kind = "L"
        elif informed:
            kind = "M"
        elif not informed:
            if buy:
                kind = "L" if probOrder(belP, ask) > random.random() else "M"
            else:
                kind = "L" if probOrder(bid, belP) > random.random() else "M"
        
    (time, amount) = s.getBuy(kind) if buy else s.getSell(kind)
    
    (buyQuan, sellQuan) = lob.getBuySellQuantity(time)
    (buyQuan, sellQuan) = updateBuySellQuantity(buyQuan, sellQuan, belP, spot)
    expP = lob.expP(sellQuan-buyQuan) # For statistics only
    
    ##################
    # Order decision #
    ##################
    if time > 0:
        if buy and belP < ask:
            shift = 3 if informed else 6
            undercut = 1/(1+np.e**(-1*(belP-bid-shift))) > r.random()
            if undercut and belP >= bid+1 and bid+1 <= ask-1:
                p = bid+1
            else:
                p = min(min(belP, bid), ask-1)
                p2 = p-1
                excess = sellQuan-buyQuan
                liqAt = lob.getLiquidity(p2) - amount
                if excess >= liqAt or p2 >= bid:
                    p = p2
        elif not buy and belP > bid:
            shift = 3 if informed else 6
            undercut = 1/(1+np.e**(-1*(ask-belP-shift))) > r.random()
            if undercut and belP <= ask-1 and ask-1 >= bid+1:
                p = ask-1
            else:
                p = max(max(belP, ask), bid+1)
                p2 = p+1
                excess = buyQuan - sellQuan
                liqAt = lob.getLiquidity(p2) + amount
                if excess >= liqAt or p2 <= ask:
                    p = p2
        # A not informed trader only does liquidity order when some hope on execution
        if not forced:
            if buy and p == bid:
                liqAt = lob.getLiquidity(p) - amount
                if sellQuan < liqAt:
                    kind = "M"
            elif not buy and p == ask:
                liqAt = lob.getLiquidity(p) + amount
                if buyQuan < liqAt:
                    kind = "M"
            
    liq = kind == "L"
        
    ###################
    # Order execution #
    ###################  
    if liq:
        with lob.capacity.request() as request:
            yield request
            p = round(p)
            kind = "X" if buy else "M"
            trading = amount*belP/p if buy else amount
            (X, M, contract) = lob.add(trading, p, name)
        yield env.timeout(time) | contract.succes
    with lob.capacity.request() as request:
        yield request
        if liq:
            (X, M, completionPer, filled) = retrieveFunds(amount, contract, name, belP, X, M)
            o.orders.append((name, spot, lob.spot(), expP, belP, amount, time, env.now-now, X, M, completionPer, "L"))
        elif buy:
            (X, M) = lob.tradeM(-amount*belP)
            completionPer = (-amount*belP)/X/belP
            o.orders.append((name, spot, lob.spot(), expP, belP, amount, time, env.now-now, X, M, completionPer, "M"))
        else:
            (X, M) = lob.trade(amount)
            completionPer = M/amount/belP
            o.orders.append((name, spot, lob.spot(), expP, belP, amount, time, env.now-now, X, M, completionPer, "M"))
    