#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 13:15:34 2022

@author: caren
"""
import random
import math
from AMM import AMM
from Trader import updateBuySellQuantity
import numpy as np

random.seed(10)
    
def trade(env, name, amm, s, o, r):
    print(name)
    now = env.now
    
    
    ###################
    # Liquidity check #
    ###################    
    lst_condition_result = [(amm.getLiquidityDown() < s.minLiquidity, True), (amm.getLiquidityUp() < s.minLiquidity, False)]
    random.shuffle(lst_condition_result)
    
    for condition, param in lst_condition_result:
        if condition:
            forced = True
            buy = param
    # TODO: Use this to force liquidity order
    
    ##################
    # Initialization #
    ##################
    spot = amm.spot()
    informed = r.random() < s.infP
    belP = round(s.trueP+(min(max(np.random.normal(0, 1), -3), 3))) if informed else spot
    buy = 1/(1+np.e**(-s.agressiveness*(belP-spot))) > r.random()
    (time, amount) = s.getBuy() if buy else s.getSell()
    
    ######################
    # Decision variables #
    ######################
    wait = time > 0
    liq = wait
    M = round(-min(amount, 0)*belP)
    X = max(amount, 0)
    
    ##################
    # Order decision #
    ##################
    (buyQuan, sellQuan) = amm.getBuySellQuantity(time)
    (buyQuan, sellQuan) = updateBuySellQuantity(buyQuan, sellQuan, belP, spot)
    expP = amm.getPriceX(sellQuan-buyQuan) # Price after waiting
    
    if buy and belP <= spot and time > 0:    
        expXn = amm.getX(-amount*belP) # Assets received when buying right now
        expXw = amm.getX(-amount*belP, expP) # Assets received after waiting 
        
        # Compute best range
        rL = max(min(math.ceil(expP), math.floor(spot)-2), math.floor(belP)-1)
        rU = max(amm.getLiquidityStart(rL), rL+2)
        
        # Compute expected fees
        if rU <= math.floor(spot):
            # Compute percentage of fees that occur in range
            sell = amm.getXPrice(rL)- amm.getXPrice(rU)
            perInRange = sell/(sellQuan - buyQuan)
            Lself = M/(math.sqrt(rU)-math.sqrt(rL))*(rU-rL)
            LOther = amm.getL(rL, rU)
            
            Xfee = amm.F/(1-amm.F)*sellQuan*perInRange*Lself/(Lself+LOther)
            Mfee = amm.F/(1-amm.F)*buyQuan*perInRange*math.sqrt(rL*rU)*Lself/(Lself+LOther)
            
            expXl = Xfee + M/math.sqrt(rL*rU) + amm.getX(Mfee, expP)
            
            if expXl <= expXn and expXn < expXw:
                liq = False
            elif expXl <= expXn and expXn >= expXw:
                liq = False
                wait = False
            elif expXl < expXw:
                liq = False
            # When rU > math.floor(spot) we know rL <= spot <= rU
            # Therefore, we are trying to get fees and so it is fine to do liquidity order.
            # TODO: Check this
    elif buy and time > 0: #(When belP > spot)
        # Trade immedeately. 
        # The price is expected to go up, thus range order would only sell X.
        # Thus buying not possible.
        # Futhermore, now is cheap.
        liq = False
        wait = False
    elif not buy and belP >= spot and time > 0:
        expMn = amm.getM(amount)
        expMw = amm.getM(amount, expP)
        
        expMn = amm.getM(amount)
        expMw = amm.getM(amount, expP)
        
        # Compute best range
        rU = min(max(math.floor(expP), math.ceil(spot)+2), math.ceil(belP)+1)
        rL = min(amm.getLiquidityStart(rU), rU-2)
        if rL >= math.ceil(spot):
            buyInRange = amm.getMPrice(rU)-amm.getMPrice(rL)
            perInRange = buyInRange/(buyQuan-sellQuan)
            Lself = math.sqrt(rL*rU)/(math.sqrt(rU)-math.sqrt(rL))*X*(rU-rL)
            LOther = amm.getL(rL, rU)
            
            Xfee = sellQuan*perInRange*Lself/(Lself+LOther)
            Mfee = buyQuan*perInRange*math.sqrt(rL*expP)*Lself/(Lself+LOther)
            
            expMl = Mfee + X*math.sqrt(rL*rU) + amm.getM(Xfee, expP)
            
            if expMl <= expMn and expMn < expMw:
                liq = False
            elif expMl <= expMn and expMn >= expMw:
                liq = False
                wait = False
            elif expMl < expMw:
                liq = False
        # When rL < math.floor(spot) we know rL <= spot <= rU
        # Therefore, we are trying to get fees and so it is fine to do liquidity order.
        # TODO: Check this
    else: #(when belP < spot)
        wait = False
        liq = False
        # Price is expected to go down. Therefore, a range order would only sell
        # since market orders are buying.
        # Thus range order not possible. Furthermore, now expensive to sell.
    
    if liq:
        with amm.capacity.request() as request:
            yield request
            kind = "X" if buy else "M"
            (X, M, contract) = amm.add(X, M, rL, rU, kind, True)
    if wait:
        yield env.timeout(time)
    with amm.capacity.request() as request:
        yield request
        if liq:
            (X, M, completionPer, filled) = retrieveFunds(amount, contract, name, belP, X, M)
            o.orders.append((name, spot, belP, amount, time, env.now-now, X, M, completionPer, "L"))
        elif buy:
            (X, M) = amm.tradeM(M)
            completionPer = -X/amount
            o.orders.append((name, spot, belP, amount, time, env.now-now, X, M, completionPer, "M"))
        else:
            (X, M) = amm.trade(X)
            completionPer = M/(amount*belP)
            o.orders.append((name, spot, belP, amount, time, env.now-now, X, M, completionPer, "M"))
    

def retrieveFunds(amount, contract, name, belP, assets=0, money=0):
    if amount < 0:
        (assets2, money, filled) = contract.retrieve(money, name)
        assets += assets2
        completionPer = -assets/amount
    else:
        (assets, money2, filled) = contract.retrieve(assets, name)
        money += money2
        completionPer = money/(amount*belP)
    return (assets, money, completionPer, filled)
    