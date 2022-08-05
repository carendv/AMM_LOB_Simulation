#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 13:15:34 2022

@author: caren
"""
import random
import math
from Trader import updateBuySellQuantity, retrieveFunds
import numpy as np

random.seed(10)
    
def trade(env, name, amm, s, o, r):
    #print(name)
    now = env.now
    forced = False
    
    ###################
    # Liquidity check #
    ###################    
    lst_condition_result = [(amm.getLiquidityDown() < s.minLiquidity, True), (amm.getLiquidityUp() < s.minLiquidity, False)]
    random.shuffle(lst_condition_result)
    
    for condition, param in lst_condition_result:
        if condition:
            forced = True
            buy = param
    
    ##################
    # Initialization #
    ##################
    spot = amm.spot()
    informed = r.random() < s.infP
    belP = round(s.trueP+(min(max(np.random.normal(0, 1), -3), 3))) if informed and not forced else spot
    buy = 1/(1+np.e**(-s.agressiveness*(belP-spot))) > r.random() if not forced else buy
    kind = "L" if forced else False
    (time, amount) = s.getBuy(kind) if buy else s.getSell(kind)
    
    ######################
    # Decision variables #
    ######################
    wait = False
    liq = time > 0
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
        rL = max(min(math.ceil(expP), math.floor(spot)-2), math.floor(belP)-1) # Use here 95% interval
        rU = max(amm.getLiquidityStart(rL), rL+2)
        
        # Compute expected fees
        if rU <= math.floor(spot) and not forced:
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
                wait = True
            elif expXl <= expXn and expXn >= expXw:
                liq = False
                wait = False
            elif expXl < expXw:
                liq = False
                wait = True
        elif rU > math.floor(spot) and not forced:
            liq=True
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
        if rL >= math.ceil(spot) and not forced:
            buyInRange = amm.getMPrice(rU)-amm.getMPrice(rL)
            perInRange = buyInRange/(buyQuan-sellQuan)
            Lself = math.sqrt(rL*rU)/(math.sqrt(rU)-math.sqrt(rL))*X*(rU-rL)
            LOther = amm.getL(rL, rU)
            
            Xfee = sellQuan*perInRange*Lself/(Lself+LOther)
            Mfee = buyQuan*perInRange*math.sqrt(rL*expP)*Lself/(Lself+LOther)
            
            expMl = Mfee + X*math.sqrt(rL*rU) + amm.getM(Xfee, expP)
            
            if expMl <= expMn and expMn < expMw:
                liq = False
                wait = True
            elif expMl <= expMn and expMn >= expMw:
                liq = False
                wait = False
            elif expMl < expMw:
                liq = False
                wait = True
        # When rL < math.floor(spot) we know rL <= spot <= rU
        # Therefore, we are trying to get fees and so it is fine to do liquidity order.
        # TODO: Check this
    else: #(when belP < spot)
        wait = False
        liq = False
        # Price is expected to go down. Therefore, a range order would only sell
        # since market orders are buying.
        # Thus range order not possible. Furthermore, now expensive to sell.
    
    
    ###################
    # Order execution #
    ###################  
    if liq:
        with amm.capacity.request() as request:
            yield request
            kind = "X" if buy else "M"
            (X, M, contract) = amm.add(X, M, rL, rU, kind, False)
        filled = False
        while env.now < now + time and not filled:
            yield env.timeout(min(100, now+time-env.now))
            if buy and amm.spot() < rL:
                filled = True
            elif not buy and amm.spot() >rU:
                filled = True
    if wait:
        yield env.timeout(time)
    with amm.capacity.request() as request:
        yield request
        if liq:
            (X, M, completionPer, filled) = retrieveFunds(amount, contract, name, belP, X, M)
            o.orders.append((name, spot, amm.spot(), expP, belP, amount, time, env.now-now, X, M, completionPer, "L"))
        elif buy:
            (X, M) = amm.tradeM(M)
            completionPer = (-amount*belP/X)/belP
            o.orders.append((name, spot, amm.spot(), expP, belP, amount, time, env.now-now, X, M, completionPer, "M"))
        else:
            (X, M) = amm.trade(X)
            completionPer = M/amount/belP
            o.orders.append((name, spot, amm.spot(), expP, belP, amount, time, env.now-now, X, M, completionPer, "M"))
    