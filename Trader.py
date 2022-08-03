#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 16:07:47 2022

@author: caren
"""
import random
import math
from AMM import AMM
from LOB import LOB, getIndex
import numpy as np
from matplotlib import pyplot as plt

random.seed(10)
    
def trade(env, name, exchange, s, o, r):
    guessedTime = 0
    now = env.now
    spot = exchange.spot()
    ask = exchange.bestAskPrice()
    bid = exchange.bestBidPrice()
    informed = r.random() < s.infP
    belP = round(s.trueP+(min(max(np.random.normal(0, 1), -3), 3))) if informed else spot
    kind = None
    
    buy = 1/(1+np.e**(-s.agressiveness*(belP-spot))) > r.random()
    
    kind = None
    changed = False

    randomDraw = random.random()
    if isinstance(exchange, LOB):
        ret = randomCheckOrderLOB(exchange, belP, s)
    elif isinstance(exchange, AMM):
        ret = randomCheckOrderAMM(exchange, belP, s)

    if isinstance(exchange, AMM) and ret:
        (time, amount, lower, upper, kindOrder, changed, guessedTime) = ret
        buy = amount < 0
    elif isinstance(exchange, LOB) and ret:
        (time, amount, p, kindOrder, changed, guessedTime) = ret
        buy = amount < 0
    else:
        #TODO: shuffle the order at which conditions are looked at.
        if buy and probOrder(belP, ask) > randomDraw: 
            kind = "L"
        elif not buy and probOrder(bid, belP) > randomDraw:
            kind = "L"
        elif buy and probOrder(belP, ask) < randomDraw:
            kind = "M"
        elif not buy and probOrder(bid, belP) < randomDraw:
            kind = "M"
        if isinstance(exchange, LOB):
            (time, amount, p, changed, kindOrder, guessedTime) = freeOrderLOB(buy, kind, belP, s, env, exchange)
        else:
            (time, amount, lower, upper, changed, kindOrder, guessedTime) = yield env.process(freeOrderAMM(buy, kind, belP, s, env, exchange))
   
    # Market order
    if time == 0:
        with exchange.capacity.request() as request:
            yield request
            if amount < 0:
                amountM = -amount*belP
                (assets, money) = exchange.tradeM(amountM)
                completionPer = -assets/amount
            else:
                (assets, money) = exchange.trade(amount)
                completionPer = money/(amount*belP)
            o.orders.append((name, spot, belP, amount, time, guessedTime, env.now-now, assets, money, completionPer, "M"))
            o.strats.append(("M", informed, changed))
            return 
        
    if isinstance(exchange, AMM):
        kind = "Assets" if amount < 0 else "Money"
        assets = max(amount, 0)
        money = max(-belP*amount, 0)
        with exchange.capacity.request() as request:
            yield request
            (assets, money, contract) = exchange.add(assets, money, lower, upper, kind)
            o.strats.append((kindOrder, informed, changed, lower, upper))
        notEnded = True
        endTime = now + time
        while endTime >= env.now and notEnded:
            yield env.timeout(min(60, endTime-env.now))
            # You can get what you want
            if buy and (contract.expX(money) + assets > -amount or endTime == env.now):
                with exchange.capacity.request() as request:
                    yield request
                    if contract.expX(money) + assets > -amount or endTime == env.now:
                        (assets, money, completionPer, filled) = retrieveFunds(amount, contract, name, belP, assets, money)
                        notEnded = False
            elif not buy and (contract.expM(assets) + money > amount*belP or endTime == env.now):
                with exchange.capacity.request() as request:
                    yield request
                    if contract.expM(assets) + money > amount*belP or endTime == env.now:
                        (assets, money, completionPer, filled) = retrieveFunds(amount, contract, name, belP, assets, money)
                        notEnded = False
    else:
        with exchange.capacity.request() as request:
            yield request
            (assets, money, contract) = exchange.add(amount, p, name)
            o.strats.append((kindOrder, informed, changed, p))
        yield contract.succes | env.timeout(time)
        with exchange.capacity.request() as request:
            yield request
            (assets, money, completionPer, filled) = retrieveFunds(amount, contract, name, belP)

    o.orders.append((name, spot, belP, amount, time, guessedTime, env.now-now, assets, money, completionPer, filled))

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
        
    

def forcedLiquidityOrderAMM(buy, belP, amm, s):
    (lower, upper, kindOrder) = bestRange(amm, buy, belP, amm.spot())
    (time, amount) = s.getBuy("L") if buy else s.getSell("L")
    return (time, amount, lower, upper, kindOrder, False, -1)
    
def freeOrderAMM(buy, kind, belP, s, env, exchange):
    (time, amount) = s.getBuy(kind) if buy else s.getSell(kind)        
    (lower, upper, kindOrder) = bestRange(exchange, buy, belP, exchange.spot())
    spot = exchange.spot()
    changed = False
    # If we have some patience, we will look if we expect the limit order 
    # to succeed in time. If we don't expect it to, we still go for market order.
    if time > 0:
        # Get the quantity we expect to be traded via market orders in waiting time
        (buyQuan, sellQuan) = exchange.getBuySellQuantity(time)
        
        # Alter these numbers toward the direction of the expected price.
        (buyQuan, sellQuan) = updateBuySellQuantity(buyQuan, sellQuan, belP, spot)

        # In the AMM, we compute the expected price after waiting.
        # If this is higher than our range, we have succesfully traded within the time
        # Else, if the price is better than the current, waiting is also the way to go.
        limit = False
        
        # If the price doesn't change, at least you'll earn fees. So limit order is better
        if sellQuan == buyQuan:
            limit = True
        
        assets = max(amount, 0)
        money = max(-belP*amount, 0)
        
        checkpoints = [0.01, 0.05, 0.1, 0.25, 0.5, 1]
        
        for i in checkpoints:
            with exchange.capacity.request() as request:
                yield request
                (L, fX, fM, expP) = exchange.getExpLiq(assets, money, lower, upper, (sellQuan-buyQuan)*i)
            
            sL = math.sqrt(lower)
            sU = math.sqrt(upper)
            
            # If price not past range, compute the total price of this range order
            expAss = L*(sU-expP)/(sU*expP) + fX/env.now*time*L*i
            expMon = L*(expP-sL) + fM/env.now*time*L*i
            
            # You expect the range order to succeed if price is past limit
            if buy and expP <= sL:
                limit = True
            elif not buy and expP >= sU:
                limit = True
            elif buy and expP > sU:
                limit = False
            elif not buy and expP < sL:
                limit = False
            elif buy and expAss+exchange.expM(expMon, expP) > exchange.expM(-belP*amount):
                limit = True
            elif not buy and expMon+exchange.exp(expAss, expP) > exchange.exp(amount):
                limit = True
                
            if limit:
                break
        
        if not limit:
            time = 0
            kindOrder = "M"
            changed = True
    return (time, amount, lower, upper, changed, kindOrder, -1)


def forcedLiquidityOrderLOB(buy, belP, lob, s):
    (time, amount) = s.getBuy("L") if buy else s.getSell("L")
    (p, kindOrder) = bestLimit(amount, belP, lob.spot(), lob.bestAskPrice(), lob.bestBidPrice()) #(assets, belP, spot, bestAsk, bestBid)
    return (time, amount, p, kindOrder, False, -1)

def freeOrderLOB(buy, kind, belP, s, env, exchange):
    (time, amount) = s.getBuy(kind) if buy else s.getSell(kind)
    spot = exchange.spot()
    ask = exchange.bestAskPrice()
    bid = exchange.bestBidPrice()
    (p, kindOrder) = bestLimit(amount, belP, spot, ask, bid)
    guessedTime = 0
    changed = False
    
    # If we have some patience, we will look if we expect the limit order 
    # to succeed in time. If we don't expect it to, we still go for market order.
    if time > 0:
        # Get the quantity we expect to be traded via market orders in waiting time
        (buyQuan, sellQuan) = exchange.getBuySellQuantity(time)
        
        # Alter these numbers toward the direction of the expected price.
        (buyQuan, sellQuan) = updateBuySellQuantity(buyQuan, sellQuan, belP, spot)
    
        # For the AMM/LOB the decision of taking the limit/range order differs.
        # In the LOB, we look at two things: 
        # First, can we get what we want by waiting?
        # Second, does the market order after waiting probably better?
        # If either has answer yes, we do the limit order
        liqAtP = exchange.getLiquidity(p) + abs(amount)
        
        # If waiting gives a better price, and waiting is within time, do it
        if buy:
            liqAtPnew = exchange.getLiquidity(p-1) + abs(amount)
            if sellQuan - buyQuan >= liqAtPnew and p-1<bid: # Not at best bid
                p -= 1
                guessedTime = time*liqAtPnew/(sellQuan - buyQuan)
            elif sellQuan >= liqAtPnew and p>=bid: # at best bid
                p-= 1
                guessedTime = time*liqAtPnew/sellQuan
            elif sellQuan - buyQuan >= liqAtP and p<bid:
                guessedTime = time*liqAtP/(sellQuan - buyQuan)
            elif sellQuan >= liqAtP and p >=bid:
                guessedTime = time*liqAtP/sellQuan
            elif belP <= ask:
                guessedTime = time
            else:
                time = 0
                kindOrder = "M"
                changed = True
                guessedTime = 0
        else:
            liqAtPnew = exchange.getLiquidity(p+1) + abs(amount)
            if buyQuan - sellQuan >= liqAtPnew and p+1 > ask:
                p += 1
                guessedTime = time*liqAtPnew/(buyQuan - sellQuan)
            elif buyQuan >= liqAtPnew and p <= ask:
                p += 1
                guessedTime = time*liqAtPnew/buyQuan
            elif buyQuan - sellQuan >= liqAtP and p > ask:
                guessedTime = time*liqAtP/(buyQuan - sellQuan)
            elif buyQuan >= liqAtP and p <= ask:
                guessedTime = time*liqAtP/buyQuan
            elif belP >= bid:
                guessedTime = time
            else:
                time = 0
                kindOrder = "M"
                changed = True
                guessedTime = 0
    return (time, amount, p, changed, kindOrder, guessedTime)
        
def bestRange(amm, buy, belP, spot):
    # Believed price is higher than market, we want to buy
    # Buy-limit undercut
    if belP > spot and buy:
        lower = math.floor(spot)
        upper = math.ceil(spot)
        if lower == upper:
            upper +=1
        kind = "L-Buy-U"
    # Believed price is higher than market, we want to sell
    # Sell limit at correct price
    elif belP > spot and not buy:
        lower = amm.getLiquidityStart(belP)
        upper = amm.s.maxPriceRange
        if upper == lower:
            lower -= 1
        kind = "L-Sell-A"
    # Believed price is lower than market, we want to sell
    # Sell limit undercut
    elif belP < spot and not buy:
        lower = math.floor(spot)
        upper = math.ceil(spot)
        if lower == upper:
            lower -= 1
        kind = "L-Sell-U"
    # Believed price is lower than market, we want to buy
    # Buy limit at correct price
    elif belP < spot and buy:
        upper = amm.getLiquidityStart(belP)
        lower = amm.s.minPriceRange
        if upper == lower:
            upper += 1
        
        kind = "L-Buy-A"
    else:
        lower = int(min(np.floor(np.percentile(amm.prices, 5)**2), np.floor(spot)-2))
        upper = int(max(np.ceil(np.percentile(amm.prices, 95)**2), np.ceil(spot)+2))
        kind = "L-row"
    
    if lower == upper:
        lower -= 1
        upper += 1
    
    lower = int(max(lower, amm.s.minPriceRange))   
    upper = int(min(upper, amm.s.maxPriceRange))
    
    
    return (lower,upper, kind)

# Assets are the number of assets we sell to the market, thus. 
# Assets < 0: Buy from market
# Assets > 0: Sell to market
def bestLimit(assets, belP, spot, bestAsk, bestBid):
    # Informed trader
    if belP != spot:
        # Informed buyer
        if assets < 0:
            # Undercut if we believe the true price to be higher.
            # We still get a good deal since it is cheap.
            if belP > bestBid:
                return (min(bestBid+1, bestAsk-1), "L-Buy-U")
            # If we believe the price to be lower, we don't want to pay extra.
            else:
                return (belP, "L-Buy-A")
        # informed seller
        else:
            # Undercut if we believe the true price to be lower.
            # We still get a good deal since it is sold more expensive.
            if belP < bestAsk:
                return (max(bestAsk-1, bestBid+1), "L-Sell-U")
            # If we believe the price to be higher, we don't want to get less.
            else:
                return (belP, "L-Sell-A")
    # Uninformed trader
    else:
        if assets < 0:
            return (bestBid, "L-Buy-Row")
        else:
            return (bestAsk, "L-Sell-Row")

def probOrder(belP, marketPrice):
    return 1/(1+np.e**(-(marketPrice-belP)/1))

def randomCheckOrderAMM(exchange, belP, s):
    lst_condition_result = [(exchange.getLiquidityDown() < s.minLiquidity, True), (exchange.getLiquidityUp() < s.minLiquidity, False)]
    random.shuffle(lst_condition_result)
    
    for condition, param in lst_condition_result:
        if condition:
            return forcedLiquidityOrderAMM(param, belP, exchange, s)

def randomCheckOrderLOB(exchange, belP, s):
    lst_condition_result = [(exchange.getLiquidityDown() < s.minLiquidity, True), (exchange.getLiquidityUp() < s.minLiquidity, False)]
    random.shuffle(lst_condition_result)
    
    for condition, param in lst_condition_result:
        if condition:
            return forcedLiquidityOrderLOB(param, belP, exchange, s)

def updateBuySellQuantity(buyQuan, sellQuan, belP, spot):
    totQuan = buyQuan + sellQuan
    p = belP/spot
    weight = (1+abs(1-p))/2
    buyQuan = buyQuan*(1-weight) + totQuan*p/(1 + p)*weight
    sellQuan = totQuan - buyQuan
    return (buyQuan, sellQuan)
    