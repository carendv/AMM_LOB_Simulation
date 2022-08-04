#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 14:38:49 2022

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
    ret = randomCheckOrderLOB(exchange, belP, s)

    if ret:
        (time, amount, p, kindOrder, changed) = ret
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
        (time, amount, p, changed, kindOrder) = freeOrderLOB(buy, kind, belP, s, env, exchange)
   
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
            o.orders.append((name, spot, belP, amount, time, env.now-now, assets, money, completionPer, "M"))
            o.strats.append(("M", informed, changed))
            return 
        
    with exchange.capacity.request() as request:
        yield request
        (assets, money, contract) = exchange.add(amount, p, name)
        o.strats.append((kindOrder, informed, changed, p))
    yield contract.succes | env.timeout(time)
    with exchange.capacity.request() as request:
        yield request
        (assets, money, completionPer, filled) = retrieveFunds(amount, contract, name, belP)

    o.orders.append((name, spot, belP, amount, time, env.now-now, assets, money, completionPer, filled))

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

def forcedLiquidityOrderLOB(buy, belP, lob, s):
    (time, amount) = s.getBuy("L") if buy else s.getSell("L")
    (p, kindOrder) = bestLimit(amount, belP, lob.spot(), lob.bestAskPrice(), lob.bestBidPrice()) #(assets, belP, spot, bestAsk, bestBid)
    return (time, amount, p, kindOrder, False)

def freeOrderLOB(buy, kind, belP, s, env, exchange):
    (time, amount) = s.getBuy(kind) if buy else s.getSell(kind)
    spot = exchange.spot()
    ask = exchange.bestAskPrice()
    bid = exchange.bestBidPrice()
    (p, kindOrder) = bestLimit(amount, belP, spot, ask, bid)
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
            elif sellQuan >= liqAtPnew and p>=bid: # at best bid
                p-= 1
            elif sellQuan - buyQuan >= liqAtP and p<bid:
                pass
            elif sellQuan >= liqAtP and p >=bid:
                pass
            elif belP <= ask:
                pass
            else:
                time = 0
                kindOrder = "M"
                changed = True
        else:
            liqAtPnew = exchange.getLiquidity(p+1) + abs(amount)
            if buyQuan - sellQuan >= liqAtPnew and p+1 > ask:
                p += 1
            elif buyQuan >= liqAtPnew and p <= ask:
                p += 1
            elif buyQuan - sellQuan >= liqAtP and p > ask:
                pass
            elif buyQuan >= liqAtP and p <= ask:
                pass
            elif belP >= bid:
                pass
            else:
                time = 0
                kindOrder = "M"
                changed = True
    return (time, amount, p, changed, kindOrder)

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
    