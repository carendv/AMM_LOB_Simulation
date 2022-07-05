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

random.seed(10)
    
def trade(env, name, exchange, s, o, r):
    print(name)
    guessedTime = -1
    now = env.now
    spot = exchange.spot()
    ask = exchange.bestAskPrice()
    bid = exchange.bestBidPrice()
    o.prices.append(spot)
    o.spread.append(ask-bid)
    belP = round(s.trueP+(min(max(np.random.normal(0, 1), -3), 3))) if r.random() < s.infP else spot
    kind = None
    
    if belP != spot:
        buy = 0.95/(1+np.e**(-0.5*(belP-spot))) > r.random()
    else:
        buy = 0.5 < r.random()
    
    kind = None

    randomDraw = random.random()
    if isinstance(exchange, LOB) and exchange.getLiquidity(900) < 40000:
        kind = "L"
        buy = True
    elif isinstance(exchange, LOB) and exchange.getLiquidity(1200) < 40000:
        kind = "L"
        buy = False
    elif isinstance(exchange, AMM) and \
        (sum(list(map(abs, exchange.dL)))/2 < math.sqrt(40000000)*2 or \
         exchange.L < 20000000):
        kind = "L"
    elif buy and probOrder(belP, ask) > randomDraw:
        kind = "L"
    elif not buy and probOrder(bid, belP) > randomDraw:
        kind = "L"
    elif buy and probOrder(belP, ask) < randomDraw:
        kind = "M"
    elif not buy and probOrder(bid, belP) < randomDraw:
        kind = "M"
    
    (time, amount) = s.getBuy(kind) if buy else s.getSell(kind)
    
    # If we have some patience, we will look if we expect the limit order 
    # to succeed in time. If we don't expect it to, we still go for market order.
    if time > 0:
        # Get the quantity we expect to be traded via market orders in waiting time
        (buyQuan, sellQuan) = exchange.getBuySellQuantity(time)
        
        # Alter these numbers toward the direction of the expected price.
        # For instance, if we expect the price to increase by 10%, there should
        # be 10% more buy orders. However, we take the average of this and the 
        # current buy/sell quantity, since the history also takes price increases
        # into account.
        totQuan = buyQuan + sellQuan
        sellQuan = sellQuan/2 + totQuan/(1 + belP/spot)/2
        buyQuan = totQuan - sellQuan
        
        # For the AMM/LOB the decision of taking the limit/range order differs.
        # In the LOB, we look at two things: 
        # First, can we get what we want by waiting?
        # Second, does the market order after waiting probably better?
        # If either has answer yes, we do the limit order
        if isinstance(exchange, LOB):
            p = bestLimit(amount, belP, spot, exchange.bestAskPrice(), exchange.bestBidPrice())
            liqAtP = exchange.getLiquidity(p) + abs(amount)
            
            # If waiting gives a better price, and waiting is within time, do it
            if buy:
                liqAtPnew = exchange.getLiquidity(p-1) + abs(amount)
                if sellQuan - buyQuan >= liqAtPnew:
                    p -= 1
                    guessedTime = time*liqAtPnew/(sellQuan - buyQuan)
                else:
                    guessedTime = time*liqAtP/sellQuan
                    if liqAtP > sellQuan or belP > ask:
                        time = 0
            else:
                liqAtPnew = exchange.getLiquidity(p+1) + abs(amount)
                if buyQuan - sellQuan >= liqAtPnew:
                    p += 1
                    guessedTime = time*liqAtPnew/(buyQuan - sellQuan)
                else:
                    guessedTime = time*liqAtP/buyQuan
                    if liqAtP > buyQuan or belP < bid:
                        time = 0
                
            # if buy and sellQuan - buyQuan >= exchange.getLiquidity(p-1) + abs(amount):
            #     p -= 1
            # if not buy:
            #     liqAtPnew = exchange.getLiquidity(p+1) + abs(amount)
            # if not buy and buyQuan - sellQuan >= exchange.getLiquidity(p+1) + abs(amount):
            #     p += 1
            
            # if buy:
            #     guessedTime = time*(liqAtP + abs(amount))/sellQuan
            # else:
            #     guessedTime = time*(liqAtP + abs(amount))/buyQuan
            
            # if (liqAtP > sellQuan and buy) \
            #     and (liqAtP > buyQuan and not buy) \
            #     and (belP > spot and buy) and (belP < spot and not buy):
            #         time = 0
        # In the AMM, we compute the expected time after waiting.
        # If this is higher than our range, we have succesfully traded within the time
        # Else, if the price is better than the current, waiting is also the way to go.
        # TODO, incorporate fees?
        elif isinstance(exchange, AMM):
            (lower, upper) = bestRange(exchange, buy, belP, spot)
            expP = exchange.getPrice(sellQuan-buyQuan)
            assets = max(amount, 0)
            money = max(-belP*amount, 0)
            with exchange.capacity.request() as request:
                yield request
                expL = exchange.getExpLiq(assets, money, lower, upper)
            
            sL = math.sqrt(lower)
            sU = math.sqrt(upper)
            
            # If price not past range, compute the total price of this range order
            expAss = expL*(sU-expP)/(sU*expP)
            expMon = expL*(expP-sL)

            if ((expP < sU and not buy) or (expP > sL and buy)) \
                and ((expP < sU and not buy and expMon+expAss*(1-exchange.F)*belP < -amount*(1-exchange.F)*belP) \
                or (expP > sL and buy and expAss+expMon*(1-exchange.F)/belP < -amount)):
                    time = 0
    
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
            return 
        
    if isinstance(exchange, AMM):
        kind = "Assets" if amount < 0 else "Money"
        (lower, upper) = bestRange(exchange, buy, belP, spot)
        with exchange.capacity.request() as request:
            yield request
            (assets, money, contract) = exchange.add(assets, money, lower, upper, kind)
        notEnded = True
        while time > 0 and notEnded:
            yield env.timeout(min(60, time))
            time -= 60
            # You can get what you want
            if contract.expX() + assets > -amount:
                notEnded = False
    else:
        with exchange.capacity.request() as request:
            yield request
            (assets, money, contract) = exchange.add(amount, p, name) #Temporarily also name
        yield contract.succes | env.timeout(time)
    with exchange.capacity.request() as request:
        yield request
        if amount < 0:
            (assets2, money, filled) = contract.retrieve(money, name)
            assets += assets2
            completionPer = -assets/amount
        else:
            (assets, money2, filled) = contract.retrieve(assets, name)
            money += money2
            completionPer = money/(amount*belP)
    o.orders.append((name, spot, belP, amount, time, guessedTime, env.now-now, assets, money, completionPer, filled))
        
def bestRange(amm, buy, belP, spot):
    # Believed price is higher than market, we want to buy
    # Buy-limit undercut
    if belP > spot and buy:
        lower = math.floor(spot)
        upper = math.ceil(spot)
    # Believed price is higher than market, we want to sell
    # Sell limit at correct price
    elif belP > spot and not buy:
        lower = math.floor(spot)
        upper = amm.s.maxPriceRange
    # Believed price is lower than market, we want to sell
    # Sell limit undercut
    elif belP < spot and not buy:
        lower = math.floor(spot)
        upper = math.ceil(spot)
    # Believed price is lower than market, we want to buy
    # Buy limit at correct price
    elif belP < spot and buy:
        lower = amm.s.minPriceRange
        upper = math.ceil(spot)
    elif not buy:
        lower = np.floor(np.percentile(amm.prices, 5)**2)#np.floor(spot)
        upper = np.ceil(np.percentile(amm.prices, 95)**2)
    else:
        lower = np.floor(np.percentile(amm.prices, 5)**2)
        upper = np.ceil(np.percentile(amm.prices, 95)**2)#np.ceil(spot)
    
    if lower == upper:
        lower = max(lower - 1, amm.s.minPriceRange)
        upper = min(upper + 1, amm.s.maxPriceRange)
    if lower > spot:
        lower = math.floor(spot)
    if upper < spot:
        upper = math.ceil(spot)
    return (lower,upper)

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
                return min(bestBid+1, bestAsk-1)
            # If we believe the price to be lower, we don't want to pay extra.
            else:
                return belP
        # informed seller
        else:
            # Undercut if we believe the true price to be lower.
            # We still get a good deal since it is sold more expensive.
            if belP < bestAsk:
                return max(bestAsk-1, bestBid+1)
            # If we believe the price to be higher, we don't want to get less.
            else:
                return belP
    # Uninformed trader
    else:
        if assets < 0:
            return bestBid
        else:
            return bestAsk

def probOrder(belP, marketPrice):
    return 1/(1+np.e**(-(marketPrice-belP)/4))
