#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 30 14:29:48 2022

@author: caren
"""
from Exchange import Exchange
import math
from Settings import minPriceRange, Settings
import simpy

class LOB(Exchange):
    def __init__(self, env, name, s):
        super().__init__(env, name, s)
        self.orders = [None] * (s.maxPriceRange - s.minPriceRange+1) # The buy side of the limit order book. Sell limit orders will end up here
        self.sell = [None] * (s.maxPriceRange - s.minPriceRange+1) # The sell side of the limit order book. Buy limit orders will end up here.
        initAsk = SellLimitOrder(s.initLiqSell, 1002, env.event(), self)        
        initBid = BuyLimitOrder(s.initLiqBuy, 998, env.event(), self)
        self.bestAsk = initAsk
        self.bestBid = initBid
        self.orders[getIndex(1002)] = initAsk
        self.orders[getIndex(998)] = initBid
    
    def __getliquidity__(self):
        totalAsked = self.getLiquidityUp()
        totalSold = self.getLiquidityDown()
        return (totalAsked, totalSold)
    
    # Get all liquidity that makes price go down, or all the assets you can sell
    def getLiquidityDown(self):
        bid = self.bestBid
        totalAsked = 0
        while bid:
            totalAsked += bid.assets
            bid = bid.next
        return totalAsked
    
    def getLiquidityUp(self):
        ask = self.bestAsk
        totalSold = 0
        while ask:
            totalSold += ask.assets
            ask = ask.next
        return totalSold
    
    def getLiquidity(self, price):
        liq = 0
        if self.bestAsk and price >= self.bestAsk.price:
            order = self.bestAsk
            while order and price >= order.price:
                liq += order.assets
                order = order.next
        elif self.bestBid.price and price <= self.bestBid.price:
            order = self.bestBid
            while order and price <= order.price:
                liq += order.assets
                order = order.next
        return liq
    
    def spot(self):
        if self.bestAsk and self.bestBid:
            return (self.bestAsk.price + self.bestBid.price)/2
        elif self.bestAsk and not self.bestBid:
            return 0
        elif self.bestBid and not self.bestAsk:
            return math.inf
        return -1
    
    def bestBidPrice(self):
        if self.bestBid:
            return self.bestBid.price 
        return self.s.minPriceRange

    def bestAskPrice(self):
        if self.bestAsk:
            return self.bestAsk.price
        return self.s.maxPriceRange
    
    def __bigSpread__(self):
        X = 1000
        (XAsk, MAsk) = self.__expMon__(self.bestAsk, X)
        (XBid, MBid) = self.__expMon__(self.bestBid, X)
        ask = MAsk/X if XAsk == 0 else self.s.maxPriceRange
        bid = MBid/X if XBid == 0 else self.s.minPriceRange
        return ask-bid
    
    def __expMon__(self, order, X):
        M = 0
        while order.next and X>0:
            if order.assets > X:
                M += X*order.price
                X = 0
            else:
                M += order.assets*order.price
                X -= order.assets
                order = order.next
        return (X, M)
    
    def add(self, assets, price, name=-1): #Name temporarily
        # Sell limit order
        if assets > 0:
            # TODO: Convert to a market order (party)
            if self.bestBid and price <= self.bestBid.price:
                return
            prev = None
            priceCheck = price
            while self.bestAsk and not prev and priceCheck >= self.bestAsk.price:
                prev = self.orders[getIndex(priceCheck)]
                priceCheck -= 1
            lo = SellLimitOrder(assets, price, self.env.event(), self, prev=prev, name=name)
            self.orders[getIndex(price)] = lo
            if not self.bestAsk:
                self.bestAsk = lo
            elif self.bestAsk.price > price:
                lo.next = self.bestAsk.next
                self.bestAsk.setprev(lo)
                self.bestAsk = lo
            return (0, 0, lo)
        # Buy limit order
        elif assets < 0:
            # TODO: Convert to a market order (partly)
            if self.bestAsk and price >= self.bestAsk.price:
                return
            prev = None
            priceCheck = price
            while self.bestBid and not prev and priceCheck <= self.bestBid.price:
                prev = self.orders[getIndex(priceCheck)]
                priceCheck +=1
            lo = BuyLimitOrder(-assets, price, self.env.event(), self, prev=prev, name=name)
            self.orders[getIndex(price)] = lo
            if not self.bestBid:
                self.bestBid = lo
            elif self.bestBid.price < price:
                lo.next = self.bestBid.next
                self.bestBid.setprev(lo)
                self.bestBid = lo
            return (0, 0, lo)
    
    # We trade money for assets, thus buying  
    def tradeM(self, amount):
        (totalBought, totalSold) = self.__getliquidity__() #Statistics
        self.allPrices.append((self.spot(), totalBought, totalSold))
        assets = 0
        while self.bestAsk and amount>0.00000001:
            (assetsAdd, amount) = self.bestAsk.get(amount)
            assets += assetsAdd
        self.marketBuyTransactions.push(assets + amount/self.s.maxPriceRange)
        
        self.addStatistics(sellVol = 0, buyVol = assets)
        
        return (assets, amount)
      
    # We assume amount>0, thus selling
    def trade(self, amount):
        orAmount = amount
        (totalBought, totalSold) = self.__getliquidity__()
        self.allPrices.append((self.spot(), totalBought, totalSold))
        money = 0
        while self.bestBid and amount>0.00000001:
            (amount, moneyAdd) = self.bestBid.get(amount)
            money += moneyAdd
        
        self.marketSellTransactions.push(orAmount - amount)
        self.addStatistics(sellVol = orAmount - amount, buyVol = 0)
           
        return (amount, money)

def getIndex(price):
    return price - minPriceRange

class LimitOrder(object):
    def __init__(self, assets, price, succes, lob, prev=None, name=-1):
        self.assetsOr = assets
        self.assets = assets
        self.price = price
        self.prev = None
        self.next = None
        self.name = name
        self.setprev(prev)
        self.succes = succes
        self.lob = lob
    
    def setAssets(self, assets):
        self.assets = max(assets, 0)
        if self.assets <= 0.00000001:
            self.cancel()
    
    def setnext(self, other):
        if other:
            self.next = other
            self.prev = other.prev
            other.prev = self
        if self.prev:
            self.prev.next = self
        return
    
    def setprev(self, other):
        if other:
            self.prev = other
            self.next = other.next
            other.next = self
        if self.next:
            self.next.prev = self
        return
    
    def cancel(self):
        if self.succes.triggered:
            return
        if self.prev:
            self.prev.next = self.next
        if self.next:
            self.next.prev = self.prev
        if self.lob.bestAsk == self:
            self.lob.bestAsk = self.next
        elif self.lob.bestBid == self:
            self.lob.bestBid = self.next
        if self.lob.orders[getIndex(self.price)] == self:
            if self.prev and self.prev.price == self.price:   
                self.lob.orders[getIndex(self.price)] = self.prev
            else:
                self.lob.orders[getIndex(self.price)] = None
        self.succes.succeed(self.assets)
        self.next = None
        self.prev = None
        self.assets = 0
    
    def retrieve(self, kind, amount, name = -1):
        if kind == "Money":
            money = (self.assetsOr-self.assets)*self.price
            (assets, money2) = self.lob.trade(self.assets + amount)
            return (assets, money+money2, self.assets==0)
        else:
            assets = self.assetsOr-self.assets
            money = self.assets*self.price
            (assets2, money) = self.lob.tradeM(money+ amount)
            return (assets+assets2, money, self.assets==0)

class BuyLimitOrder(LimitOrder):
    def __init__(self, assets, price, succes, lob, prev=None, name=-1):
        super().__init__(assets, price, succes, lob, prev, name)
    
    def get(self, assets):
        assetBuy = min(self.assets, assets)
        money = self.price * assetBuy
        self.setAssets(self.assets - assetBuy)
        return (assets-assetBuy, money)
    
    # Since this was a buy limit order, we know they want assets
    def retrieve(self, amount, name=-1):
        filled = self.assets == 0
        assets = self.assetsOr-self.assets
        money = self.assets*self.price
        (assets2, money) = self.lob.tradeM(money)
        self.setAssets(0)
        return (assets+assets2, money, filled)
        
class SellLimitOrder(LimitOrder):
    def __init__(self, assets, price, succes, lob, prev=None, name=-1):
        super().__init__(assets, price, succes, lob, prev, name)
    
    def get(self, money):
        assets = min(money / self.price, self.assets)
        money -= assets*self.price
        self.setAssets(self.assets-assets)
        return(assets, money)
    
    # Since this was a buy limit order, we know they want money
    def retrieve(self, amount, name=-1):
        filled = self.assets == 0
        money = (self.assetsOr-self.assets)*self.price
        (assets, money2) = self.lob.trade(self.assets)
        self.setAssets(0)
        return (assets, money+money2, filled)