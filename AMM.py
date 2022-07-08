#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 10:27:39 2022

@author: caren
"""
import math
from Exchange import Exchange
import simpy
from Settings import Settings
from collections import deque
from matplotlib import pyplot as plt
import random
import numpy as np

class AMM(Exchange):
    def __init__(self, env, name, s):
        super().__init__(env, name, s)
        s.initLiqBuy = s.initLiqSell
        
        # Compute liquidity of pool, given initial buy liquidity and price
        pb = math.sqrt(s.AMMmax)
        a = (1-s.trueP/s.AMMmax)
        b = -2*s.initLiqSell*s.trueP/pb
        c = -s.trueP*s.initLiqSell**2
        self.L = round((-b+math.sqrt(b**2-4*a*c))/(2*a))
        
        # Compute lower bound of initial range order given initial sell and pool liquidity
        a = round( ((s.trueP*(s.initLiqSell+self.L/pb) - s.initLiqBuy*s.trueP)/self.L)**2 )
        
        # Update price for rounding errors
        self.sP = math.sqrt( (s.initLiqBuy*s.trueP + self.L*math.sqrt(a))/(s.initLiqSell + self.L/pb) )
        
        
        # Variables needed for updates
        self.dL = [0] * (s.maxPriceRange - s.minPriceRange+1)
        self.nR = [0] * (s.maxPriceRange - s.minPriceRange+1)
        
        # Variables needed to keep track of fees
        self.feeGrowthX = [0] * (s.maxPriceRange - s.minPriceRange+1)
        self.feeGrowthM = [0] * (s.maxPriceRange - s.minPriceRange+1)
        self.F = s.fee 
        
        # Initialize first liquidity
        self.dL[a-s.minPriceRange] = self.L
        self.dL[s.AMMmax-s.minPriceRange] = -self.L
        self.nR[a-s.minPriceRange] += 1 # Number of references to price tick
        self.nR[s.AMMmax-s.minPriceRange] += 1 # Number of references to price tick
        self.index = a-s.minPriceRange
        
        # Additional statistics calculation needed for fee calculations
        self.prices = deque(maxlen=s.lookNTransactionsBack)
        self.prices.append(self.sP)
        self.tradePrices = []
        
    def __getVirtual__(self):
        X = self.L/self.sP
        M = self.L*self.sP
        return (X, M)
    
    def __getReal__(self):
        (Xvir, Mvir) = self.__getVirtual__()
        (_, pb) = self.__getTickWindowUp__()
        (_, pa) = self.__getTickWindowDown__()
        X = Xvir-self.L/pb
        M = Mvir-self.L*pa
        return (X,M)
    
    def spot(self):
        return self.sP*self.sP
    
    def __setIndex__(self, index):
        self.index = index
        self.__setL__()
    
    def __setL__(self):
        L1 = round(np.sum(self.dL[:self.index+1]))
        L2 = -round(np.sum(self.dL[self.index+1:]))
        self.L = min(L1, L2)
        
    def bestBidPrice(self):
        return self.exp(1)

    def bestAskPrice(self):
        return self.exp(-1)
    
    
    # Trade X for M, as seen from the AMM. This means:
    # A negative amount of X means that we buy it from the AMM
    # A positive amount of X means that we sell it to the AMM
    # A negative amount of M returned means that we pay it to the AMM
    # A positive amount of M returned means that we get it from the AMM
    def trade(self, X):
        orX = X
        self.prices.append(self.sP)
        self.allPrices.append(self.sP)
        
        # Only full units of assets are accepted, though we do not check
        M = 0
        (inL, inU, pL, pU) = self.__getTickWindow__()
        
        
        # If the trade crosses it's upper price (only possible when buying),
        # split up the trade to go to the crossing
        while X < 0:
            if self.sP < pU and self.L == 0:
                break
            elif self.sP == pU and inU == len(self.nR)-1:
                break
            elif self.sP == pU:
                self.__setIndex__(inU)
            else:
                pN = self.__getPricedX__(X) if self.index > -1 else math.inf
                
                # Check if we don't cross a border
                # If we don't, just execute the order
                # If we do (computed price higher than border), split up
                if pN <= pU:
                    dM = self.__getdMpN__(pN)/(1-self.F)
                    M -= dM
                    self.sP = pN
                    self.feeGrowthM[inU] += abs(dM*self.F/self.L)
                    self.marketBuyTransactions.push(X)
                    X = 0
                else: 
                    # Update all returns of user for partial swap
                    dM = self.__getdMpN__(pU)/(1-self.F)
                    M -= dM
                    dX = self.__getdXpN__(pU)
                    X -= dX
                    
                    # Update liquidity
                    self.sP = pU
                    if dM > 0:
                        self.feeGrowthM[inU] += abs(dM*self.F/self.L)
                        self.marketBuyTransactions.push(dX)
                    self.__setIndex__(inU)
            (inL, inU, pL, pU) = self.__getTickWindow__()
        
        # If the trade crosses it's lower price (only possible when selling),
        # split up the trade to go to the crossing
        while X > 0 and (self.index >= 0 or self.L>0):
            if self.sP > pL and self.L == 0:
                break
            elif self.sP == pL and self.index == 0:
                break
            elif self.sP == pL:
                self.__setIndex__(inL)
            else:
                pN = self.__getPricedX__(X*(1-self.F)) if inU < len(self.nR) else 0
                
                if pN >= pL:
                    M -= self.__getdMpN__(pN)
                    self.sP = pN
                    self.feeGrowthX[inU] += abs(X*self.F/self.L)
                    self.marketSellTransactions.push(X)
                    X=0
                else:
                    M -= self.__getdMpN__(pL)
                    dX = self.__getdXpN__(pL)/(1-self.F)
                    X -= dX
                    self.sP = pL
                    if inU < len(self.nR) and not dX == 0:
                        self.feeGrowthX[inU] += abs(dX*self.F/self.L)
                        self.marketSellTransactions.push(dX)
                    self.__setIndex__(inL)
            (inL, inU, pL, pU) = self.__getTickWindow__()
        
        self.tradePrices.append(M/(orX-X)) if not orX-X == 0 else None
        if orX > 0:
            self.addStatistics(sellVol = orX-X, buyVol = 0)
        else: 
            self.addStatistics(sellVol = 0, buyVol = X-orX)
        
        return (X, M)  
        
    # Trade M for X, as seen from the AMM. This means:
    # A negative amount of M means that we sell X untill we get this
    # A positive amount of M means that we buy as much X as we can
    # A negative amount of X returned means that we pay it to the AMM
    # A positive amount of X returned means that we get it from the AMM
    def tradeM(self, M):
        orM = M
        self.prices.append(self.sP)
        self.allPrices.append(self.sP)
        
        # Only full units of assets are accepted, though we do not check
        X = 0
        (inL, inU, pL, pU) = self.__getTickWindow__()
        
        
        # The price lowers.
        # This is only possible when we sell X, and thus retrieve money (M<0)
        # From the AMM.
        # If the lower border is crossed, trade up to the crossing and split
        while M < 0:
            if self.sP > pL and self.L == 0:
                break
            elif self.sP == pL and self.index == 0:
                break
            elif self.sP == pL:
                self.__setIndex__(inL)
            else:
                pN = self.__getPricedM__(M) if self.index > -1 else math.inf
                
                # Check if we don't cross a border
                # If we don't, just execute the order
                # If we do (computed price higher than border), split up
                if pN >= pL:
                    M = 0
                    dX = self.__getdXpN__(pN)/(1-self.F)
                    X -= dX
                    self.sP = pN
                    self.feeGrowthX[inU] += abs(dX*self.F/self.L)
                    self.marketSellTransactions.push(dX)
                else: 
                    # Update all returns of user for partial swap
                    dX = self.__getdXpN__(pL)/(1-self.F)
                    X -= dX
                    M -= self.__getdMpN__(pL)
                    
                    # Update liquidity
                    self.sP = pL
                    if dX > 0:
                        self.feeGrowthX[inU] += abs(dX*self.F/self.L)
                        self.marketSellTransactions.push(dX)
                    self.__setIndex__(inL)
            (inL, inU, pL, pU) = self.__getTickWindow__()
        
        # If the trade crosses it's upper price (only possible when selling),
        # split up the trade to go to the crossing
        while M > 0 and (inU < len(self.nR) or self.L>0):
            if self.sP < pU and self.L == 0:
                break
            elif self.sP == pU and inU == len(self.nR)-1:
                break
            elif self.sP == pU:
                self.__setIndex__(inU)
            else:
                pN = self.__getPricedM__(M*(1-self.F)) if inU < len(self.nR) else 0
                
                if pN <= pU:
                    dX = self.__getdXpN__(pN)
                    X -= dX
                    self.sP = pN
                    self.feeGrowthM[inU] += abs(M*self.F/self.L)
                    self.marketBuyTransactions.push(-dX)
                    M=0
                else:
                    dX = self.__getdXpN__(pU)
                    X -= dX
                    dM = self.__getdMpN__(pU)/(1-self.F)
                    M -= dM
                    self.sP = pU
                    if inU < len(self.nR) and not dM == 0:
                        self.feeGrowthM[inU] += abs(dM*self.F/self.L)
                        self.marketBuyTransactions.push(dX)
                    self.__setIndex__(inU)
            (inL, inU, pL, pU) = self.__getTickWindow__()
        
        self.tradePrices.append((orM-M)/(X)) if not orM-M == 0 else None
        
        if orM > 0:
            self.addStatistics(sellVol = 0, buyVol = X)
        else: 
            self.addStatistics(sellVol = X, buyVol = 0)
        
        return (X, M)  
    
    # Expected money for the number of assets (X).
    # X < 0: buying from AMM
    # X > 0: Selling to AMM
    # M: Money you get from the AMM
    def exp(self, X): 
        p = self.sP
        L = self.L
        index = self.index
        M = 0
        while X<0 and L>0:
            (index, pU) = self.__getTickWindowUp__(index)
            pN = 1/(1/p + X/L)
            if pN < pU and pN > 0:
                M += L*(pN-p)
                M = M/(1-self.F)
                X=0
                break
            else:
                M += L*(pU-p)
                X -= (1/pU - 1/p)*L
                p = pU
                L+= self.dL[index] 
        X = X*(1-self.F)
        while X > 0 and L>0:
            (indexL, pL) = self.__getTickWindowDown__(index) # Easy fix for the fact that index saved is lower index
            pN = 1/(1/p + X/L)
            if pN > pL:
                M -= L*(pN-p)
                X = 0
                break
            elif index < 0:
                M -= L*(pL-p)
                X -= (1/pL - 1/p)*L
                p = pL
                break
            else:
                M -= L*(pL-p)
                X -= (1/pL - 1/p)*L
                p = pL
                L -= self.dL[index]
                index = indexL
        
        if X < 0:
            return math.inf
        elif X > 0:
            return 0
        else:
            return M
        
    # TODO: Check function
    def getPrice(self, assets, index = None, sP = None, L = None):
        if not index:
            index = self.index
        if not sP:
            sP = self.sP
        if not L:
            L = self.L
            
        (inL, pL) = self.__getTickWindowDown__(index)
        (inU, pU) = self.__getTickWindowUp__(index)
        
        while assets > 0 and (index >= 0 or L>0):
            if sP > pL and L == 0:
                return sP
            elif sP == pL and index <= 0:
                return sP
            elif sP == pL:
                L = round(L-self.dL[index])
                index = inL
            else:
                pN = self.__getPricedX__(assets*(1-self.F)) if inU < len(self.nR) else 0
                
                if pN >= pL:
                    return pN
                else:
                    assets -= self.__getdXpN__(pL)/(1-self.F)
                    sP = pL
                    index = inL
                    L = round(L-self.dL[index])
            (inL, pL) = self.__getTickWindowDown__(index)
            
        if assets > 0:
            return math.sqrt(self.s.minPriceRange)
        
        while assets < 0:
            if sP < pU and L == 0:
                return sP
            elif sP == pU and inU == len(self.nR)-1:
                return sP
            elif sP == pU:
                L += self.dL[pU]
                index = pU
            else:
                pN = self.__getPricedX__(assets) if index > -1 else math.inf
                if pN <= pU:
                    return pN
                else: 
                    assets -= self.__getdXpN__(pU)
                    sP = pU
                    L += self.dL[inU]
                    index = inU
            (inU, pU) = self.__getTickWindowUp__(index)
        
        return math.sqrt(self.s.maxPriceRange)
    
    # TODO: check function
    def getExpLiq(self, assets, money, lower, upper):
        f = self.F
        index = self.index
        sP = self.sP
        self.F = 0
        (_, _, C) = self.add(assets, money, lower, upper)
        L = C.L
        C.retrieve()
        self.F = f
        self.sP = sP
        self.__setIndex__(index)
        return L
    
    def __getTickWindow__(self):
        (inU, pU) = self.__getTickWindowUp__()
        (inL, pL) = self.__getTickWindowDown__()
        return (inL, inU, pL, pU)
    
    # Returns pair: index upperbound, sqrt price upper bound
    def __getTickWindowUp__(self, inU=None):
        if not inU:
            inU = self.index
        
        inU +=1
        while inU < len(self.nR)-1 and self.nR[inU] == 0:
            inU += 1
        pu = math.sqrt(inU+self.s.minPriceRange)
        return (inU, pu)
    
    # Returns pair: index next lowerbound, sqrt price lower bound
    def __getTickWindowDown__(self, index = None):
        if not index:
            index = self.index
        
        inL = index
        
        if inL == 0:
            return (0, 0)
        
        inL -=1
        while inL >= 0 and self.nR[inL] == 0:
            inL -= 1
        
        if self.index == len(self.nR)-1:
            return (inL, math.inf)
        else:
            return (inL, math.sqrt(index+self.s.minPriceRange))
    
    def __getPricedX__(self, dX): return 1/(1/self.sP + dX/self.L)
    def __getPricedM__(self, dM): return dM/self.L + self.sP
    def __getdMpN__(self, pN): return self.L*(pN-self.sP)
    def __getdXpN__(self, pN): return (1/pN - 1/self.sP)*self.L
    
    # In this simulation we assume that we can only add X or M, not both
    # Return X, M, L
    def add(self, X, M, a, b, kind=None): 
        # Compute some prices we need
        pa = math.sqrt(a) # Get sqrt(pa)
        pb = math.sqrt(b) # Get sqrt(pb)
        (fX, fM) = self.__retrieveFeesPerUnitInRange__(pa, pb)
        
        if pb <= self.sP: # Only money
            L = round(M/(pb-pa))
            self.__add__(pa, pb, L)
            return (X, 0, Contract(pa, pb, L, fX, fM, self, kind))
        elif pa >= self.sP: # Only asset
            L = round(X*pa*pb/(pb-pa))
            self.__add__(pa, pb, L)
            return (0, M, Contract(pa, pb, L, fX, fM, self, kind))
        
        (inL, inU, pm, pp) = self.__getTickWindow__()
        
        # If no liquidity to trade, then check if you can get liquidity
        if self.L == 0:
            if self.sP == pm and X == 0:
                return (X, M, Contract(pa, pb, 0, fX, fM, self, kind))
            if pm == math.inf and X == 0:
                return (X, M, Contract(pa, pb, 0, fX, fM, self, kind))
            if self.sP == pp and M == 0:
                return (X, M, Contract(pa, pb, 0, fX, fM, self, kind))
            if self.sP == pp and X == 0:
                self.__setIndex__(inU)
                return self.add(X, M, a, b, kind)
            if self.sP == pm and M == 0:
                self.__setIndex__(inL)
                return self.add(X, M, a, b, kind)
                
        R = M/X if X else math.inf
        Rp = self.__Rp__(X, M, pp)
        Rm = self.__Rm__(X, M, pm)
        rc = self.__rab__(self.sP, pa, pb) # r at current price
        rp = self.__rab__(pp, pa, pb) # r at p+
        rm = self.__rab__(pm, pa, pb) # r at p-
        
        lF = self.L/(1-self.F) # Helper variable
        
        # Case in which we have too little X and the border isn't crossed
        if R > rc and Rp <= rp and self.L>0:
            a = X + self.L/self.sP - lF/pb
            b = lF*(1+self.sP/pb)-self.L-pa*X+M/pb-self.L*pa/self.sP
            c = self.L*pa-lF*self.sP-M
            sp1 = (-b + math.sqrt(b*b-4*a*c))/(2*a)
            dM = lF*(sp1 - self.sP)/(1-self.F)
            (dX, dM2) = self.tradeM(dM) # Already takes in fee
            X += dX
            M -= dM + dM2
        # Case in which we have too much X and the border isn't crossed
        elif R < rc and Rm >= rm and self.L>0: # Derivation, see P3
            a = X + lF/self.sP - self.L/pb 
            b = -lF-(X+lF/self.sP)*pa+self.L+(M+self.L*self.sP)/pb # Different than source (derived)
            c = lF*pa - self.L*self.sP-M # Different than source (derived)
            sp1 = (-b + math.sqrt(b*b-4*a*c))/(2*a)
            dX = lF*(1/sp1 - 1/self.sP)
            (dX2, dM) = self.trade(dX)
            M += dM
            X -= dX + dX2
        # Case in which we have too little X and the border is crossed
        elif R > rc and Rp >= rp:
            dM = self.__getdMpN__(pp)/(1-self.F)
            (dX, dM2) = self.tradeM(dM)
            X += dX
            M -= dM + dM2 
            self.sP = pp
            self.__setIndex__(inU)
            return self.add(X, M, a, b, kind)
        # Case in which we have too much X and the border is crossed
        elif R < rc and Rm < rm:
            dX = self.__getdXpN__(pm)/(1-self.F)
            (dX2, dM) = self.trade(dX)
            M += dM
            X = X - dX + dX2
            self.sP = pm
            self.__setIndex__(inL)
            return self.add(X, M, a, b, kind)
        
        # Compute the liquidity and left over assets/money
        Lx = round(X*self.sP*pb/(pb-self.sP))
        Lm = round(M/(self.sP-pa))
        if Lx == Lm:
            L = Lx
            self.__add__(pa, pb, Lx)
            return (0, 0, Contract(pa, pb, Lx, fX, fM, self, kind))
        elif Lx < Lm:
            M -= Lx*(self.sP-pa)
            self.__add__(pa, pb, Lx)
            return (0, M, Contract(pa, pb, Lx, fX, fM, self, kind))
        else:
            X -= Lm*(pb-self.sP)/self.sP/pb
            self.__add__(pa, pb, Lm)
            return (X, 0, Contract(pa, pb, Lm, fX, fM, self, kind))
    
    # Function that updates all references to new liquidity position
    def __add__(self, pa, pb, L):  
        if L <=0:
            return
        Il = self.index
        (Iu, _) = self.__getTickWindowUp__()
        pa = round(pa*pa)
        a = pa - self.s.minPriceRange
        pb = round(pb*pb)
        b = pb  - self.s.minPriceRange
        self.nR[a] += 1
        self.nR[b] += 1
        self.dL[a] = round(self.dL[a] + L)
        self.dL[b] = round(self.dL[b] - L)
        if Il < a and pa <= self.spot() and  pb >= self.spot():
            self.__setIndex__(a)
        elif Il < b and pb <= self.spot():
            self.__setIndex__(b)
        else:
            self.__setL__()
    
    def __rab__(self, p, pa, pb): 
        ph = min(max(p, pa), pb)
        return (ph-pa)/(1/ph - 1/pb) if not ph==pb else math.inf
    
    def __Rp__(self, X, M, pp): 
        return (M + max(self.L/(1-self.F)*(self.sP-pp), -M))/(X + self.L*(1/self.sP - 1/pp))
    
    def __Rm__(self, X, M, pm):
        d = max(self.L/(1-self.F)*(1/self.sP - 1/pm), -X)
        return (M + self.L*(self.sP-pm))/(X+d) if not d==-X else math.inf
    
    def retrieve(self, c):
        if c.L <= 0:
            return (0, 0)
        (fXn, fMn) = self.__retrieveFeesPerUnitInRange__(c.pa, c.pb)
        a = c.a-self.s.minPriceRange
        b = c.b-self.s.minPriceRange
        self.nR[a] -= 1
        self.nR[b] -= 1
        self.dL[a] = round(self.dL[a] - c.L)
        self.dL[b] = round(self.dL[b] + c.L)
        if self.nR[a] == 0 and not (self.dL[a] == 0):
            self.dL[a] = 0
        if self.nR[b] == 0 and not (self.dL[b] == 0) == 0:
            self.dL[b] = 0
        
        if self.sP < c.pa:
            X = c.L*(c.pb-c.pa)/(c.pb*c.pa)
            M = 0
        elif self.sP > c.pb:
            X = 0
            M = c.L*(c.pb-c.pa)
        else:
            X = c.L*(c.pb-self.sP)/(c.pb*self.sP)
            M = c.L*(self.sP-c.pa)
            self.__setL__()
        X += (fXn-c.fX)*c.L
        M += (fMn-c.fM)*c.L
        return (X, M)
    
    def __retrieveFeesPerUnitInRange__(self, pa, pb):
        a = math.floor(pa*pa)
        b = min(round(pb*pb), len(self.nR)-1+self.s.minPriceRange)
        X = 0
        M = 0
        while b>a:
            X += self.feeGrowthX[b - self.s.minPriceRange]
            M += self.feeGrowthM[b - self.s.minPriceRange]
            b -=1
        return (X, M)

class Contract(object):
    def __init__(self, pa, pb, L, fX, fM, amm, kind):
        self.name = round(random.random()*10000)
        self.pa = pa
        self.a = round(pa*pa)
        self.pb = pb
        self.b = round(pb*pb)
        self.L = L
        self.fX = fX
        self.fM = fM
        self.amm = amm
        self.kind = kind
        
    def expM(self):
        (X, M) = self.__expRet__()
        return M + self.amm.exp(X)
    
    def expX(self):
        (X, M) = self.__expRet__()
        return X + M/self.amm.spot()
    
    def __expRet__(self):
        (fXn, fMn) = self.amm.__retrieveFeesPerUnitInRange__(self.pa, self.pb)
        
        if self.amm.sP < self.pa:
            X = self.L*(self.pb-self.pa)/(self.pb*self.pa)
            M = 0
        elif self.amm.sP > self.pb:
            X = 0
            M = self.L*(self.pb-self.pa)
        else:
            X = self.L*(self.pb-self.amm.sP)/(self.pb*self.amm.sP)
            M = self.L*(self.amm.sP-self.pa)
        X += (fXn-self.fX)*self.L
        M += (fMn-self.fM)*self.L
        return (X, M)    

    def retrieve(self, amount = 0, name = -1):
        (asset, money) = self.amm.retrieve(self)
        if self.kind == "Assets":
            filled = money == 0
            (asset2, money) = self.amm.tradeM(money+amount)
            asset +=asset2
        elif self.kind == "Money":
            filled = asset == 0
            (asset, money2) = self.amm.trade(asset+amount)
            money += money2
        else:
            filled = math.nan
        self.L = 0
        return (asset, money, filled)
    
class BuyRangeOrder(Contract):
    def __init__(self, pa, pb, L, fX, fM, amm):
        self.__super__(pa, pb, L, fX, fM, amm)
    
    def retrieve(self, amount = 0, name = -1):
        # We want to retrieve assets. Extra amount is thus money
        (asset, money) = self.amm.retrieve(self)
        filled = money == 0
        (asset2, money) = self.amm.tradeM(money+amount)
        asset +=asset2
        return (asset, money, filled) 
    
class SellRangeOrder(Contract):
    def __init__(self, pa, pb, L, fX, fM, amm):
        self.__super__(pa, pb, L, fX, fM, amm)
    
    def retrieve(self, amount = 0, name = -1):
        # We want to retrieve moeny. Extra amount is thus asset
        (asset, money) = self.amm.retrieve(self)
        filled = asset == 0
        (asset, money2) = self.amm.trade(asset+amount)
        money += money2
        return (asset, money, filled)
        
        
# # Buy comparison. Price goes to 1020. Stepping into the pool is better option
# amm = AMM(simpy.Environment(), 1, Settings())
# (X, M, contract) = amm.add(0, 1020, 999, 1020)
# for i in range(400):
#     amm.trade(0.5)
#     amm.trade(-0.75)
# #amm.sP = math.sqrt(1021)
# print(amm.sP*amm.sP)
# (X, M) = amm.retrieve(contract)
# print(1020/(amm.__tradeM__(M)+X))

# amm = AMM(simpy.Environment(), 1, Settings())
# X = amm.__tradeM__(1020)
# print(1020/X)

# # Sell comparison. Price goes to 1020. 
# amm = AMM(simpy.Environment(), 1, Settings())
# (X, M, contract) = amm.add(1, 0, 1018, 1019)
# for i in range(400):
#     amm.trade(0.5)
#     amm.trade(-0.75)
# plt.scatter(range(len(amm.tradePrices)), amm.tradePrices)
# plt.show()
# listEven = amm.tradePrices[::2]
# listOdd = amm.tradePrices[1::2]
# plt.scatter(range(len(listOdd)), [a-b for a,b in zip(listOdd, listEven)])
# plt.show()
# print(amm.sP*amm.sP)
# (X, M) = amm.retrieve(contract)
# (X, M2) = amm.trade(X)
# print((M+M2)/1)

# amm = AMM(simpy.Environment(), 1, Settings())
# amm.sP = math.sqrt(1020)
# print(amm.trade(1)[1]/1)

# # Sell comparison. Price goes to 1020. 
# numTrades = 10000
# amm = AMM(simpy.Environment(), 1, Settings())
# for i in range(numTrades):
#     tr = random.randint(-3, 3)
#     while tr == 0:
#         tr = random.randint(-3, 3)
#     amm.trade(tr)
# plt.scatter(range(numTrades), [i*i for i in amm.allPrices])
# plt.show()

