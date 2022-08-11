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
        
        # Compute liquidity of pool, given initial buy liquidity and price
        pb = math.sqrt(s.trueP+s.initialBidAsk)
        a = (1-s.trueP/s.trueP+s.initialBidAsk)
        b = -2*s.minLiquidity*s.trueP/pb
        c = -s.trueP*s.minLiquidity**2
        self.L = round((-b+math.sqrt(b**2-4*a*c))/(2*a))
        
        # Compute lower bound of initial range order given initial sell and pool liquidity
        a = round( ((s.trueP*(s.minLiquidity+self.L/pb) - s.minLiquidity*s.trueP)/self.L)**2 )
        
        # Make sure initial liquidity slowly disappears from market.
        env.process(initRangeOrder(Contract(math.sqrt(a), pb, self.L, 0, 0, self, None), s).start())
        
        # Update price for rounding errors
        self.sP = math.sqrt( (s.minLiquidity*s.trueP + self.L*math.sqrt(a))/(s.minLiquidity + self.L/pb) )
        
        
        # Variables needed for updates
        self.dL = [0] * (s.maxPriceRange - s.minPriceRange+1)
        self.nR = [0] * (s.maxPriceRange - s.minPriceRange+1)
        
        # Variables needed to keep track of fees
        self.feeGrowthX = [0] * (s.maxPriceRange - s.minPriceRange+1)
        self.feeGrowthM = [0] * (s.maxPriceRange - s.minPriceRange+1)
        self.F = s.fee 
        
        # Initialize first liquidity
        self.dL[a-s.minPriceRange] = self.L
        self.dL[s.trueP+s.initialBidAsk-s.minPriceRange] = -self.L
        self.nR[a-s.minPriceRange] += 1 # Number of references to price tick
        self.nR[s.trueP+s.initialBidAsk-s.minPriceRange] += 1 # Number of references to price tick
        self.index = a-s.minPriceRange
        
        # Additional statistics calculation needed for fee calculations
        self.prices = deque(maxlen=s.lookNTransactionsBack)
        self.prices.append(self.sP)
        self.tradePrices = []
    
    def spot(self):
        return self.sP*self.sP
    
    def bestBidPrice(self):
        return self.exp(1)

    def bestAskPrice(self):
        return self.exp(-1)
    
    def __setIndex__(self, index):
        self.index = index
        self.__setL__()
    
    def __setL__(self):
        L1 = round(np.sum(self.dL[:self.index+1]))
        L2 = -round(np.sum(self.dL[self.index+1:]))
        self.L = min(L1, L2)
    
    def __bigSpread__(self):
        X = 1000
        MBid = self.exp(X)
        MAsk = self.exp(-X)
        bid = max(MBid/X, self.s.minPriceRange)
        ask = min(MAsk/X, self.s.maxPriceRange)
        return ask-bid
    
    def __quoteSize__(self):
        pL = self.spot()-1
        pU = self.spot()+1
        return (self.getXPrice(pU), self.getXPrice(pL))
    
    
    # Trade X for M, as seen from the AMM. This means:
    # A negative amount of X means that we buy it from the AMM
    # A positive amount of X means that we sell it to the AMM
    # A negative amount of M returned means that we pay it to the AMM
    # A positive amount of M returned means that we get it from the AMM
    def trade(self, X, record=True):
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
            elif self.sP == pU and inU == len(self.nR):
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
                    if record:
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
                    if dM > 0 and record:
                        self.feeGrowthM[inU] += abs(dM*self.F/self.L)
                        self.marketBuyTransactions.push(dX)
                    self.__setIndex__(inU)
            (inL, inU, pL, pU) = self.__getTickWindow__()
        
        # If the trade crosses it's lower price (only possible when selling),
        # split up the trade to go to the crossing
        while X > 0 and (self.index > 0 or self.L>0):
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
                    if record:
                        self.feeGrowthX[inU] += abs(X*self.F/self.L)
                        self.marketSellTransactions.push(X)
                    X=0
                else:
                    M -= self.__getdMpN__(pL)
                    dX = self.__getdXpN__(pL)/(1-self.F)
                    X -= dX
                    self.sP = pL
                    if inU < len(self.nR) and not dX == 0 and record:
                        self.feeGrowthX[inU] += abs(dX*self.F/self.L)
                        self.marketSellTransactions.push(dX)
                    self.__setIndex__(inL)
            (inL, inU, pL, pU) = self.__getTickWindow__()
        
        if record:
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
    def tradeM(self, M, record=True):
        orM = M
        if record:
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
                    if record:
                        self.feeGrowthX[inU] += abs(dX*self.F/self.L)
                        self.marketSellTransactions.push(dX)
                else: 
                    # Update all returns of user for partial swap
                    dX = self.__getdXpN__(pL)/(1-self.F)
                    X -= dX
                    M -= self.__getdMpN__(pL)
                    
                    # Update liquidity
                    self.sP = pL
                    if dX > 0 and record:
                        self.feeGrowthX[inU] += abs(dX*self.F/self.L)
                        self.marketSellTransactions.push(dX)
                    self.__setIndex__(inL)
            (inL, inU, pL, pU) = self.__getTickWindow__()
        
        # If the trade crosses it's upper price (only possible when selling),
        # split up the trade to go to the crossing
        while M > 0 and (inU < len(self.nR) or self.L>0):
            if self.sP < pU and self.L == 0:
                break
            elif self.sP == pU and inU == len(self.nR):
                break
            elif self.sP == pU:
                self.__setIndex__(inU)
            else:
                pN = self.__getPricedM__(M*(1-self.F)) if inU < len(self.nR) else 0
                
                if pN <= pU:
                    dX = self.__getdXpN__(pN)
                    X -= dX
                    self.sP = pN
                    if record:
                        self.feeGrowthM[inU] += abs(M*self.F/self.L)
                        self.marketBuyTransactions.push(-dX)
                    M=0
                else:
                    dX = self.__getdXpN__(pU)
                    X -= dX
                    dM = self.__getdMpN__(pU)/(1-self.F)
                    M -= dM
                    self.sP = pU
                    if inU < len(self.nR) and not dM == 0 and record:
                        self.feeGrowthM[inU] += abs(dM*self.F/self.L)
                        self.marketBuyTransactions.push(dX)
                    self.__setIndex__(inU)
            (inL, inU, pL, pU) = self.__getTickWindow__()
        
        if record:
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
    def exp(self, X, p=None): 
        if not p:
            p = self.sP
            L = self.L
            index = self.index
        else:
            (index, _) = self.__getTickWindowDown__(int(np.ceil(p**2)-self.s.minPriceRange))
            L = round(np.sum(self.dL[:index+1]))
        M = 0
        while X<0 and L>0:
            (indexU, pU) = self.__getTickWindowUp__(index)
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
                L+= self.dL[indexU]
                index = indexU
        X = X*(1-self.F)
        while X > 0 and L>0:
            (indexL, pL) = self.__getTickWindowDown__(index) 
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
        
    def expM(self, M, p=None):
        if not p:
            p = self.sP
            L = self.L
            index = self.index
        else:
            (index, _) = self.__getTickWindowDown__(int(np.ceil(p**2)-self.s.minPriceRange))
            L = round(np.sum(self.dL[:index+1]))
            
        X = 0
        M = M*(1-self.F)
        (inU, pU) = self.__getTickWindowUp__(index)
        while M > 0 and (inU < len(self.nR) or L>0):
            if L == 0 and index == inU:
                L += self.dL[inU]
                index = inU
                break
            elif L == 0:
                M = 0
                break
            pN = M/L + p
            if pN < pU and pN > 0:
                X -= (1/pN - 1/p)*L
                M = 0
            else:
                X -= (1/pU - 1/p)*L
                M -= L*(pU-p)
                p = pU
                L += self.dL[inU]
                index = inU
                (inU, pU) = self.__getTickWindowUp__(index)
        
        return X
    
    def getPrice(self, assets, index = None, sP = None, L = None):
        if assets == 0:
            return self.sP
        
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
                pN = 1/(1/sP + assets*(1-self.F)/L)
                
                if pN >= pL:
                    return pN
                else:
                    assets -= (1/pL - 1/sP)*L
                    sP = pL
                    L = round(L-self.dL[index])
                    index = inL
            (inL, pL) = self.__getTickWindowDown__(index)
            
        if assets > 0:
            return math.sqrt(self.s.minPriceRange)
        
        while assets < 0:
            if sP < pU and L == 0:
                return sP
            elif sP == pU and inU == len(self.nR):
                return sP
            elif sP == pU:
                L += self.dL[pU]
                index = pU
            else:
                if (1/pU -1/sP)*L > assets:
                    pN = pU+1
                else:
                    pN = 1/(1/sP + assets/L)
                if pN <= pU:
                    return pN
                else: 
                    assets -= (1/pU - 1/sP)*L
                    sP = pU
                    L += self.dL[inU]
                    index = inU
            (inU, pU) = self.__getTickWindowUp__(index)
        
        return math.sqrt(self.s.maxPriceRange)
    
    def getExpLiq(self, assets, money, lower, upper, trading):
        index = self.index
        sP = self.sP
        (_, _, C) = self.add(assets, money, lower, upper, record = False)
        L = C.L
        fX = C.fX
        fM = C.fM
        expP = self.getPrice(trading)
        C.retrieve()
        self.sP = sP
        self.__setIndex__(index)
        return (L, fX, fM, expP)
    
    def __getTickWindow__(self):
        (inU, pU) = self.__getTickWindowUp__()
        (inL, pL) = self.__getTickWindowDown__()
        return (inL, inU, pL, pU)
    
    # Returns pair: index upperbound, sqrt price upper bound
    # When at upper bound, return upper bound+1
    def __getTickWindowUp__(self, inU=None):
        if not inU:
            inU = self.index
        
        inU +=1
        while inU < len(self.nR) and self.nR[inU] == 0:
            inU += 1
        pu = math.sqrt(inU+self.s.minPriceRange)
        return (inU, pu)
    
    # Returns pair: index next lowerbound, sqrt price lower bound
    # Returns -1 when there is no next lower bound
    def __getTickWindowDown__(self, index = None):
        if index == None:
            index = self.index
        
        inL = index
        
        if inL == 0:
            return (0, math.sqrt(self.s.minPriceRange))
        
        inL -=1
        while inL > 0 and self.nR[inL] == 0:
            inL -= 1
        
        if index == len(self.nR):
            return (inL, math.inf)
        else:
            return (inL, math.sqrt(index+self.s.minPriceRange))
    
    def __getPricedX__(self, dX): return 1/(1/self.sP + dX/self.L)
    def __getPricedM__(self, dM): return dM/self.L + self.sP
    def __getdMpN__(self, pN): return self.L*(pN-self.sP)
    def __getdXpN__(self, pN): return (1/pN - 1/self.sP)*self.L
    
    # Return X, M, L
    def add(self, X, M, a, b, kind=None, record=True): 
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
                return self.add(X, M, a, b, kind, record)
            if self.sP == pm and M == 0:
                self.__setIndex__(inL)
                return self.add(X, M, a, b, kind, record)
                
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
            (dX, dM2) = self.tradeM(dM, record) # Already takes in fee
            X += dX
            M -= dM + dM2
        # Case in which we have too much X and the border isn't crossed
        elif R < rc and Rm >= rm and self.L>0: # Derivation, see P3
            a = X + lF/self.sP - self.L/pb 
            b = -lF-(X+lF/self.sP)*pa+self.L+(M+self.L*self.sP)/pb # Different than source (derived)
            c = lF*pa - self.L*self.sP-M # Different than source (derived)
            sp1 = (-b + math.sqrt(b*b-4*a*c))/(2*a)
            dX = lF*(1/sp1 - 1/self.sP)
            (dX2, dM) = self.trade(dX, record)
            M += dM
            X -= dX + dX2
        # Case in which we have too little X and the border is crossed
        elif R > rc and Rp >= rp:
            dM = self.__getdMpN__(pp)/(1-self.F)
            (dX, dM2) = self.tradeM(dM, record)
            X += dX
            M -= dM + dM2 
            self.sP = pp
            self.__setIndex__(inU)
            return self.add(X, M, a, b, kind, record)
        # Case in which we have too much X and the border is crossed
        elif R < rc and Rm < rm:
            dX = self.__getdXpN__(pm)/(1-self.F)
            (dX2, dM) = self.trade(dX, record)
            M += dM
            X = X - dX + dX2
            self.sP = pm
            self.__setIndex__(inL)
            return self.add(X, M, a, b, kind, record)
        
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
    
    def getLiquidityDown(self):
        index = self.index
        sP = self.sP
        L = self.L
        lower = math.sqrt(self.s.minPriceRange)
        (inL, pL) = self.__getTickWindowDown__(index)
        
        X = 0
        
        while sP > lower and (L > 0 or sP == pL) and index>=0:
            X += (1/pL - 1/sP)*L
            L = round(L-self.dL[index])
            index = inL
            sP = pL
            (inL, pL) = self.__getTickWindowDown__(index)
            
        return X
    
    def getLiquidityUp(self):
        index = self.index
        sP = self.sP
        L = self.L
        upper = math.sqrt(self.s.maxPriceRange)
        (inU, pU) = self.__getTickWindowUp__(index)
        
        X = 0
        while sP < upper and (L > 0 or sP == pU):
            X -= (1/pU - 1/sP)*L
            index = inU
            L = round(L+self.dL[index])
            sP = pU
            (inU, pU) = self.__getTickWindowUp__(index)
        return X
    
        
    # Computes the point from there is no liquidity 
    # between current price and this price
    # Returns price when there is liquidity over whole range
    def getLiquidityStart(self, p2):
        dLCum = np.cumsum(self.dL)
        spot = self.spot()
        if p2 > spot:
            lower = int(np.ceil(self.spot()))-self.s.minPriceRange
            upper = int(np.floor(p2))-self.s.minPriceRange
        else:
            lower = int(np.ceil(p2))-self.s.minPriceRange
            upper = int(np.floor(self.spot()))-self.s.minPriceRange
            
        try:
            index = np.where(dLCum[lower:upper+1]==0)[0][0]+lower
        except:
            if p2 > spot:
                index = int(np.floor(p2))-self.s.minPriceRange
            else:
                index = int(np.ceil(p2))-self.s.minPriceRange
        return index+self.s.minPriceRange
    
    def getPriceX(self, X):
        (sP, index, L) = (self.sP, self.index, self.L)
        
        # While there are assets left that are sold
        while (X > 0):
            (inL, pL) = self.__getTickWindowDown__(index)
            if ((inL <= 0 or pL<sP) and L==0):
                return sP**2
            elif round(pL-sP,14) == 0:
                L -= self.dL[index]
                index = inL
                continue
            dP = min((1/pL - 1/sP), X/L)
            X = round(X-dP*L,8)
            sP = 1/(1/sP+dP)
            L -= self.dL[index]
            index = inL
        # While there are assets left that are bought
        while (X < 0):
            (inU, pU) = self.__getTickWindowUp__(index)
            if (inU >= len(self.nR) or (L == 0 and sP < pU)):
                return sP**2
            elif round(sP-pU, 14)==0:
                L =+ self.dL[inU]
                index = inU
                continue
            dP = max((1/pU - 1/sP), X/L)
            X = round(X-dP*L,8)
            sP = 1/(1/sP+dP)
            index = inU
            L += self.dL[inU]
        return sP**2
    
    # Function that computes how much X has to be added to get to the price.
    # We assume price <= spot when called by the trader
    def getXPrice(self, price):
        gP = math.sqrt(price)
        (sP, index, L) = (self.sP, self.index, self.L)
        X = 0
        
        while gP < sP:
            (inL, pL) = self.__getTickWindowDown__(index)
            if ((inL <= 0 or pL>sP) and L==0):
                return sP**2
            elif round(pL-sP, 14) == 0:
                L -= self.dL[index]
                index = inL
                continue
            border = max(pL, gP)
            dP = 1/border - 1/sP
            X = round(X+dP*L,8)
            sP = sP + dP if border==pL else gP
            L == L - self.dL[index] if border == pL else L
            index = inL if border == pL else index
        
        
        while gP > sP:
            (inU, pU) = self.__getTickWindowUp__(index)
            if ((inU >= len(self.nR)) or (sP<pU and L==0)):
                return X
            elif (L==0 and round(sP-pU, 14)==0):
                index = inU
                L += self.dL[index]
                continue
            border = min(pU, gP)
            dP = border - sP
            X -= (1/(sP+dP)-1/sP)*L
            sP = sP+dP if border == pU else gP
            index = inU if border == pU else index
            L = L + self.dL[inU] if border == pU else L
        return X
    
    # Compute amount of X one gets for positive amount of M
    def getX(self, M, price=None):
        if not price:  
            (sP, index, L) = (self.sP, self.index, self.L)
        else:
            sP = math.sqrt(price)
            index = math.floor(price)-self.s.minPriceRange
            L = round(np.sum(self.dL[:index+1]))
        X = 0
        while M>0:
            (inU, pU) = self.__getTickWindowUp__(index)
            if ((inU >= len(self.nR)) or (sP<pU and L==0)):
                return X
            elif (L==0 and round(sP-pU, 14)==0):
                index = inU
                L += self.dL[index]
                continue
            dP = min(pU - sP, M/L)
            X -= (1/(sP+dP)-1/sP)*L
            M = M-dP*L if dP == pU - sP else 0
            sP = sP+dP
            index = inU
            L += self.dL[inU]
        return X
    
    # Get the price after trading M
    def getPriceM(self, M):
        (sP, index, L) = (self.sP, self.index, self.L)
        
        # While there is M left that is sold (X is bought)
        while (M > 0):
            (inU, pU) = self.__getTickWindowUp__(index)
            if round(sP-pU, 14)==0:
                L += self.dL[inU]
                index = inU
                continue
            elif ((sP < pU and L==0) or inU == len(self.nR)):
                return sP**2
            
            dP = min(M/L, pU-sP)
            M -= dP*L
            sP += dP
            if (sP == pU):
                index = inU
                L += self.dL[inU]
        # While there is M left that is bought (X is sold)
        while (M < 0):
            (inL, pL) = self.__getTickWindowDown__(index)
            if round(sP-pL, 14)==0:
                L -= self.dL[index]
                index = inL
                continue
            elif ((inL <= 0 or pL==sP) and L==0):
                return sP**2
            dP = max(M/L, pL-sP)
            M -= dP*L
            sP += dP
            if (sP == pL):
                L -= self.dL[index]
                index = inL
                
        return sP**2
    
    # Get the amount of M that has to be added to get to the price
    # We assume price >= spot
    def getMPrice(self, price):
        gP = math.sqrt(price)
        (sP, index, L) = (self.sP, self.index, self.L)
        M = 0
        
        while gP > sP:
            (inU, pU) = self.__getTickWindowUp__(index)
            if round(sP-pU, 14)==0:
                L += self.dL[inU]
                index = inU
                continue
            elif ((sP < pU and L==0) or inU == len(self.nR)):
                return math.inf
            
            dP = min(gP-sP, pU-sP)
            M += dP*L
            sP += dP
            if (sP == pU):
                index = inU
                L += self.dL[inU]
        return M
    
    def getM(self, X, price=None):
        if not price:
            (sP, index, L) = (self.sP, self.index, self.L)
        else:
            sP = math.sqrt(price)
            index = math.floor(price)-self.s.minPriceRange
            L = round(np.sum(self.dL[:index+1]))
        M = 0
        while X>0:
            (inL, pL) = self.__getTickWindowDown__(index)
            if ((inL <= 0 or pL<sP) and L==0):
                return M
            elif (L==0 and round(sP-pL, 14)==0):
                L-= self.dL[index]
                index = inL
                continue
            dP = min(1/pL - 1/sP, X/L)
            X -= dP*L
            M -= (1/(dP + 1/sP) - sP)*L 
            sP = 1/(dP + 1/sP)
            L -= self.dL[index]
            index = inL
        return M
    
    def getFirstLiqChange(self, oP):
        if (oP < self.spot()):
            index = math.ceil(oP)-self.s.minPriceRange + 1
            while self.dL[index] == 0 and index<self.index:
                index += 1
        else:
            index = math.floor(oP)-self.s.minPriceRange - 1
            while self.dL[index] == 0 and index > self.index:
                index -= 1
        return index
    
    def getL(self, rL, rU):
        l = rL - self.s.minPriceRange
        u = rU - self.s.minPriceRange
        L = np.sum(self.dL[:l+1])
        tot = 0
        while l < u:
            tot += L
            l += 1
            L += self.dL[l]
        
        return tot

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
        
    def expM(self, assets = 0):
        (X, M) = self.__expRet__()
        return M + self.amm.exp(X+assets)
    
    def expX(self, money = 0):
        (X, M) = self.__expRet__()
        return X + self.amm.expM(M+money)
    
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
        (X, M) = self.amm.retrieve(self)
        if self.kind == "X":
            filled = M == 0
            (X2, M) = self.amm.tradeM(M+amount)
            X += X2
        elif self.kind == "M":
            filled = X == 0
            (X, M2) = self.amm.trade(X+amount)
            M += M2
        else:
            filled = math.nan
        self.L = 0
        return (X, M, filled)      
    
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

class initRangeOrder():
    def __init__(self, contract, settings):
        self.c = contract
        self.lower = self.c.a
        self.upper = self.c.b
        self.amm = self.c.amm
        self.s = settings
        self.end = settings.shockWait
        self.initL = self.c.L
    
    def start(self):
        while self.c.L > 0.01*self.initL:
            yield self.amm.env.timeout(100)
            with self.amm.capacity.request() as request:
                yield request
                L = self.c.L
                Lshould = self.initL*(self.end-self.amm.env.now)/self.end
                per = max(Lshould/L, 0.01)
                (assets, money, _) = self.c.retrieve()
                assets = max(assets*per,0)
                money = max(money*per,0)
                if assets > 0 and money > 0:
                    (_, _, self.c) = self.amm.add(assets, money, self.lower, self.upper)
                else: 
                    self.c.L = 0









