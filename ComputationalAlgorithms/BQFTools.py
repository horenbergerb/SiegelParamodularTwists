#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 19:06:48 2020

@author: captdishwasher
"""

#So i'm not totally sure i like the structure of a class for BQF and a class for SL2Z


def gcd(a, b):
    while(b):
        a, b = b, a%b
    return a

class BQF:
    def __init__(self, a: int, b: int, c: int):
        self.a = a
        self.b = b
        self.c = c

    def print_bqf(self):
        print("{}x^2 + {}xy + {}y^2".format(self.a, self.b, self.c))
        
    def evaluate(self, x: int, y: int):
        return self.a*x**2 + self.b*x*y + self.c*y**2

    def transform(self, matrix, in_place = False):
        if in_place:
            new_bqf = matrix.mult(self)
            self.a = new_bqf.a
            self.b = new_bqf.b
            self.c = new_bqf.c
            return
        return matrix.mult(self)
    
    def discriminant(self):
        return self.b**2-(4*self.a*self.c)

    #WARNING: NOT PROPER EQUIVALENCE
    def same_coefficients(self, bqf):
        if self.a == bqf.a and self.b == bqf.b and self.c == bqf.c:
            return True
        return False
    
    def is_primitive(self):
        return self.same_coefficients(self.get_primitive_form())
    
    def get_primitive_form(self):
        ab_gcd = gcd(self.a,self.b)
        bc_gcd = gcd(self.b,self.c)
        total_gcd = gcd(ab_gcd, bc_gcd)
        if total_gcd > 1:
            return BQF(self.a/total_gcd, self.b/total_gcd, self.c/total_gcd)
        else:
            return self
        
    #positive definite bqfs only represent either positive or negative integers
    def is_pos_def(self):
        return self.discriminant() < 0
        
    def is_reduced(self):
        if not self.is_primitive():
            return False
        
        if abs(self.b) <= self.a and self.a <= self.c:
            if not (abs(self.b) == self.a or self.a == self.c) or self.b >= 0:
                    return True

        return False

    def reduce(self):
        
        if not self.is_pos_def():
            raise Exception("BQF is not positive definite")
            return
        
        #first we swap a and b to ensure a < c
        if self.a > self.c:
            swapper = SL2Z(0,-1,1,0)
            self.transform(swapper, in_place=True)
        
        #next we make abs(b) as small as possible
        #we know we can make b->(2am+b) for any m
        #using transformation
        #|1 m|
        #|0 1|
        
        #determine whether to add or subtract from b
        m = -1*(self.a/abs(self.a)*self.b/abs(self.b))
        subtractor = SL2Z(1, m, 0, 1)
        
        #subtract off of b until the next subtraction would make the magnitude larger
        while abs(self.b) > abs(self.b+2*m*self.a):
            self.transform(subtractor, in_place=True)
            self.print_bqf()
            
        #finally, edge cases b<0 and a=-b or a=c are not reduced yet
        if self.a == -1*self.b:
            b_sign_flip = SL2Z(1, 1, 0, 1)
            self.transform(b_sign_flip, in_place=True)
        if self.a == self.c:
            b_sign_flip = SL2Z(0, -1, 1, 0)
            self.transform(b_sign_flip, in_place=True)
            
                
class SL2Z:
    # (alpha  beta)
    # (rho    delta)
    def __init__(self, alpha: int, beta: int, rho: int, delta: int):
        self.alpha = alpha
        self.beta = beta
        self.rho = rho
        self.delta = delta
        
        if alpha*delta - beta*rho != 1:
            raise Exception("SL2Z matrix does not have det 1")
            
    def mult(self, bqf):
        new_a = bqf.evaluate(self.alpha, self.rho)
        new_b = (2*bqf.a*self.alpha*self.beta)+(bqf.b*self.alpha*self.delta)+(bqf.b*self.rho*self.beta)+(2*bqf.c*self.rho*self.delta)
        new_c = bqf.evaluate(self.beta, self.rho)
        return BQF(new_a, new_b, new_c)
    
    def print_mat(self):
        print("|{}\t {}|\n|{}\t {}|".format(self.alpha, self.beta, self.rho, self.delta))
            
        
test = BQF(2,7,23)
test.print_bqf()
print(test.evaluate(3, 2))


mat = SL2Z(2, 1, 1, 1)
mat.print_mat()
test.transform(mat, in_place = True)
test.print_bqf()
print(test.is_reduced())
test.reduce()
test.print_bqf()
print(test.is_reduced())