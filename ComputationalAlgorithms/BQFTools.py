#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 19:06:48 2020

@author: captdishwasher
"""

#TO DO:
#test the SPMF and BQF functions rigorously
#use LMFDB JSON file
#maybe crunch some actual coefficients

#finish spmf stuff
    #add level to SPMF so p doesn't divide N
    #be more rigorous about setting weight
#ensure all matrices are multiplying properly
#be better about edge cases. I throw exceptions at points that need to be fleshed out
#throw exceptions when bqfs or matrices have non-integer values

#BIG: talk to kirk about wiring this into his coefficient calculating code


#returns gcd of a, b
def gcd(a, b):
    while(b):
        a, b = b, a%b
    return a

#returns whether a divides b
def divides(a, b):
    return (b%a == 0)

def divides_exactly(a, b):
    return (divides(a,b) and not divides(a, b/a))

#This is NOT OPTIMIZED (not sure how much better it could be)
#returns 1 if a is quad residue mod p and a is not 0 mod p
#returns -1 if a is non-quadratic residue mod p
#returns 0 if a is 0 mod p
def legendre(a, p):
    #returns legendre (a/p)
    if (a%p == 0):
        return 0
    for x in range (1, p):
        #check if input is quad residue mod p
        if((x*x)%p == a%p):
            return 1
    return -1

#a siegel paramodular form is modeled here
#as a collection of bqfs and their associated coefficients
#the list of BQFs and coeffs is not sorted
#but there are seach functions for getting what you need out of them
class SPMF:
    def __init__(self, weight = None):
        #[[BQF1, coeff1], [BQF2, coeff2], ...]
        self.coeffs = []
        self.weight = weight
        self.level

    #made for LMFDB
    def load_from_json(self, json_file_name):
        with open(filename) as f:
        data = json.load(f)
        self.weight = int(data["weight"])
        
        for disc in data["Fourier_coefficients"]:
            for cur_BQF in data["Fourier_coefficients"][disc]:
                cleaned_BQF = BQF.strip('()').replace(" ", "")
                self.add_coeff(BQF(int(cleaned_BQF.split()[0]), int(cleaned_BQF.split()[1]), int(cleaned_BQF.split()[2])), int(data["Fourier_coefficients"][disc][cur_BQF]))
    
    #should have a better sol'n than raising error on duplicates            
    def add_coeff(self, BQF, coeff):
        if any(x[0].equals(BQF) for x in self.coeffs):
            raise Exception("Same BQF added to SPMF twice")
        self.coeffs.append([BQF, coeff])
    
    #search function by discriminant
    def get_by_disc(self, disc):
        return [x for x in self.coeffs if x.discriminant()==disc]
    
    #gets a list of discriminants
    def get_all_discs(self):
        return set(x[0].discriminant() for x in self.coeffs)
    
    #finds a particular bqf and coeff matching a given bqf
    def get_by_BQF(self, BQF):
        for x in coeffs:
            if x[0].same_coefficients(BQF):
                return x
        return None
    
    #parses all the bqfs and coeffs into an LMFDB-style dict. Ready to be dumped to JSON
    def get_coeff_dict(self):
        discs = self.get_all_discs()
        output = dict()
        for disc in discs:
            cur_coeffs = self.get_by_disc(disc)
            for coeff in cur_coeffs:
                output[disc]["({}, {}, {})".format(coeff[0].a, coeff[0].b, coeff[0].c)] = coeff[1]
        return output
        
    #twist described in case 1 of Dr. Jennifer Johnson-Leung's paper
    def twist_case_1(p, level):
        #need to implement level of SPMFs for this check:
        #if divides(p, self.N):
            #raise Exception("p cannot divide N")
        
        #note: dirichlet character is the legendre mod p
        
        #does new SPMF have the same weight? I think it does
        new_SPMF = SPMF(weight=self.weight)
        
        for cur_BQF in self.coeffs:
            #conditions for case 1 from jen's paper
            if divides(p, cur_BQF[0].b) or not divides(p**4, cur_BQF[0].a):
                continue
            #first term of the twist function
            term_1 = p**(1-self.weight)
            #second term
            term_2 = legendre(cur_BQF[0].b, p)
            
            #the sum terms from 1 to p-1
            #NOTE: we are going to handle the p**-1 separately since we need to be careful with division
            total_sum = 0
            for beta in range(1, p):
                sum_term_1 = legendre(beta, p)
                #notice we left out p**-1
                new_transform = Matrix(1, -1*beta, 0, p)
                #sum term 2 is the A^tSA process
                sum_term_2 = cur_BQF[0].transform(new_transform)
                #checking divis conditions (see 9.2.2 in Overview of Siegel Paramodular Forms)
                if not (divides(p**2, sum_term_2.b) and divides(p**4, sum_term_2.c)):
                    raise Exception("twist case 1 did not have sufficient divisibility wrt p")
                sum_term_2.b = sum_term_2.b/p**2
                sum_term_2.c = sum_term_2.c/p**4
                #now we just need the coeff of this new bQF
                sum_term_coeff = self.get_by_bqf(sum_term_2)
                if sum_term_coeff == None:
                    raise Exception("twist case 1 needs coefficients we don't have")
                #calculating the completed term for the sum
                total_sum += sum_term_1*sum_term_coeff[1]
            
            #putting the sum together with the rest of the twist
            new_coeff = term_1*term_2*total_sum
            
            #adding to our new SPMF
            new_SPMF.add_coeff(cur_BQF, new_coeff)
            
        return new_SPMF

#models a bqf and the things you'd want to do with it
class BQF:
    def __init__(self, a: int, b: int, c: int):
        self.a = a
        self.b = b
        self.c = c
        self.disc = self.b**2-(4*self.a*self.c)

    #aesthetic function
    def print_bqf(self):
        print("{}x^2 + {}xy + {}y^2".format(self.a, self.b, self.c))
        
    #plugs in two values and calculates output
    def evaluate(self, x: int, y: int):
        return self.a*x**2 + self.b*x*y + self.c*y**2

    #transforms by a matrix
    def transform(self, matrix, in_place = False):
        if in_place:
            new_bqf = matrix.mult(self)
            self.a = new_bqf.a
            self.b = new_bqf.b
            self.c = new_bqf.c
            self.disc = self.b**2-(4*self.a*self.c)
            return
        return matrix.mult(self)
    
    #shows you the discriminant
    def discriminant(self):
        return self.disc

    #WARNING: NOT PROPER EQUIVALENCE
    #compares two bqfs to see if a1==a2, b1==b2, c1==c2
    def same_coefficients(self, bqf):
        if self.a == bqf.a and self.b == bqf.b and self.c == bqf.c:
            return True
        return False
    
    #boolean check if primitive
    def is_primitive(self):
        return self.same_coefficients(self.get_primitive_form())
    
    #returns equivalent primitive form
    def get_primitive_form(self):
        ab_gcd = gcd(self.a,self.b)
        bc_gcd = gcd(self.b,self.c)
        total_gcd = gcd(ab_gcd, bc_gcd)
        if total_gcd > 1:
            return BQF(self.a/total_gcd, self.b/total_gcd, self.c/total_gcd)
        else:
            return self
        
    #positive definite bqfs only represent either positive or negative integers
    #this returns boolean. true if bqf is positive definite
    def is_pos_def(self):
        return self.discriminant() < 0
        
    #checks if the bqf is reduced
    def is_reduced(self):
        if not self.is_primitive():
            return False
        
        if abs(self.b) <= self.a and self.a <= self.c:
            if not (abs(self.b) == self.a or self.a == self.c) or self.b >= 0:
                    return True

        return False

    #actually reduces the bqf in place
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
        new_b = (2*bqf.a*self.alpha*self.beta)+(bqf.b*self.alpha*self.delta)+(bqf.b*self.rho*self.beta)+(2*bqf.c*self.beta*self.delta)
        new_c = bqf.evaluate(self.beta, self.delta)
        return BQF(new_a, new_b, new_c)
    
    def print_mat(self):
        print("|{}\t {}|\n|{}\t {}|".format(self.alpha, self.beta, self.rho, self.delta))
    
    def transpose(self, in_place=False):
        if in_place:
            self.rho, self.beta = self.beta, self.rho
            return
        else:
            return SL2Z(self.alpha, self.rho, self.beta, self.delta)
    
#general matrix class for use in twisting and other things
class Matrix:
    # (alpha  beta)
    # (rho    delta)
    def __init__(self, alpha: int, beta: int, rho: int, delta: int):
        self.alpha = alpha
        self.beta = beta
        self.rho = rho
        self.delta = delta
            
    def mult(self, bqf):
        new_a = bqf.evaluate(self.alpha, self.rho)
        new_b = (2*bqf.a*self.alpha*self.beta)+(bqf.b*self.alpha*self.delta)+(bqf.b*self.rho*self.beta)+(2*bqf.c*self.beta*self.delta)
        new_c = bqf.evaluate(self.beta, self.delta)
        return BQF(new_a, new_b, new_c)
    
    def print_mat(self):
        print("|{}\t {}|\n|{}\t {}|".format(self.alpha, self.beta, self.rho, self.delta))
    
    def transpose(self, in_place=False):
        if in_place:
            self.rho, self.beta = self.beta, self.rho
            return
        else:
            return SL2Z(self.alpha, self.rho, self.beta, self.delta)

    
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