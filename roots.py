# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 21:11:26 2022

@author: lcuev
"""
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np

sin = np.sin
cos = np.cos

def newtons(guess,f,df):
    # edit range of runs to change n/o iterations of newtons method
    runs = range(10)
    for n,run in enumerate(runs):
        f_prime = df(guess)
        guess = guess + -1 * f(guess) * f_prime.invize()
        if f(guess).norm < 0.01:
            break
        
    return guess,n


def norm3(i,j,k):
    return i * i + j * j + k * k

def norm4(r,i,j,k):
    return r * r + i * i + j * j + k * k

# use this function to display quaternions in output console
def show(qat):
    print(qat.r, '+', qat.i, 'i +' , qat.j, 'j +', qat.k, 'k')

# quaternion class
class quat:
    def __init__(self,r,i,j,k):
        self.r = r
        self.i = i
        self.j = j
        self.k = k
        self.norm = norm4(r,i,j,k)
    
    def invize(self):
        ret = quat(self.r/self.norm,-self.i/self.norm,-self.j/self.norm,-self.k/self.norm)
        return ret
        
    def __mul__(self,opp):
        r = 0; i = 0; j = 0; k = 0;
        if type(opp) == quat:
            r = opp.r * self.r - opp.i * self.i - opp.j * self.j - opp.k * self.k
            i = opp.r * self.i + opp.i * self.r - opp.j * self.k + opp.k * self.j
            j = opp.r * self.j + opp.i * self.k + opp.j * self.r - opp.k * self.i
            k = opp.r * self.k - opp.i * self.j + opp.j * self.i + opp.k * self.r
        else:
            r = opp * self.r
            i = opp * self.i
            j = opp * self.j
            k = opp * self.k
        ret = quat(r,i,j,k)
        return ret
    def __add__(self,opp):
        r = 0; i = 0; j = 0; k = 0;
        if type(opp) == quat:
            r = self.r + opp.r
            i = self.i + opp.i 
            j = self.j + opp.j 
            k = self.k + opp.k 
        else:
            r = opp + self.r
            i = self.i
            j = self.j
            k = self.k
        ret = quat(r,i,j,k)
        return ret
    def __eq__(self,opp):
        ret = False
        error = 0.5
        if self.r > opp.r - error and self.r < opp.r + error\
        and self.i > opp.i - error and self.i < opp.i + error\
        and self.j > opp.j - error and self.j < opp.j + error\
        and self.k > opp.k - error and self.k < opp.k + error:
            ret = True
        return ret
    __rmul__ = __mul__
    __radd__ = __add__

# order 'self.N' polynomial class created from roots
class poly:
    
    def __init__(self,roots):
        self.roots = roots
        self.N = len(roots)
        self.psums = self.psums()
        self.coeffs = self.coeffs()
          
    def psums(self):
        psums = []
        for num in range(0,self.N+1):
            summ = 0
            for root in self.roots:
                fact = 1
                for k in range(num):
                    fact *= root
                summ += fact
            psums += [summ]
        return psums
    
    def coeffs(self):
        coeffs = [1]
        for k in range(1,self.N+1):
            summ = 0
            for i in range(1,k+1):
                summ += (-1) ** (i-1) * coeffs[k - i] * self.psums[i]
            coeffs += [1 / k * summ]
        for k in range(self.N + 1):
            coeffs[k] *= (-1) ** k
        return coeffs
        
    def f(self,x):
        ret = 1
        for r in self.roots:
            ret *= (x + -1 * r)
        return ret
    
    def df(self,x):
        ret = 0
        for k,coeff in enumerate(self.coeffs):
            fact = 1
            power = len(self.coeffs) - k - 1
            for i in range(power - 1):
                fact *= x
            ret += (power) * coeff * fact
        
        return ret
    
    def closest(self,guess):
        guess_mats = [guess.r,guess.i,guess.j,guess.k]
        root_mats = [[0,0,0,0] for i in range(self.N)]
        for i,root in enumerate(roots):
            root_mats[i] = [root.r,root.i,root.j,root.k]
        top = 0
        for l,coor in enumerate(root_mats[0]):
            top += (coor - guess_mats[l])**2 
        ret = 0
        for i,root in enumerate(root_mats):
            temp = 0
            for l,coor in enumerate(root):
                temp += (coor - guess_mats[l])**2 
            if temp < top:
                top = temp 
                ret = i
        return ret



# put in desired roots here (real or quaternion)
# quaternions take form quat(real,imaginary,jmaginary,kmaginary)
roots = [quat(1,0,0,0),quat(0,0.5,0,0),quat(0,-0.5,0,0), quat(-1,0,0.5,0), quat(-1,0,-0.5,0)]

# creates equation object with above roots
equation = poly(roots)

# gets functions represinting equation and derivitive of equation respectively
f,df = equation.f, equation.df

# scale of window in real imarinary plane
real_scale = 4
imag_scale = 4

# scale of output image
x_res = 200
y_res = 200

# image arrays
vals = [[0 for x in range(x_res)] for y in range(y_res)]
sets = [[0 for x in range(x_res)] for y in range(y_res)]

# render loop
for x in tqdm(range(x_res)):
    for y in range(y_res):
        # real and imaginary window created from x and y position in image
        r = (x / (x_res - 1) - 0.5) * real_scale
        i = (y / (y_res - 1) - 0.5) * imag_scale
        
        # initial guess fed from real and imaginary position
        guess = quat(r,i,0,0)

        # apply newtons method to guess
        guess,n = newtons(guess,f,df)
        
        # find root guess is closest to after application 
        sett = equation.closest(guess)
        
        # assign x and y values
        sets[x][y] = sett

# save image (file name can be edited to avoid replacement)
# example: fname = 'Anim/new_pic.png'
fname = 'Anim/frame' + str(np.random.random()%100) + '.png'
plt.imsave(fname,sets,cmap = 'Accent')


        


