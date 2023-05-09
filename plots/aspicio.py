#!/usr/bin/python
import numpy as np
import aspictools as ast

#latex dictonnary

modeldict = {'hi':r'$\mathrm{HI}$'}


paramdict = {'logEGeV':r'$\log\left(T_{\mathrm{reh}}/\mathrm{GeV}\right)$'}


class model(object):
    def __init__(self, name='',label=''):
        self.setname(name)
        self.label = label or self.setlabel(name)

    def setname(self, name):
        self.name = name

    def getname(self):
        return self.name

    def setlabel(self,name):
        return modeldict[name]

    def getlabel(self):
        if self.label:
            return self.label
        

class param(object):
    def __init__(self, name='', label='', value=0.0, index=None):
        self.setname(name)
        self.label = label or self.setlabel(name)
        self.value = value
        self.index = index

    def setname(self, name):
        self.name = name
   
    def getname(self):
        return self.name

    def setvalue(self, value):
        self.value = value

    def getvalue(self):
        return self.value

    def setlabel(self,name):
        return paramdict[name]

    def getlabel(self):
        if self.label:
            return self.label        


def read_modeldict(filename):
    global modeldict
    name, texform = np.loadtxt(filename,dtype='str',unpack=True,usecols=[0,1])

    for i in range(len(name)):
        modeldict[str(name[i])] = r'$'+str(texform[i])+'$'

    


def read_paramdict(filename):
    global paramdict
    name, texform = np.loadtxt(filename,dtype='str',unpack=True,usecols=[0,1])

    for i in range(len(name)):
        paramdict[str(name[i])] = r'$'+str(texform[i])+'$'



def read_aspicfile(filename):
    fixParams = []
    varParams = []
    paramNames = []

    f = open(filename,mode='rt')
    words = f.readline().split()

    while words[0] == '#':        
        if len(words) == 3:
            fixParams.append( param(name=words[1],value=np.float(words[2])) )
        else:
            for word in words[1:]:
                paramNames.append(word)

        words = f.readline().split()


    for line in f:
        values  = line.split()
        params = [param(name = paramNames[i], value = np.float(values[i])) for i in range(len(values))]
        varParams.append(params)

        
    f.close()
            
    return fixParams, varParams


    

def read_contour(filename):
    x, y =  np.loadtxt(filename,unpack=True,usecols=[0,1])
    return x,y


def read_nsr(filename):
    ns, r, N, lnRhoReh =  np.loadtxt(filename,unpack=True,usecols=[0,1,2,3])
    EGeV=ast.log_energy_reheat_ingev(lnRhoReh)
    return ns,r,N,EGeV


def aspicfile_to_nsr(filename):
    fixParams, varParams = read_aspicfile(filename)

    ns = []
    r = []
    N = []
    lnRhoReh = []
    
    for par in varParams:
        ns.append(par[0].getvalue())
        r.append(par[1].getvalue())
        N.append(par[2].getvalue())
        lnRhoReh.append(par[3].getvalue())

    ns = np.array(ns)
    r = np.array(r)
    N = np.array(N)
    lnRhoReh = np.array(lnRhoReh)
            
    EGeV = ast.log_energy_reheat_ingev(lnRhoReh)

    return ns,r,N,EGeV


def read_sr21(filename):
    eps2, eps1, N, lnRhoReh =  np.loadtxt(filename,unpack=True,usecols=[0,1,2,3])
    EGeV=ast.log_energy_reheat_ingev(lnRhoReh)
    return eps2,eps1,N,EGeV
