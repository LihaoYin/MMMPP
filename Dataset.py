# encoding:utf-8

import numpy as np
import csv
import math
from scipy.stats import norm
import torch
import random
import copy

# use_cuda = torch.cuda.is_available()
# device = torch.device('cuda:0' if use_cuda else 'cpu')

def fpca(X, l, thres = 0.01, pmax = 20):
    values, vectors = np.linalg.eig(X)
    i = 0
    while i < min(len(values), 20) and values[i] > thres*values[0]:
        i += 1
    return values[0:i]*l, vectors[:,0:i]/math.sqrt(l)

def cov_from_eig(values, vectors):
    Q = vectors
    R = Q.T
    L = np.diag(values)
    return Q.dot(L).dot(R)

class Dataset:
    def __init__(self, dataset_path, label_path = '', model = 'single_level', naccounts = 0, mdays = 1, event_types = 1):
        self.labels = {}
        if label_path:
            self.naccounts = self.read_labels(label_path)
        else:
            self.naccounts = naccounts
        self.mdays = mdays
        self.event_types = event_types
        self.accounts = {}
        self.Process = self.read_data(dataset_path)
        self.paras = {}
        self.elems = {}
        if model == 'multi_level':
            self.process_multilevel()
        if model == 'single_level':
            self.process_singlelevel()

        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    def read_labels(self, label_path):
        field_types = [('accountNumber', int),
                       ('label', int)]
        naccounts = 0
        with open(label_path) as f:
            for row in csv.DictReader(f):
                row.update((key, conversion(row[key])) for key, conversion in field_types)
                self.labels[row['accountNumber']] = row['label']
                naccounts += 1
        return naccounts

    def read_data(self, dataset_path):
        field_types = [('accountNumber', int),
                       ('timeStamps', float),
                       ('days', int),
                       ('category', int)]
        process = np.empty(self.naccounts * self.mdays, dtype=object)
        process[...] = [[] for _ in range(self.naccounts * self.mdays)]
        process.shape = (self.naccounts, self.mdays)
        i = 0
        with open(dataset_path) as f:
            for row in csv.DictReader(f):
                row.update((key, conversion(row[key])) for key, conversion in field_types)
                if row['accountNumber'] in self.accounts:
                    id = self.accounts[row['accountNumber']]
                else:
                    id = i
                    self.accounts[row['accountNumber']] = id
                    i += 1
                process[id, row['days']-1].append(row['timeStamps'])



        return process

    def pkernel(self, q, kernel = "gaussian", mean = 0, sd = math.sqrt(0.2)):
        if kernel == "gaussian":
            return norm.cdf(q,mean,sd)

    def dkernel(self, t, bwd, kernel = "gaussian", mean = 0, sd = math.sqrt(0.2)):
        if kernel == "gaussian":
            return norm.pdf(t/bwd,mean,sd)/bwd


    def process_multilevel(self, lbd = 0, ubd = 1, bwd = 0.1, ngrids = 100):
        l = (ubd-lbd)/ngrids
        grids = np.array([x*l+lbd for x in range(ngrids+1)])
        grids_mid = np.array([x*l+l/2+lbd for x in range(ngrids)])

        A = np.zeros((ngrids,ngrids))
        B = np.zeros((ngrids,ngrids))
        C = np.zeros((ngrids,ngrids))
        A2 = np.zeros((ngrids,ngrids))

        edge = self.pkernel((grids_mid-lbd)/bwd) - self.pkernel((grids_mid-ubd)/bwd)
        edge2 = np.outer(edge, edge)

        self.elems['a'] = torch.empty(self.naccounts,ngrids,ngrids)
        self.elems['b'] = torch.empty(self.naccounts,ngrids)
        self.elems['v'] = torch.empty(self.naccounts,ngrids)
        self.elems['y'] = torch.empty(self.naccounts,ngrids)
        self.elems['ct'] = torch.empty(self.naccounts,ngrids)

        for i in range(self.naccounts):
            for j in range(self.mdays):
                process = self.Process[i,j]
                if process:
                    tmp = self.dkernel(np.subtract.outer(process, grids_mid), bwd)
                    tmp1 = np.sum(tmp, axis=0)
                    tmp2 = np.outer(tmp1, tmp1)
                    A2 += tmp2
                    A += tmp2 - np.matmul(tmp.T, tmp)

        for i in range(self.naccounts):
            process = list(np.concatenate(self.Process[i,]).flat)
            tmp = self.dkernel(np.subtract.outer(process, grids_mid), bwd)
            tmp1 = np.sum(tmp, axis=0)
            tmp2 = np.outer(tmp1, tmp1)
            if process:
                B += tmp2

            tmp2 -= np.matmul(tmp.T, tmp)
            self.elems['a'][i,:,:] = torch.tensor(tmp2/edge2)
            self.elems['b'][i,:] = torch.tensor(tmp1/edge/self.naccounts)
            ct = np.histogram(process, bins=grids)[0] + 1
            v = (ubd-lbd)/ngrids/ct
            self.elems['y'][i,:] = torch.tensor((ct-1)/v)
            self.elems['v'][i,:] = torch.tensor(v)
            self.elems['ct'][i,:] = torch.tensor(ct)

        B -= A2

        for j in range(self.mdays):
            process = list(np.concatenate(self.Process[:,j]).flat)
            if process:
                tmp = self.dkernel(np.subtract.outer(process, grids_mid), bwd)
                tmp1 = np.sum(tmp, axis=0)
                tmp2 = np.outer(tmp1, tmp1)
                C += tmp2

        C -= A2

        process = list(np.concatenate(list(np.concatenate(self.Process).flat)).flat)
        tmp = self.dkernel(np.subtract.outer(process, grids_mid), bwd)
        tmp1 = np.sum(tmp, axis=0)
        tmp2 = np.outer(tmp1, tmp1)

        D = tmp2 - B - C - A2

        A = A/self.naccounts/self.mdays/edge2
        B = B/self.naccounts/self.mdays/(self.mdays-1)/edge2
        C = C/self.naccounts/self.mdays/(self.naccounts-1)/edge2
        D = D/self.naccounts/self.mdays/(self.naccounts-1)/(self.mdays-1)/edge2

        self.paras['Cov.y'] = np.log(C)-np.log(D)
        self.paras['Cov.z'] = np.log(A)+np.log(D)-np.log(B)-np.log(C)

        values, vectors = fpca(self.paras['Cov.y'], l)
        self.paras['Cov.y'] = cov_from_eig(values, vectors)
        values, vectors = fpca(self.paras['Cov.z'], l)
        self.paras['Cov.z'] = cov_from_eig(values, vectors)


    def process_singlelevel(self, lbd = 0, ubd = 1, bwd = 0.1, ngrids = 100):
        l = (ubd-lbd)/ngrids
        grids = np.array([x*l+lbd for x in range(ngrids+1)])
        grids_mid = np.array([x*l+l/2+lbd for x in range(ngrids)])

        edge = self.pkernel((grids_mid-lbd)/bwd) - self.pkernel((grids_mid-ubd)/bwd)
        edge2 = np.outer(edge, edge)

        self.elems['a'] = torch.empty(self.naccounts,ngrids,ngrids)
        self.elems['b'] = torch.empty(self.naccounts,ngrids)
        self.elems['v'] = torch.empty(self.naccounts,ngrids)
        self.elems['y'] = torch.empty(self.naccounts,ngrids)
        self.elems['ct'] = torch.empty(self.naccounts,ngrids)
        for i in range(self.naccounts):
            process = list(np.concatenate(self.Process[0,]).flat)
            tmp = self.dkernel(np.subtract.outer(process, grids_mid), bwd)
            tmp1 = np.sum(tmp, axis=0)
            tmp2 = np.outer(tmp1, tmp1)
            tmp2 = tmp2 - np.matmul(tmp.T, tmp)
            self.elems['a'][i,:,:] = torch.tensor(tmp2/edge2/self.naccounts)
            self.elems['b'][i,:] = torch.tensor(tmp1/edge/self.naccounts)
            ct = np.histogram(process, bins=grids)[0] + 1
            v = (ubd-lbd)/ngrids/ct
            self.elems['y'][i,:] = torch.tensor((ct-1)/v)
            self.elems['v'][i,:] = torch.tensor(v)
            self.elems['ct'][i,:] = torch.tensor(ct)

    def convergence(self, new, log_likelihood, reltol, abstol):
        delta1 = 0
        selected_keys = ('mu.x', 'Cov.x')
        for key in selected_keys:
            delta1 = max(delta1, torch.max(torch.abs(new[key]-self.paras[key])))
        if len(log_likelihood) > 1:
            delta2 = abs((log_likelihood[-2]-log_likelihood[-1])/log_likelihood[-2])
        else:
            delta2 = 0

        return delta1 < abstol and delta2 < reltol

    def ES_initial(self, nclusters, ngrids = 100):
        self.paras['pi'] = np.full(nclusters, 1/nclusters)
        self.paras['mu.x'] = torch.empty(self.nclusters, ngrids)
        self.paras['Cov.x'] = torch.empty(self.nclusters, ngrids, ngrids)
        omega = random.choices(range(nclusters), k=self.naccounts)
        omega = np.array(omega)
        for c in range(nclusters):
            ids = np.where(omega == c)[0]
            self.paras['mu.x'][c,...] = torch.sum(self.elems['b'].index_select(0, torch.tensor(ids)), axis=0)
            self.paras['Cov.x'][c,...] = torch.sum(self.elems['a'].index_select(0, torch.tensor(ids)), axis=0)

            self.paras['Cov.x'][c,...] = torch.log(self.paras['Cov.x'][c,...]*self.paras['pi'][c]/ (torch.outer(self.paras['mu.x'][c,],self.paras['mu.x'][c,])) )
            self.paras['mu.x'][c,...] = torch.log(self.paras['mu.x'][c,...]/self.paras['pi'][c]) - torch.diag(self.paras['Cov.x'][c,...])/2


    def MSMPP(self, nclusters, nsampling = 2000, maxIter = 50, ngrids = 100, reltol = 10**-4, abstol = 10**-2):
        new = copy.deepcopy(self.paras)
        log_likelihood = []
        for iter in range(maxIter):
            X = torch.empty(nclusters, nsampling, ngrids)
            for c in range(nclusters):
                values, vectors = fpca(self.paras['Cov.x'][c,...], l)
                xi_x = torch.randn(len(values), nsampling) * values[:,None]
                X[c,...] = (torch.matmul(torch.tensor(vectors), xi_x) + self.paras['mu.x'][c,:][:,None]).T
            f = torch.matmul(X, (self.elems['v']*self.elems['y']).T) - torch.matmul(torch.exp(X), (self.elems['v']*self.elems['ct']).T)
            f_max = torch.amax(f,dim=(0,1))
            omega = torch.mean(torch.exp(f - f_max[None,None,:]),axis=1) * self.paras['pi'][:,None]
            log_likelihood.append( torch.sum(torch.log(torch.sum(omega,axis=0)) + f_max) )
            omega = omega - torch.log(torch.sum(omega,axis=0))[None,:]
            id = torch.argmax(omega, axis=0)
            omega = omega/torch.sum(omega,axis=0)[None,:]
            new['pi'] = torch.mean(omega, axis=1)
            new['mu.x'] = torch.matmul(omega.float(), self.elems['b'])
            new['Cov.x'] = torch.matmul(self.elems['a'].permute(1,2,0),omega.T.float()).permute(2,0,1)

            for c in range(nclusters):
                new['Cov.x'][c,...] = torch.log(new['Cov.x'][c,...]*new['pi'][c]/ (torch.outer(new['mu.x'][c,],new['mu.x'][c,])) )
                new['mu.x'][c,...] = torch.log(new['mu.x'][c,...]/new['pi'][c]) - torch.diag(new['Cov.x'][c,...])/2

            if self.convergence(new, log_likelihood, reltol, abstol):
                break

            self.paras = copy.deepcopy(new)

        return id, log_likelihood[-1]


    def MSMPP_GPU(self, nclusters, nsampling = 2000, maxIter = 50, ngrids = 100, reltol = 10**-4, abstol = 10**-2):
        new = copy.deepcopy(self.paras)
        log_likelihood = []
        for iter in range(maxIter):
            X = torch.empty(nclusters, nsampling, ngrids).to(self.device)
            for c in range(nclusters):
                values, vectors = fpca(self.paras['Cov.x'][c,...], l)
                xi_x = torch.randn(len(values), nsampling) * values[:,None]
                X[c,...] = (torch.matmul(torch.tensor(vectors), xi_x) + self.paras['mu.x'][c,:][:,None]).T
            f = torch.matmul(X, (self.elems['v']*self.elems['y']).T) - torch.matmul(torch.exp(X), (self.elems['v']*self.elems['ct']).T)
            f_max = torch.amax(f,dim=(0,1))
            omega = torch.mean(torch.exp(f - f_max[None,None,:]),axis=1) * self.paras['pi'][:,None]
            log_likelihood.append( torch.sum(torch.log(torch.sum(omega,axis=0)) + f_max) )
            omega = omega - torch.log(torch.sum(omega,axis=0))[None,:]
            id = torch.argmax(omega, axis=0)
            omega = omega/torch.sum(omega,axis=0)[None,:]
            new['pi'] = torch.mean(omega, axis=1)
            new['mu.x'] = torch.matmul(omega.float(), self.elems['b'])
            new['Cov.x'] = torch.matmul(self.elems['a'].permute(1,2,0),omega.T.float()).permute(2,0,1)

            for c in range(nclusters):
                new['Cov.x'][c,...] = torch.log(new['Cov.x'][c,...]*new['pi'][c]/ (torch.outer(new['mu.x'][c,],new['mu.x'][c,])) )
                new['mu.x'][c,...] = torch.log(new['mu.x'][c,...]/new['pi'][c]) - torch.diag(new['Cov.x'][c,...])/2

            if self.convergence(new, log_likelihood, reltol, abstol):
                break

            self.paras = copy.deepcopy(new)
            torch.cuda.empty_cache()
            
        return id, log_likelihood[-1]
