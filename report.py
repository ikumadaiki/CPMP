from gurobipy import *
from asyncio.constants import SENDFILE_FALLBACK_READBUFFER_SIZE
import random
from random import seed
from scipy.stats import multivariate_normal
import numpy as np

class CuttingPlaneAlgorithm:
    def __init__(self, ship_num, stack_num, height, n_init, size, size_, beta, mean, cov, OC, zk, aaa):
        self.ship_num=  ship_num
        self.stack_num = stack_num
        self.height = height
        self.n_init = n_init
        self.beta = beta
        self.size = size
        self.size_ = size_
        self.mean = mean
        self.cov = cov
        self.OC = OC
        self.zk = zk
        self.O1 = [i+1 for i in range(self.ship_num)]
        self.S = [i+1 for i in range(self.stack_num)]
        self.H = [i+1 for i in range(self.height)]
        self.P = [i+1 for i in range(self.ship_num)]
        self.f = stack_num*height-sum(self.n_init)
        self.aaa = aaa
        # 期待値と分散共分散行列の準備
        data_1 = np.random.multivariate_normal(mean, cov, size=self.size)

        O_ = np.argsort(data_1)
        O_ = O_+np.ones((size,ship_num)).astype(int)
        global order
        order = O_.tolist()

    def get_penalty(self,result):
        self.result=result
        self.penalty=[]
        for k in order:
            OO=k
            a=0
            for j in range(self.stack_num):
                for i in range(1,self.height):
                    for i_ in range(i+1,self.height+1):
                        if self.result[i-1][j]!=0 and self.result[i_-1][j]!=0:
                            if OO.index(self.result[i-1][j])>OO.index(self.result[i_-1][j]):
                                # print(j+1,height-i+1,O)
                                a+=1
                                # print("penalty!")
                                    # print(i,i_,j+1,O)
                                break
            self.penalty.append(a)
        self.zk=[i for i, x in enumerate(self.penalty) if x-self.alpha_tmp>0]
        
        for zk_ in self.zk:
            if not zk_ in self.aaa:
                self.aaa.append(zk_)
                self.OC.append(order[zk_])
        return (self.zk,self.OC,self.penalty)

    def solve(self):
        def add_cutting_plane(model,where):
            if where == GRB.callback.MIPSOL:
                c_var,x_var,x_tmp,d_var={},{},{},{}
                for var in self.model._vars:
                    if "u" in var.VarName:
                        u_var=var
                    if "alpha" in var.VarName:
                        alpha_var=var
                        self.alpha_tmp=model.cbGetSolution(var)
                    if "c_" in var.VarName:
                        s=int(var.VarName.split("_")[-3])
                        h=int(var.VarName.split("_")[-2])
                        i=int(var.VarName.split("_")[-1])
                        c_var[s,h,i]=var
                    if "d_" in var.VarName:
                        i=int(var.VarName.split("_")[-1])
                        d_var[i]=var
                    if "x_" in var.VarName:
                        s=int(var.VarName.split("_")[-3])
                        h=int(var.VarName.split("_")[-2])
                        p=int(var.VarName.split("_")[-1])
                        x_var[s,h,p]=var
                        x_tmp[s,h,p]=model.cbGetSolution(var)
                global result
                result=self.get_result(x_tmp)
                print(result)
                zk,OC,penalty=self.get_penalty(result)
                a=[]
                dic={}
                for i in OC:
                    # print(i)
                    if not i in a:
                        a.append(i)
                for j in a:
                    dic[tuple(j)]=OC.count(j)
                print("Len_dic=",len(dic))
                print(len(zk),len(OC))
                for s in self.S:
                    for h in range(2,len(self.H)+1):
                        for h_ in range(1,h):
                            for i,o in enumerate(dic):
                                for j,p in enumerate(o):
                                    model.cbLazy(c_var[s,h,i]>=dic[o]*(quicksum(x_var[s,h,k] for k in o[j:])-quicksum(x_var[s,h_,k] for k in o[j:])))
                for i in range(len(dic)):
                    model.cbLazy(d[i]>=quicksum(c_var[s,h,i] for s in self.S for h in self.H if h!=1)-alpha_var)
                model.cbLazy(u_var>=alpha_var+quicksum(d[i] for i in range(len(dic)))/((1-self.beta)*len(order)))
        self.model=Model("CuttingPlaneCVaRMinimization")
        x,c,d={},{},{}
        for s in self.S:
            for h in self.H:
                for p in self.P:
                    x[s,h,p]=self.model.addVar(vtype="B",name="x_"+str(s)+"_"+str(h)+"_"+str(p))

        for s in self.S:
            for h in range(2,len(self.H)+1):
                for i,o in enumerate(order):
                    c[s,h,i]=self.model.addVar(vtype="C",lb=0,name="c_"+str(s)+"_"+str(h)+"_"+str(i))

        for i in range(len(order)):
            d[i]=self.model.addVar(vtype="C",lb=0)
        alpha=self.model.addVar(vtype="C",lb=0,name="alpha")
        u=self.model.addVar(vtype="C",lb=0,name="u")

        self.model.update()
        self.model._vars=self.model.getVars()

        for p in self.P:
            self.model.addConstr(quicksum(x[s,h,p] for s in self.S for h in self.H)==self.n_init[p-1])

        for s in self.S:
            for h in self.H:
                self.model.addConstr(quicksum(x[s,h,p] for p in self.P)<=1)

        for s in self.S:
            for h in range(1,len(self.H)):
                self.model.addConstr(quicksum(x[s,h+1,p] for p in self.P)<=quicksum(x[s,h,p] for p in self.P))
        # if len(self.OC)==0:
        #     self.OC=[self.O1]
        # else:
        #     if self.O1 in self.OC:
        #         self.OC=self.OC.remove(self.O1)
        if len(self.OC)>self.size:
            self.OC=self.OC[self.size:]
        # print("OC=",self.OC)
        for s in self.S:
            for h in range(2,len(self.H)+1):
                for i,o in enumerate([self.O1]):
                    for j,p in enumerate(o):
                        self.model.addConstr(quicksum(x[s,h,k] for k in o[j:])<=quicksum(x[s,h-1,k] for k in o[j:]))

        for i in range(len([self.O1])):
            self.model.addConstr(d[i]>=quicksum(c[s,h,i] for s in self.S for h in self.H if h!=1)-alpha)

        self.model.addConstr(u>=alpha+quicksum(d[i] for i in range(len(order)))/((1-self.beta)*len(order)))
        self.model.setObjective(u)
        self.model._vars=self.model.getVars()
        self.model.params.LazyConstraints = 1
        self.model.params.TimeLimit = 3660
        if self.f>=self.height:
            self.model.optimize(add_cutting_plane)
        else:
            print("実行不可")

    def get_optimal_val(self):
        if self.model.Status==GRB.OPTIMAL:
            self.LB=self.model.ObjVal
            return self.model.ObjVal
        else:
            return None

    def get_optimal_sol(self):
        if self.model.Status==GRB.OPTIMAL:
            x_opt={}
            c_opt={}
            for var in self.model._vars:
                if "x_" in var.VarName:
                    # print(var)
                    s=int(var.VarName.split("_")[-3])
                    h=int(var.VarName.split("_")[-2])
                    p=int(var.VarName.split("_")[-1])
                    x_opt[s,h,p]=var.X
                if "c_" in var.VarName:
                    s=int(var.VarName.split("_")[-3])
                    h=int(var.VarName.split("_")[-2])
                    i=int(var.VarName.split("_")[-1])
                    c_opt[s,h,i]=var.X
                if "alpha" in var.VarName:
                    alpha_opt=var.X
                if "u" in var.VarName:
                    u_opt=var.X
            return (x_opt,c_opt,alpha_opt,u_opt)
        else:
            return None

    def get_result(self,x_tmp):
        self.x_tmp=x_tmp
        EPS=1.e-6
        self.result=np.zeros((self.height,self.stack_num))
        for (s,h,p) in self.x_tmp:
            if self.x_tmp[s,h,p]>EPS:
                self.result[self.height-h][s-1]=int(p)
        self.result=self.result.astype(int)
        return self.result
    
    # def get_penalty(self,result):
    #     if self.model.Status==GRB.OPTIMAL:
    #         self.result=result
    #         self.penalty=[]
    #         x,c,alpha,u=self.get_optimal_sol()
    #         for k in self.O:
    #             OO=k
    #             a=0
    #             for j in range(self.stack_num):
    #                 for i in range(1,self.height):
    #                     for i_ in range(i+1,self.height+1):
    #                         if self.result[i-1][j]!=0 and self.result[i_-1][j]!=0:
    #                             if OO.index(self.result[i-1][j])>OO.index(self.result[i_-1][j]):
    #                                 # print(j+1,height-i+1,O)
    #                                 a+=1
    #                                 # print("penalty!")
    #                                     # print(i,i_,j+1,O)
    #                                 break
    #             self.penalty.append(a)
    #         self.zk=[i for i, x in enumerate(self.penalty) if x-alpha>0]
    #         for zk_ in self.zk:
    #             self.OC.append(self.O[zk_])
    #         return (self.zk,self.OC,self.penalty)

def robust(ship_num,stack_num,height,n,Gamma,size_,mean,cov):
  # nと初期配置を変更しなければいけない
  num=sum(n)

  import numpy as np

  O1=[i+1 for i in range(ship_num)]

  Q=[i+1 for i in range(stack_num)]
  L=[i+1 for i in range(height)]
  P=[i+1 for i in range(ship_num)]
  I=[i+1 for i in range(num)]
  f=stack_num*height-len(I)
  a=1
  gamma=[]
  for i in n:
      for j in range(1,i+1):
          gamma.append(a)
      a+=1
  
  m=Model("BI")

  # 変数の定義
  alpha,beta={},{}
  for i in I:
      for q in Q:
          alpha[i,q]=m.addVar(vtype="B")
          beta[i,q]=m.addVar(vtype="B")
  J=[]
  for i in I:
      J.append([])
      for j in I:
          if gamma[i-1]<gamma[j-1]:
              if gamma[j-1]-gamma[i-1]<=Gamma:
                  J[i-1].append(j)

  for q in Q:
      m.addConstr(quicksum((alpha[i,q]+beta[i,q]) for i in I)<=len(L))

  for i in I:
      m.addConstr(quicksum((alpha[i,q]+beta[i,q]) for q in Q)==1)

  for i in I:
      for j in J[i-1]:
          for q in Q:
              m.addConstr(alpha[i,q]+alpha[j,q]+beta[j,q]<=1)

  m.setObjective(quicksum(beta[i,q] for i in I for q in Q))

  if f>=height:
      m.optimize()

  EPS=1.e-6

  if m.Status == GRB.OPTIMAL:
      print("====================================================")

      EPS=1.e-6
      a=[]
      for q in Q:
          a.append([])
      for (i,q) in alpha:
          if alpha[i,q].X>EPS:
              a[q-1].append(gamma[i-1])
      
      for (i,q) in beta:
          if beta[i,q].X>EPS:
              a[q-1].append(gamma[i-1])

      for q in Q:
          a[q-1]=sorted(a[q-1],reverse=True)

      global result_r
      result_r=np.zeros((height,stack_num))
      for q in Q:
          for i,r in enumerate(a[q-1]):
              result_r[height-i-1][q-1]=r
          # print(i,r)

      result_r=result_r.astype(int)
      
      print(result_r)
      print("the objective function", m.objVal)

      from asyncio.constants import SENDFILE_FALLBACK_READBUFFER_SIZE
      import random
      from random import seed
      import numpy as np
      from scipy.stats import multivariate_normal

      np.random.seed()
      data_1 = np.random.multivariate_normal(mean, cov, size=size_)

      O_=np.argsort(data_1)
      O_=O_+np.ones((size_,ship_num)).astype(int)
      OR=O_.tolist()


      global penalty_r
      penalty_r=[]
      for k in OR:
          O=k
          a=0
          for j in range(stack_num):
              for i in range(1,height):
                  for i_ in range(i+1,height+1):
                      if result_r[i-1][j]!=0 and result_r[i_-1][j]!=0:
                          if O.index(result_r[i-1][j])>O.index(result_r[i_-1][j]):
                              a+=1
                              # if Gamma ==2:
                                # print(j+1,height-i+1,O)
                              # print("penalty!")
                              # print(i,i_,j+1,O)
                              break
          penalty_r.append(a)
      
      penalty_r=np.sort(penalty_r)
      penalty_r=penalty_r[round(0.75*size_):]
      # print(penalty_r)

      # import matplotlib.pyplot as plt
      # plt.boxplot(penalty_r)

import seaborn as sns
sns.set()
def boxplot(result):
    data_2 = np.random.multivariate_normal(self.mean, self.cov, self.size, self.size_)

    OR_=np.argsort(data_2)
    OR_=OR_+np.ones((size_,len(inst[2]))).astype(int)
    OR=OR_.tolist()
    global penalty_cvar
    penalty_cvar=[]
    for k in OR:
        OR2=k
        a=0
        for j in range(inst[0]):
            for i in range(1,inst[1]):
                for i_ in range(i+1,inst[1]+1):
                    if result[i-1][j]!=0 and result[i_-1][j]!=0:
                        if OR2.index(result[i-1][j])>OR2.index(result[i_-1][j]):
                            # print(j+1,height-i+1,O)
                            a+=1
                            # print("penalty!")
                                # print(i,i_,j+1,O)
                            break
        penalty_cvar.append(a)
    penalty_cvar=np.sort(penalty_cvar)
    penalty_cvar=penalty_cvar[round(0.75*size_):]
    pcvar=[]
    pcvar+=list(penalty_cvar)
    data=(tuple(pcvar),)
    l=["75%-CVaR"]
    for g in range(1,6):
        prob=[]
        ship_num=len(inst[2])
        robust(ship_num,inst[0],inst[1],inst[2],g,size_,mean,cov)
        prob+=list(penalty_r)
        #   print(prob)
        data+=(tuple(prob),)
        l.append("Gamma="+str(g))
        # print("data=",data)
    import matplotlib.pyplot as plt
    plt.boxplot(data,labels=l)
    plt.show()
    # plt.hist(penalty_cvar,bins=50)
    # plt.show()

def hist(penalty_cvar):
  import matplotlib.pyplot as plt
  fig, ax = plt.subplots(2,3,figsize=(20,10),tight_layout=True)
  for g in range(1,6):
    prob=[]
    ship_num=len(inst[2])
    robust(ship_num,inst[0],inst[1],inst[2],g,size_,mean,cov)
    prob+=list(penalty_r)
    # print(prob)
    # data+=(tuple(prob),)
    # l.append("Gamma="+str(g))
    # print("data=",data)
  # plt.boxplot(data,labels=l)
  # plt.show()

    ax[int(g/3)][g%3].hist(prob,bins=5,range=(0,5))
    ax[int(g/3)][g%3].set_title("Gamma="+str(g))
    plt.xlabel("荷繰りの数",fontname="MS Gothic")
  ax[0][0].hist(penalty_cvar,bins=5,range=(0,5))
  ax[0][0].set_title("75%-CVaR")
  plt.xlabel("荷繰りの数",fontname="MS Gothic")
  plt.show()
  print("len=",len(prob))