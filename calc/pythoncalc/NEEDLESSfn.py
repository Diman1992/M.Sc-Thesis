import numpy as np

#massen
mu, mc, mt = 2, 1275, 173070
md, ms, mb = 5, 95, 4180

#params
eps = 0.2

print("massenratios ", mt/mt, mc/mt, mu/mt)
print("epsilons ", eps**0, eps**4, eps**6)

#matrix
Yu = np.matrix([[eps**6, eps**5, eps**3],[ eps**5, eps**4, eps**2],[ eps**3, eps**2, eps**0]])
Yd = np.matrix([[eps**4, eps**3, eps**2],[ eps**3, eps**2, eps**1],[ eps**2, eps**1, eps**0]])
YuInv = np.linalg.inv(Yu)
#YdInv = np.linalg.inv(Yd)
print(YuInv)

V1 = np.dot(Yd, YuInv)
print(V1)
