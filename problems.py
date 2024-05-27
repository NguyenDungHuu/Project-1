import numpy as np
import sympy
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import cm
from matplotlib.ticker import LinearLocator

#N = sympy.Symbol('N')
R = sympy.Symbol('R')
t = sympy.Symbol('t')

alpha = 110
beta = 3

b = 1.25
a = 0.2



T = [0]
li = [(beta/alpha)*(t/alpha)**(beta-1)]



for i in range(1, 12):
    T.append(alpha*((-sympy.log(R)/b**i) + (a*T[i-1]/alpha)**beta)**(1/beta)-a*T[i-1])

# for i in range(1, 12):
#     T[i] = T[i].subs(R, 0.7)
# xt = range(0, 12)
# plt.scatter(xt, T)
# plt.show()



for j in range(1, 12):
    li.append(b**j*(beta/alpha)*((t+a*T[j-1])/alpha)**(beta-1))

def Sum_of_T(N, r):
    sum = T[0]
    for i in range(1, N+2):
        sum += T[i]
    return sum.subs(R, r)
def C_pm_m(N, r):
    C_pm_s = 400
    gamma = 1.2
    sum = 0
    for i in range(1, N+1):
        sum += C_pm_s + gamma*li[i].subs(t, T[i])
    return sum.subs(R, r)

def T_pm_m(N):
    t_pm_m = 1
    return N*t_pm_m

def C_cm_m(N, r):
    C_cm_s = 150
    return C_cm_s*(N+1)*-sympy.log(R).subs(R, r)

def T_cm_m(N, r):
    t_cm_m = 0.5
    return t_cm_m*(N+1)*-sympy.log(R).subs(R, r)

# def p():
#     p_0 = 0.008
#     u = 0.075
#     sigma = 20
#     theta = 1.1
#     p = [0]
#     for i in range(1,len(li)):
#         p.append(p_0 + u*(1 - sympy.exp(-sigma*li[i]**theta)))
#     return p

# def p_avr(r):
#     p_0 = 0.008
#     u = 0.075
#     sigma = 20
#     theta = 1.1
#     p_avr = [0]
#     for i in range(1,len(li)):
#         p_avr.append(((p_0 + u*(1 - sympy.exp(-sigma*li[i].subs(t, 0).subs(R, r)**theta))) + (p_0 + u*(1 - sympy.exp(-sigma*li[i].subs(t, T[i]).subs(R, r)**theta))))/2)
#     return p_avr

p_0 = 0.008
u = 0.075
sigma = 20
theta = 1.1
def C_d(N, r):
    C_d_v = 10
    v = 100
    pro = 0
    for i in range(1, N+2):
        pro += T[i]*v*((p_0 + u*(1 - sympy.exp(-sigma*li[i].subs(t, 0).subs(R, r)**theta))) + (p_0 + u*(1 - sympy.exp(-sigma*li[i].subs(t, T[i]).subs(R, r)**theta))))/2
    return C_d_v*pro.subs(R, r)

C_r = 2000
T_r = 2

def ETC(n, r):
    ts = (C_pm_m(n, r) + C_cm_m(n, r) + C_r + C_d(n, r))
    ms = (Sum_of_T(n, r) + T_pm_m(n) + T_cm_m(n, r) + T_r)
    return ts/ms




# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# X = np.arange(1, 8, 1)
# Y = np.arange(0.4, 1, 0.1)
# #
# X, Y = np.meshgrid(X, Y)
# Z = np.zeros(X.shape)
#
#
# for i in range(0, X.shape[0]):
#     for j in range(0, X.shape[1]):
#         Z[i, j] = ETC(X[i,j], Y[i,j])
#
#
#
# surf = ax.plot_surface(X, Y, Z, cmap= cm.jet, linewidth = 0.2, rstride= 1, cstride = 1)
# ax.set_xlabel("N")
# ax.set_ylabel("R")
# ax.set_zlabel("ETC")
# #
# #
# plt.show()

xi = [180]
omega = 1.2
for i in range(1, 12):
    xi.append(xi[0] + omega*li[i])



Xi = [0]
h = sympy.Symbol('h')

for i in range(1, 12):
    Xi.append(omega*b**i*((t+a*T[i-1])/alpha)**beta)



def E_o(N, r):
    sum = 0
    for i in range(1, N+2):
        sum += 180 * T[i].subs(R, r) + Xi[i].subs(t, T[i]).subs(R, r) - Xi[i].subs(t, 0).subs(R, r)
    return sum



def E_i(N):
    E_iv = 25
    t_pm = 1
    return E_iv*N*t_pm

def E_w(N, r):
    E_wv = 15
    return E_wv*(N+1)*-sympy.log(R).subs(R, r)

def E_p(N, r):
    E_pv = 50
    E_pr = 300
    return E_pv*(N+(N+1)*-sympy.log(R)).subs(R, r) + E_pr

def Y_y(N, r):
    v = 100
    t_i = 0.4
    sum = 0
    sub = 0
    for i in range(1, N+2):
        sum += T[i]*v

    for j in range(1, N+2):
        sub += v*((p_0 + u*(1 - sympy.exp(-sigma*li[j].subs(t, 0).subs(R, r)**theta))) + (p_0 + u*(1 - sympy.exp(-sigma*li[j].subs(t, T[j]).subs(R, r)**theta))))/2*T[j]*(1-t_i)

    return (sum - sub).subs(R, r)

def EE(n, r):
    ts = Y_y(n, r)
    ms = (E_o(n, r) + E_i(n) + E_w(n, r) + E_p(n, r))
    return ts/ms



# X = np.arange(1, 11, 1)
# Y = np.arange(0.1, 1, 0.1)
# # S = []
# # for i in range(0, len(X)):
# #     print("& "+str(round(EE(X[i], 0.9), 6)) + " ", end='')
#
# # for i in range(0, len(X)):
# #     for j in range(0, len(Y)):
# #         S.append(round(EE(X[i], Y[j]), 6))
# # Z = np.array(S)
# # print(max(Z))
# # # print(min(S))
# # # #
# # # Z1 = []
# # # Z2 = []
# # # M = 100
# # # for i in range(0, len(X)):
# # #     for j in range(0, len(Y)):
# # #         Z1.append(ETC(X[i], Y[j]))
# # # for i in range(0, len(Z1)):
# # #     if Z1[i] <= M:
# # #         M = Z1[i]
# # # print(M)
# X, Y = np.meshgrid(X, Y)
# Z1 = np.zeros(X.shape)
# Z2 = np.zeros(X.shape)
# # #
# # for i in range(0, X.shape[0]):
# #     for j in range(0, X.shape[1]):
# #         print(X[i, j], Y[i,j], ETC(X[i, j], Y[i,j]), EE(X[i, j], Y[i,j]))
# #
# #
# for i in range(0, X.shape[0]):
#     for j in range(0, X.shape[1]):
#         Z1[i, j], Z2[i, j] = ETC(X[i,j], Y[i,j]) ,EE(X[i,j], Y[i,j])
# #
# #
# fig = plt.figure()
#
# ax = fig.add_subplot(1, 2, 1, projection='3d')
# ax.plot_surface(X, Y, Z1, cmap= cm.jet, linewidth = 0.2, rstride= 1, cstride = 1)
# ax.set_xlabel("N")
# ax.set_ylabel("R")
# plt.title("ETC")
#
# ax = fig.add_subplot(1, 2, 2, projection='3d')
# ax.plot_surface(X, Y, Z2, cmap= cm.jet, linewidth = 0.2, rstride= 1, cstride = 1)
# ax.set_xlabel("N")
# ax.set_ylabel("R")
# plt.title("EE")
#
# plt.show()






