from problems import ETC, EE
import numpy as np
import sympy
from pymoo.optimize import minimize
from pymoo.algorithms.moo.nsga2 import NSGA2
# from pymoo.operators.crossover.sbx import SBX
# from pymoo.operators.mutation.pm import PM
from pymoo.operators.sampling.rnd import FloatRandomSampling
from pymoo.core.problem import ElementwiseProblem
from pymoo.core.problem import Problem
from pymoo.visualization.scatter import Scatter

class MyProblem(ElementwiseProblem):

    def __init__(self):
        super().__init__(n_var=1,
                         n_obj=2,
                         n_ieq_constr=0,
                         xl=0.6,
                         xu=0.8)

    def _evaluate(self, r, out, *args, **kwargs):

        R = sympy.Symbol('R')
        t = sympy.Symbol('t')

        alpha = 110
        beta = 3

        b = 1.25
        a = 0.95

        T = [0]
        li = [(beta / alpha) * (t / alpha) ** (beta - 1) * sympy.exp(-(t / alpha) ** beta)]

        for i in range(1, 12):
            T.append(alpha * (-sympy.log((sympy.log(R) / b ** i) + sympy.exp(-(a * T[i - 1] / alpha) ** beta))) ** (
                        1 / beta) - a * T[i - 1])

        for j in range(1, 12):
            li.append(b ** j * (beta / alpha) * ((t + a * T[j - 1]) / alpha) ** (beta - 1) * sympy.exp(
                -((t + a * T[j - 1]) / alpha) ** beta))

        def Sum_of_T():
            sum = T[0]
            for i in range(1, 4):
                sum += T[i]
            return sum

        def C_pm_m():
            C_pm_s = 400
            gamma = 1.2
            sum = 0
            for i in range(1, 3):
                sum += C_pm_s + gamma * li[i].subs(t, T[i])
            return sum

        def T_pm_m():
            t_pm_m = 1
            return 2 * t_pm_m

        def C_cm_m(R):
            C_cm_s = 150
            return C_cm_s * 3 * -sympy.log(R)

        def T_cm_m(R):
            t_cm_m = 0.5
            return t_cm_m * 3 * -sympy.log(R)

        def p_avr():
            p_0 = 0.008
            u = 0.075
            sigma = 20
            theta = 1.1
            p_avr = [0]
            for i in range(1, len(li)):
                p_avr.append(((p_0 + u * (1 - sympy.exp(-sigma * li[i].subs(t, 0) ** theta))) + (
                            p_0 + u * (1 - sympy.exp(-sigma * li[i].subs(t, T[i]) ** theta)))) / 2)
            return p_avr

        def C_d():
            C_d_v = 10
            v = 100
            pro = 0
            for i in range(1, 4):
                pro += T[i] * p_avr()[i] * v
            return C_d_v * pro

        C_r = 2000
        T_r = 2

        def ETC(r):
            ts = (C_pm_m() + C_cm_m(R) + C_r + C_d()).subs(R, r)
            ms = (Sum_of_T() + T_pm_m() + T_cm_m(R) + T_r).subs(R, r)
            return ts / ms


        omega = 1.2

        Xi = [0]
        h = sympy.Symbol('h')

        for i in range(1, 12):
            Xi.append(- omega * b ** i * sympy.exp(-h))

        def E_o():
            sum = 0
            for i in range(1, 4):
                sum += 180 * T[i] + Xi[i].subs(t, T[i]) - Xi[i].subs(t, 0)
            return sum

        def E_i():
            E_iv = 25
            t_pm = 1
            return E_iv * 2 * t_pm

        def E_w():
            E_wv = 15
            return E_wv * 3 * -sympy.log(R)

        def E_p():
            E_pv = 50
            E_pr = 300
            return E_pv * (5 * -sympy.log(R)) + E_pr

        def Y_y():
            v = 100
            t_i = 0.4
            sum = 0
            sub = 0
            for i in range(1, 4):
                sum += T[i] * v

            for j in range(1, 4):
                sub += v * p_avr()[j] * T[j] * (1 - t_i)

            return (sum - sub)

        def EE(r):
            ts = Y_y().subs(R, r)
            ms = (E_o() + E_i() + E_w() + E_p()).subs(R, r)
            return ts / ms

        f1 = ETC(r)
        f2 = EE(r)


        out["F"] = [f1, -f2]

problem = MyProblem()

algorithm = NSGA2(
    pop_size=100,
    crossover=0.9,
    mutation=0.1,
    eliminate_duplicates=True
)

termination = ("n_gen", 100)

res = minimize(problem,
               algorithm,
               termination,
               seed=1,
               # save_history=True,
               verbose=False
)


plot = Scatter()
plot.add(problem.pareto_front(), plot_type="line", color="black", alpha=0.7)
plot.add(res.F, facecolor="none", edgecolor="red")
plot.show()
