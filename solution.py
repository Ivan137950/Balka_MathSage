from sage.all import *


class Equasions():
    def __init__(self):
        E, J, P0, M5, q12, q21, x = var('E J P0 M5 q12 q21 x')
        self.P0 = P0    
        self.M5 = M5
        self.E = E
        self.J = J
        self.x = x
        self.q12 = q12
        self.q21 = q21
        self.k = (q21 - q12) / 9
        self.points = [0, 2, 7, 10, 15, 23, 12, 12, 21, 21] 
        self.coefs = [self.x - p for p in self.points]

        self.M = self.M_fun()
        self.Q = self.Q_fun()
        self.theta = self.theta_fun()     
        self.w = self.w_fun()
        self.sp = self.sum_powers()
        self.sm_1 = self.sum_momentums_1()
        self.sm_2 = self.sum_momentums_2()

    def __str__(self):
        str1 = f'w(x) = {self.w}'
        str2 = f'theta(x) = {self.theta}'
        str3 = f'M(x) = {self.M}'
        str4 = f'Q(x) = {self.Q}'
        str5 = f'Moms_1(x) = {self.sm_1}'
        str6 = f'Moms_2(x) = {self.sm_2}'
        str7 = f'Pows(x) = {self.sum_powers()}'
        return f'THE SYSTEM OF THE EQUASIONS: \n\n {str1} \n\n {str2} \n\n {str3} \n\n {str4}  \n\n {str5} \n\n {str6} \n\n {str7} \n---------------------------------------------------------------------------------------------------------'
    

    def r_x(self):
        r = var('r')
        return r
    
    def q_x(self):
        r = var('r')
        return r ** 2/2
    
    def equation(self):
        return [self.r_x() for _ in range(6)] + [self.q_x()] + [integrate(self.q_x(), r)] + [self.q_x()] + [integrate(self.q_x(), r)]
    

    def M_fun(self):        
        (R2, R7, R10, R15, R23) = var('R2 R7 R10 R15 R23')
        alphas = [self.P0, R2, R7, R10, R15, R23, self.q12, -self.k,-self.q21, self.k]
        eq = 0
        for i in range(len(alphas)):
            eq += alphas[i] * self.equation()[i](r = self.coefs[i]) * unit_step(self.coefs[i])
        eq += self.M5
        return eq
    
    def Q_fun(self):
        (R2, R7, R10, R15, R23) = var('R2 R7 R10 R15 R23')
        alphas = [self.P0, R2, R7, R10, R15, R23, self.q12, -self.k,-self.q21, self.k]
        eq = 0
        for i in range(len(alphas)):
            eq += alphas[i] * diff(self.equation()[i], r)(r = self.coefs[i]) * unit_step(self.coefs[i])
        return eq
    
    def theta_fun(self):
        (R2, R7, R10, R15, R23) = var('R2 R7 R10 R15 R23')
        alphas = [self.P0, R2, R7, R10, R15, R23, self.q12, -self.k,-self.q21, self.k]
        eq = 0
        for i in range(len(alphas)):
            eq += alphas[i] * integrate(self.equation()[i], r)(r = self.coefs[i]) * unit_step(self.coefs[i])
        eq += self.M5
        theta = var('theta')
        eq /= (self.J * self.E)
        eq += theta
        return eq
    
    def w_fun(self):
        (R2, R7, R10, R15, R23) = var('R2 R7 R10 R15 R23')
        alphas = [self.P0, R2, R7, R10, R15, R23, self.q12, -self.k,-self.q21, self.k]
        eq = 0
        for i in range(len(alphas)):
            eq += alphas[i] * integrate(integrate(self.equation()[i], r), r)(r = self.coefs[i]) * unit_step(self.coefs[i])
        eq += self.M5
        theta, w = var('theta w')
        eq /= (self.J * self.E)
        eq += (theta * x + w)
        return eq
    
    def sum_powers(self):
        return self.Q(x = 25)

    def sum_momentums_1(self):
        return self.M(x = 25)
    
    def sum_momentums_2(self):
        (R2, R7, R10, R15, R23) = var('R2 R7 R10 R15 R23')
        return 2 * R2 + 7 * R7 + 10 * R10 + 15 * R15 + 23 * R23 + self.M5 + 1/27*(self.q12 - self.q21)*(12) ** 3 - 1/27*(self.q12 - self.q21)*(21) ** 3 + self.q12*(12) ** 2 - self.q21*(21) ** 2
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------


class SystemSolution():
    def __init__(self, points, p=1, m=1, q12=1, q21=2):
        self.equation = Equasions()
        self.w = self.equation.w
        self.points = points
        self.summa = self.equation.sp
        self.moms = self.equation.M
        self.p = p
        self.m = m
        self.q12 = q12
        self.q21 = q21
        self.system = []
        self.fullfill_system()

    def __str__(self):
       strs = "SYSTEM: \n\n" + '\n'.join(map(str, simplify(self.system)))
       strs += '\n\n' + "SOLUTION: \n\n" + '\n'.join(map(str, self.solve_system()))
       strs += "\n\n ---------------------------------------------------------------------- \n\n"
       return strs
    
    def fullfill_system(self):
        for p in self.points:
            self.system.append(self.w(x = p, M6 = self.m, P0 = self.p, q12 = self.q12, q21 = self.q21) == 0)
        self.system.append(self.moms(x = 26, M6 = self.m, P0 = self.p, q12 = self.q12, q21 = self.q21) == 0)
        self.system.append(self.equation.sm(M6 = self.m, P0 = self.p, q12 = self.q12, q21 = self.q21) == 0)
        self.system.append(self.summa(M6 = self.m, P0 = self.p, q12 = self.q12, q21 = self.q21) == 0)

    
    def solve_system(self):
        #print('\n'.join(map(str, simplify(self.system))))
        solution = solve(self.system, R3, R8, R11, R16, R24, theta, w0)
        return solution
    

    #def float_solution(self):
    #    keys = ['R3', 'R8', 'R11', 'R16', 'R24', 'theta', 'w0']
    #    solution = []
    #    ss = sol.solve_system()
    #    for s in map(str, ss):
    #        float_s = s.split('==')[1][2:-1]
    #        solution.append(float(float_s.split('/')[0]) / float(float_s.split('/')[1]))
    #    return {keys[i]: solution[i] for i in range(len(keys))}



# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------


#P0, M6, K24,  x = var('P0 M6 K24 x')

eq = Equasions()
print(eq)

#points = [3, 8, 11, 16]
#p, m, q1, q2 = 0.1, 0.02, 0.02, 0.04 
#sol = SystemSolution(points, p, m, q1, q2)
##solut = sol.float_solution()
#print(sol)
#print(solut)

#values = [{'P0': p, 
#           'M6': m, 
#           'q12': q1, 
#           'q21': q2, 
#           'R3': solut['R3'], 
#           'R8': solut['R8'], 
#           'R11': solut['R11'], 
#           'R16': solut['R16'], 
#           'R24': solut['R24'], 
#           'theta': solut['theta'], 
#           'w0': solut['w0'], 
#           'x': xi}
#           for xi in range(0, 26)
#           ]

#w_points = [eq.w(**(val)) for val in values]
#theta_points = [eq.theta(**(val)) for val in values]
#m_points = [eq.M(**(val)) for val in values]
#q_points = [eq.Q(**(val)) for val in values]
#x_points = [val['x'] for val in values]

#print(w_points)
#print(theta_points)
#print(m_points)
#print(q_points)

#import matplotlib.pyplot as plt

#labels = ['w', '$\Theta$', 'm', 'q']
#for i, dataset in enumerate([w_points, theta_points, m_points, q_points]):
#    plt.plot(x_points, [0 for _ in range(26)], 'k')
#    plt.plot([3, 8, 11, 16, 24], [0 for _ in range(5)], 'ro')
#    plt.plot(x_points, dataset, 'b')
#    plt.xlabel('x')
#    plt.ylabel(labels[i])
#    plt.title(labels[i])
#    plt.legend()
#    plt.grid()
#    plt.show()