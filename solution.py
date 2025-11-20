from sage.all import *
import matplotlib.pyplot as plt

# Составление системы уравнений
class Equasions():  
    def __init__(self):

        E, J, K, P0, M5, q12, q21, x = var('E J K P0 M5 q12 q21 x')

        self.P0 = P0    # Известная сила
        self.M5 = M5    # Известный момент
        self.K = K      # Коэффициент упругости
        self.E = E      # Модуль упругости
        self.J = J      # Момент инерции 
        self.x = x
        self.q12 = q12  # Распределённая сила
        self.q21 = q21  # Распределённая сила
        self.k = (q21 - q12) / 9
        self.points = [0, 2, 7, 10, 15, 23, 12, 12, 21, 21] 
        self.coefs = [self.x - p for p in self.points]

        self.M = self.M_fun()           # момент
        self.Q = self.Q_fun()           # сила
        self.theta = self.theta_fun()   # угол
        self.w = self.w_fun()           # прогиб
        self.R23 = self.R23_fun()       # R23 = K * w(23); можно явно вывести через остальные R_i, p
        self.sp = self.sum_powers()
        self.sm = self.sum_momentums()

    def __str__(self):
        str1 = f'w(x) = {self.w}'
        str2 = f'theta(x) = {self.theta}'
        str3 = f'M(x) = {self.M}'
        str4 = f'Q(x) = {self.Q}'
        str5 = f'Moms_1(x) = {self.sm}'
        str6 = f'Pows(x) = {self.sum_powers()}'
        str7 = f'R23 = {self.R23}'
        return f'THE SYSTEM OF THE EQUASIONS: \n\n {str1} \n\n {str2} \n\n {str3} \n\n {str4}  \n\n {str5} \n\n {str6} \n\n {str7} \n\n---------------------------------------------------------------------------------------------------------'
    

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
    
    def R23_fun(self):
        return self.K * self.w(x = 23)
    
    def sum_powers(self):
        return self.Q(x = 25)

    def sum_momentums(self):
        return self.M(x = 25)
    
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

# Решение системы уравнений
class SystemSolution():
    def __init__(self, p = 0, m = 0, k = 1, e = 1, j = 1, q1 = 0, q2 = 0):
        self.P0, self.M5, self.K, self.E, self.J, self.q12, self.q21 = p, m, k, e, j, q1, q2
        self.points = [2, 7, 10, 15]
        self.equasion = Equasions()
        self.R23 = self.equasion.R23
        self.w = self.equasion.w(R23 = self.R23)
        self.theta = self.equasion.theta(R23 = self.R23)
        self.M = self.equasion.M(R23 = self.R23)
        self.Q = self.equasion.Q(R23 = self.R23)
        self.sp = self.equasion.sp(R23 = self.R23)
        self.sm = self.equasion.sm(R23 = self.R23)

        self.system = self.make_system()

    def __str__(self):
        return "THE TRANSFORMED EQUASIONS:\n\n" + "\n\n".join(map(str, self.system)) + "-----------------------------------------------------------------------------------------------------------------"
    

    def make_system(self):
        eq = []
        for p in self.points:
            temp = self.w(R23 = self.R23, x = p)
            eq.append(temp(P0 = self.P0, M5 = self.M5, K = self.K, E = self.E, J = self.J, q12 = self.q12, q21 = self.q21))
        temp = self.sm(R23 = self.R23, x = 25)
        eq.append(temp(P0 = self.P0, M5 = self.M5, K = self.K, E = self.E, J = self.J, q12 = self.q12, q21 = self.q21))
        temp = self.sp(R23 = self.R23, x = 25)
        eq.append(temp(P0 = self.P0, M5 = self.M5, K = self.K, E = self.E, J = self.J, q12 = self.q12, q21 = self.q21))
        return eq
    

    def solve(self):
        return solve(self.system)
    
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

class DrawGraphics():
    def __init__(self, system, p = 1, m = 4, k = 0.4, e = 200000, j = 1, q1 = 0, q2 = 0):
        self.M = m
        self.P = p
        self.system = system
        self.solution = self.system.solve()
        self.solve_dict = {str(s).split('==')[0]: str(s).split('==')[1] for s in self.solution[0]}
        for key in self.solve_dict:
            nums = self.solve_dict[key].split('/')
            self.solve_dict[key] = float(nums[0][2:]) / float(nums[1][:-1])
        
        self.w = self.system.w(**self.solve_dict, P0 = p, M5 = m, K = k, E = e, J = j, q12 = q1, q21 = q2)
        self.theta = self.system.theta(**self.solve_dict, P0 = p, M5 = m, K = k, E = e, J = j, q12 = q1, q21 = q2)
        self.M = self.system.M(**self.solve_dict, P0 = p, M5 = m, K = k, E = e, J = j, q12 = q1, q21 = q2)
        self.Q = self.system.Q(**self.solve_dict, P0 = p, M5 = m, K = k, E = e, J = j, q12 = q1, q21 = q2)
        
    def make_graphics(self):
        print(self.solve_dict)

        X = [x for x in range(26)]
        W = [self.w(x = x) for x in X]
        THETA = [self.theta(x = x) for x in X]
        M = [self.M(x = x) for x in X]
        Q = [self.Q(x = x) for x in X]

        _, axes = plt.subplots(nrows=4, ncols=1, figsize=(8, 12)) 

        for i, (graphic, name) in enumerate(zip([Q, M, THETA, W], ['Q', 'M', '$\Theta$', 'W'])):
            ax = axes[i]
            
            ax.axhline(y=0, color='black', linewidth=0.8)
            
            ax.plot([2, 7, 10, 15], [0]*4, 'ro', label='R')   
            ax.plot([5], [0], 'ko', label='M5 ')             
            ax.plot([0], [0], 'go', label='P0 ')             
            ax.plot([23], [0], 'yo', label='K23')
            ax.plot([12, 21], [0]*2, 'm', label='q12-q21')               
            
            ax.plot(X, graphic, 'b-', label=name)
            
            ax.set_xlabel('x')
            ax.set_ylabel(name)
            ax.set_title(f'{name}(x)', fontsize=12)
            ax.legend(loc="upper right")
            ax.grid(True)

        plt.tight_layout() 
        plt.show()


# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

eq = Equasions()
print(eq)

system_sol = SystemSolution(p = 10, m = 47, k = 0.4, e = 2000, j = 1, q1 = 1, q2 = 2)
print(system_sol)
solution = system_sol.solve()
print("Solution: \n\n", solution)
DrawGraphics(system_sol, p = 10, m = 47, k = 0.4, e = 2000, j = 1, q1 = 1, q2 = 2).make_graphics()
