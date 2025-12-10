from sage.all import *
import matplotlib.pyplot as plt
import numpy as np
import subprocess


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
        self.k = (q21 - q12) / 9.
        self.points = [0., 2., 7., 10., 15., 23., 12., 12., 21., 21.] 
        self.coefs = [self.x - p for p in self.points]

        self.w = self.w_fun()           # прогиб
        self.theta = self.theta_fun()   # угол
        self.M = self.M_fun()           # момент
        self.Q = self.Q_fun()           # сила

        self.R23 = self.R23_fun()       # R23 = K * w(23); можно явно вывести через остальные R_i, p
        self.sp = self.Q_25()
        self.sm = self.Mom_25()

    def __str__(self):
        str1 = f'w(x) = {self.w}'
        str2 = f'theta(x) = {self.theta}'
        str3 = f'M(x) = {self.M}'
        str4 = f'Q(x) = {self.Q}'
        str5 = f'Moms_1(x) = {self.sm}'
        str6 = f'Pows(x) = {self.Q_25()}'
        str7 = f'R23 = {self.R23}'
        return f'THE SYSTEM OF THE EQUASIONS: \n\n {str1} \n\n {str2} \n\n {str3} \n\n {str4}  \n\n {str5} \n\n {str6} \n\n {str7} \n\n---------------------------------------------------------------------------------------------------------'
    
    def w_fun(self):
        (R2, R7, R10, R15, R23) = var('R2 R7 R10 R15 R23')
        theta, w = var('theta w')
        eq =  (self.x) ** 3 * self.P0 / 6 
        eq += self.M5 * (self.x - 5) ** 2 / 2 * unit_step(self.x - 5)
        for point, coef in zip([2., 7., 10., 15., 23.], [R2, R7, R10, R15, R23]):
             eq += ((self.x - point) ** 3 / 6 * unit_step(self.x - point)) * coef
        eq += (self.q12 * (self.x - 12) ** 4 / 24 + self.k * (self.x - 12) ** 5 / 120) * unit_step(self.x - 12)
        eq -= (self.q21 * (self.x - 21) ** 4 / 24 + self.k * (self.x - 21) ** 5 / 120) * unit_step(self.x - 21)
        eq  /= (self.J * self.E)
        eq += w  + theta * self.x
        return eq

    def theta_fun(self):
        (C2, C5, C7, C10, C15, C23, C12, C21) = var('C2 C5 C7 C10 C15 C23 C12 C21')
        eq2c = {unit_step(self.x - a): C for a, C in zip([2., 5., 7., 10., 15., 23., 12., 21.], 
                                                        [C2, C5, C7, C10, C15, C23, C12, C21])}
        c2eq = {eq2c[k]: k for k in eq2c}
        
        expr_with_C = self.w.subs(eq2c)
        diff_expr = diff(expr_with_C, self.x)
        final_result = diff_expr.subs(c2eq)
        return final_result

    def M_fun(self):        
        (C2, C5, C7, C10, C15, C23, C12, C21) = var('C2 C5 C7 C10 C15 C23 C12 C21')
        eq2c = {unit_step(self.x - a): C for a, C in zip([2., 5., 7., 10., 15., 23., 12., 21.], 
                                                         [C2, C5, C7, C10, C15, C23, C12, C21])}
        c2eq = {eq2c[k]: k for k in eq2c}
        
        expr_with_C = self.theta.subs(eq2c)
        diff_expr = diff(expr_with_C, self.x)
        final_result = diff_expr.subs(c2eq)
        return (self.E * self.J) * final_result
        
    def Q_fun(self):
        (C2, C5, C7, C10, C15, C23, C12, C21) = var('C2 C5 C7 C10 C15 C23 C12 C21')
        eq2c = {unit_step(self.x - a): C for a, C in zip([2., 5., 7., 10., 15., 23., 12., 21.], 
                                                         [C2, C5, C7, C10, C15, C23, C12, C21])}
        c2eq = {eq2c[k]: k for k in eq2c}
        
        expr_with_C = self.M.subs(eq2c)
        diff_expr = diff(expr_with_C, self.x)
        final_result = diff_expr.subs(c2eq)
        return final_result

    def R23_fun(self):
        return self.K * self.w(x = 23.)
    
    def Q_25(self):
        return self.Q(x = 25.)

    def Mom_25(self):
        return self.M(x = 25.)
    
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

# Решение системы уравнений
class SystemSolution():
    def __init__(self, p = 0, m = 0, k = 1, e = 1, j = 1, q1 = 0, q2 = 0):
        self.P0, self.M5, self.K, self.E, self.J, self.q12, self.q21 = p, m, k, e, j, q1, q2
        self.points = [2., 7., 10., 15.]
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
            eq.append(temp(P0 = self.P0, M5 = self.M5, K = self.K, E = self.E, J = self.J, q12 = self.q12, q21 = self.q21 == 0))
        temp = self.sm(R23 = self.R23, x = 25.)
        eq.append(temp(P0 = self.P0, M5 = self.M5, K = self.K, E = self.E, J = self.J, q12 = self.q12, q21 = self.q21) == 0)
        temp = self.sp(R23 = self.R23, x = 25.)
        eq.append(temp(P0 = self.P0, M5 = self.M5, K = self.K, E = self.E, J = self.J, q12 = self.q12, q21 = self.q21) == 0)
        return eq
    

    def solve(self):
        return solve(self.system)
    
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

#Точки для рисования графиков
class GraphicsData():
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


    def __str__(self):
        s = "COEFS:\n\n"
        s += str(self.solve_dict)
        s += "Equations:\n\n"
        s += 'w =' + str(self.w) + "\n\n"
        s += 'theta =' + str(self.theta) + "\n\n"
        s += 'M =' + str(self.M) + "\n\n"
        s += 'Q =' + str(self.Q) + "\n\n"
        return s
        
    def make_graphics_dataset(self):

        X = [x/1000 for x in range(25000)]
        W = [self.w(x = x) for x in X]
        THETA = [self.theta(x = x) for x in X]
        M = [self.M(x = x) for x in X]
        Q = [self.Q(x = x) for x in X]
        print(self.w(x = 2.), self.w(x = 7.), self.w(x = 10.), self.w(x = 15.))

        np.savetxt('W_sage.txt', W)
        np.savetxt('Theta_sage.txt', THETA)
        np.savetxt('M_sage.txt', M)
        np.savetxt('Q_sage.txt', Q)
        np.savetxt('X_sage.txt', X)


# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

#  Решение SAGE
class SageSolution():
    def __init__ (self, params):
        self.params = params

    def run_sage_solution(self):
        try:
            # Составление системы
            eq = Equasions()
            print(eq)
            
            # Решение системы
            system_sol = SystemSolution(**self.params)
            print(system_sol)
            solution = system_sol.solve()
            print("Solution: \n\n", solution)

            # Данные для графиков
            gd = GraphicsData(system_sol, **self.params)
            gd.make_graphics_dataset()

            print("Решение в Sage успешно получено")
        except:
            print("не удалось выполнить операцию")
            

#  Решение МКЭ
class CppSolution():
    def __init__(self, params):
        self.params = params

    
    def run_cpp_solution(self):
        try:
            # Запуск C++ программы с параметрами
            result = subprocess.run(
                ["./functions.exe"] + list(map(str, self.params.values())),
                capture_output=True,
                text=True
            )

            print("Решение в cpp: \n", result.stdout)
            print("Решение в cpp успешно получено")
        
        except:
            print("не удалось выполнить операцию")


# Графики
class Graphic():
    def __init__(self):
        pass
    
    def draw_graphics(self):
        W_sage = np.loadtxt('W_sage.txt')
        Theta_sage = np.loadtxt('Theta_sage.txt')
        M_sage = np.loadtxt('M_sage.txt')
        Q_sage = np.loadtxt('Q_sage.txt')
        X_sage = np.loadtxt('X_sage.txt')

        W = np.loadtxt('W.txt')
        Theta = np.loadtxt('Theta.txt')
        M = np.loadtxt('M.txt')
        Q = np.loadtxt('Q.txt')
        X = [i * 25 / (len(W)) for i in range(len(W))]


        _, axes = plt.subplots(nrows=4, ncols=1, figsize=(8, 12)) 

        for i, (graphic, name) in enumerate(zip([[-Q, Q_sage], [-M, M_sage], [-Theta, Theta_sage], [-W, W_sage]], ['Q', 'M', '$\Theta$', 'W'])):
            ax = axes[i]
            graphic_cpp, graphic_sage = graphic[0], graphic[1]
            
            ax.axhline(y=0, color='black', linewidth=0.8)
    
            ax.plot([12, 21], [0]*2, 'fuchsia', label='q12-q21', linewidth=4)      # Распределенная нагрузка            
            ax.plot([2, 7, 10, 15], [0]*4, color = 'r', marker = '^', linestyle='', label='R')   # Шарниры
            ax.plot([5], [0], color = 'k', marker = 'o', label='M5 ')              # Момент 
            ax.plot([0], [0], color = 'g', marker = 'D', label='P0 ')              # Сила
            ax.plot([23], [0], color = 'y', marker = 's', label='K23')             # Упругое закрепление

            ax.plot(X, graphic_cpp, 'b', label=name + " MFE") # Решение МКЭ
            ax.plot(X_sage, graphic_sage, linestyle='--', color='peru', label=name + " SAGE") # Решение Sage
         
            ax.set_xlabel('x')
            ax.set_ylabel(name)
            ax.set_title(f'{name}(x)', fontsize=12)
            ax.legend(loc="upper right")
            ax.grid(True)

        plt.tight_layout() 
        plt.show()
