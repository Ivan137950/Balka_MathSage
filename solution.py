import s 

params = {'p': 100., 
          'm': 700., 
          'k': 0.01, 
          'e': 2000., 
          'j': 9., 
          'q1': 10., 
          'q2': 20.
          }   

    
# Решение SAGE
Sage_sol = s.SageSolution(params)
Sage_sol.run_sage_solution()

# Решение CPP
solution_cpp = s.CppSolution(params)
solution_cpp.run_cpp_solution()

#Графики
g = s.Graphic()
g.draw_graphics()
