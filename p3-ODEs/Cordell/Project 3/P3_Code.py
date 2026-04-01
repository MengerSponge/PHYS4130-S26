#############################################################

# going to pull everything from the nb and organize it like this

#############################################################

# 1) Select SHO or SHO with damping
# 2) Select numerical method 
# 3) Plot phase space
# 4) Plot Energy vs Time
# 5) Recurse or end program

#############################################################


from Functions import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.integrate import odeint

print("Select system:")
print("1. Simple Harmonic Oscillator (SHO)")
print("2. Damped Harmonic Oscillator")
system_choice = int(input("Enter your choice (1 or 2): "))

##Fix RK4

print("\nSelect numerical method:")
print("1. Euler's Method")
print("2. RK2")
print("3. RK4")
print("4. Verlet Integration")
print("5. Scipy's ODEINT")
method_choice = int(input("Enter your choice (1, 2, 3, 4, or 5): "))


if system_choice == 1:
    x0 = float(input("Enter initial position: "))
    v0 = float(input("Enter initial velocity: "))
    tmin = 0
    tmax = float(input("Enter final time: "))
    nts = int(input("Enter number of time steps (Cannot control nts for RK4 because Scipy chooses it): "))
    if method_choice == 1:
        t, x, v = SHO_solver_Euler(x0, v0, tmin, tmax, nts, SHO_deriv)
    elif method_choice == 2:
        t, x, v = SHO_solver_RK2(x0, v0, tmin, tmax, nts, SHO_deriv)
    elif method_choice == 3:
        t, x, v = solve_ivp(fun_SHO, (tmin, tmax), [x0, v0], t_eval=np.linspace(tmin, tmax, nts, endpoint=False), method='RK45').t, solve_ivp(fun_SHO, (tmin, tmax), [x0, v0], t_eval=np.linspace(tmin, tmax, nts, endpoint=False), method='RK45').y[0], solve_ivp(fun_SHO, (tmin, tmax), [x0, v0], t_eval=np.linspace(tmin, tmax, nts, endpoint=False), method='RK45').y[1]
    elif method_choice == 4:
        t, x, v = verlet_solver(x0, v0, tmin, tmax, nts, A_verlet_SHO)
    elif method_choice == 5:
        t, x, v = SHO_solver_ODEINT(x0, v0, tmin, tmax, nts, SHO_deriv)
        pass
    else:
        print("Invalid method choice for simple harmonic oscillator.")

if system_choice == 2:
    x0 = float(input("Enter initial position: "))
    v0 = float(input("Enter initial velocity: "))
    tmin = 0
    tmax = float(input("Enter final time: "))
    nts = int(input("Enter number of time steps: "))
    if method_choice == 1:
        t, x, v = SHO_solver_Euler(x0, v0, tmin, tmax, nts, SHO_deriv_damped)
    elif method_choice == 2:
        t, x, v = SHO_solver_RK2(x0, v0, tmin, tmax, nts, SHO_deriv_damped)
    elif method_choice == 3:
        t, x, v = solve_ivp(fun, (tmin, tmax), [x0, v0], t_eval=np.linspace(tmin, tmax, nts, endpoint=False), method='RK45').t, solve_ivp(fun, (tmin, tmax), [x0, v0], t_eval=np.linspace(tmin, tmax, nts, endpoint=False), method='RK45').y[0], solve_ivp(fun, (tmin, tmax), [x0, v0], t_eval=np.linspace(tmin, tmax, nts, endpoint=False), method='RK45').y[1]
    elif method_choice == 4:
        t, x, v = verlet_solver(x0, v0, tmin, tmax, nts, A_verlet_damped)
    elif method_choice == 5:
        t, x, v = SHO_solver_ODEINT(x0, v0, tmin, tmax, nts, SHO_deriv_damped)
    else:
        print("Invalid method choice for damped harmonic oscillator.")

print("\nSelect display:")
print("1. Phase Space Plot")
print("2. Energy vs Time Plot")
print("3. Error at some time")
print("4. Number of time steps to reach a target error")
print("5. loglog of error vs time steps")
plot_choice = int(input("Enter your choice (1, 2, 3, 4, or 5): "))

if plot_choice == 1:
    plt.plot(x, v)
    plt.xlabel("Position")
    plt.ylabel("Velocity")
    if system_choice == 1:
        plt.title("Phase Space of Simple Harmonic Oscillator")
    else:
        plt.title("Phase Space of Damped Harmonic Oscillator")
    plt.grid()
    plt.show()

elif plot_choice == 2:
    energy = H(x, v)
    plt.plot(t, energy)
    plt.xlabel("Time")
    plt.ylabel("Energy")
    if system_choice == 1:
        plt.title("Energy vs Time for Simple Harmonic Oscillator")
    else:
        plt.title("Energy vs Time for Damped Harmonic Oscillator")
    plt.grid()
    plt.show()



########################################################################3

# Error at some time

elif plot_choice == 3:
    tmax_error = float(input("Enter time at which to calculate error: "))
    if system_choice == 1:
        t_verlet, x_verlet, v_verlet = verlet_solver(x0, v0, tmin, tmax, nts, A_verlet_SHO)

        t_eval = np.linspace(tmin, tmax, nts, endpoint=False)
        sol = solve_ivp(fun, (tmin, tmax), [x0, v0], t_eval=t_eval, method='RK45')
        t_rk4 = sol.t
        x_rk4 = sol.y[0]
        v_rk4 = sol.y[1]

        t_rk2, x_rk2, v_rk2 = SHO_solver_RK2(x0, v0, tmin, tmax, nts, SHO_deriv)
        t_odeint, x_odeint, v_odeint = SHO_solver_ODEINT(x0, v0, tmin, tmax, nts, SHO_deriv)
        t_euler, x_euler, v_euler = SHO_solver_Euler(x0, v0, tmin, tmax, nts, SHO_deriv)
        
        relative_error_verlet_at_tmax = relative_error(x_verlet[-1], analytical_SHO(x0, v0, t_verlet[-1]))
        relative_error_rk4_at_tmax = relative_error(x_rk4[-1], analytical_SHO(x0, v0, t_rk4[-1]))
        relative_error_rk2_at_tmax = relative_error(x_rk2[-1], analytical_SHO(x0, v0, t_rk2[-1]))
        relative_error_odeint_at_tmax = relative_error(x_odeint[-1], analytical_SHO(x0, v0, t_odeint[-1]))
        relative_error_euler_at_tmax = relative_error(x_euler[-1], analytical_SHO(x0, v0, t_euler[-1]))
        print(f"Relative error at t={tmax} for Verlet: {relative_error_verlet_at_tmax:.6e}")
        print(f"Relative error at t={tmax} for RK4: {relative_error_rk4_at_tmax:.6e}")
        print(f"Relative error at t={tmax} for RK2: {relative_error_rk2_at_tmax:.6e}")
        print(f"Relative error at t={tmax} for ODEINT: {relative_error_odeint_at_tmax:.6e}")
        print(f"Relative error at t={tmax} for Euler: {relative_error_euler_at_tmax:.6e}")
    elif system_choice == 2:
        t_verlet, x_verlet, v_verlet = verlet_solver(x0, v0, tmin, tmax, nts, A_verlet_damped)

        t_eval = np.linspace(tmin, tmax, nts, endpoint=False)
        sol = solve_ivp(fun, (tmin, tmax), [x0, v0], t_eval=t_eval, method='RK45')
        t_rk4 = sol.t
        x_rk4 = sol.y[0]
        v_rk4 = sol.y[1]

        t_rk2, x_rk2, v_rk2 = SHO_solver_RK2(x0, v0, tmin, tmax, nts, SHO_deriv_damped)
        t_odeint, x_odeint, v_odeint = SHO_solver_ODEINT(x0, v0, tmin, tmax, nts, SHO_deriv_damped)
        t_euler, x_euler, v_euler = SHO_solver_Euler(x0, v0, tmin, tmax, nts, SHO_deriv_damped)
        
        relative_error_verlet_at_tmax = relative_error(x_verlet[-1], analytical_damped_SHO(x0, v0, t_verlet[-1]))
        relative_error_rk4_at_tmax = relative_error(x_rk4[-1], analytical_damped_SHO(x0, v0, t_rk4[-1]))
        relative_error_rk2_at_tmax = relative_error(x_rk2[-1], analytical_damped_SHO(x0, v0, t_rk2[-1]))
        relative_error_odeint_at_tmax = relative_error(x_odeint[-1], analytical_damped_SHO(x0, v0, t_odeint[-1]))
        relative_error_euler_at_tmax = relative_error(x_euler[-1], analytical_damped_SHO(x0, v0, t_euler[-1]))
        print(f"Relative error at t={tmax} for Verlet: {relative_error_verlet_at_tmax:.6e}")
        print(f"Relative error at t={tmax} for RK4: {relative_error_rk4_at_tmax:.6e}")
        print(f"Relative error at t={tmax} for RK2: {relative_error_rk2_at_tmax:.6e}")
        print(f"Relative error at t={tmax} for ODEINT: {relative_error_odeint_at_tmax:.6e}")
        print(f"Relative error at t={tmax} for Euler: {relative_error_euler_at_tmax:.6e}")



################################################################3

# nts until target error


elif plot_choice == 4:
    target_error = float(input("Enter target relative error: "))
    if system_choice == 1:
        analytical_solution = analytical_SHO
        nts = [2**k for k in range(1, 20)]  # powers of 2
        for nts in nts:
            t_verlet, x_verlet, v_verlet = verlet_solver(x0, v0, tmin, tmax, nts, A_verlet_SHO)
            relative_error_verlet_at_tmax = relative_error(x_verlet[-1], analytical_solution(x0, v0, t_verlet[-1]))
            if relative_error_verlet_at_tmax < target_error:
                print(f"Verlet method achieves relative error < {target_error} at t={tmax} with nts={nts}")
                break
        nts = [2**k for k in range(1, 20)]  # powers of 2
        for nts in nts:
            t_rk2, x_rk2, v_rk2 = SHO_solver_RK2(x0, v0, tmin, tmax, nts, SHO_deriv)
            relative_error_rk2_at_tmax = relative_error(x_rk2[-1], analytical_solution(x0, v0, t_rk2[-1]))
            if relative_error_rk2_at_tmax < target_error:
                print(f"RK2 method achieves relative error < {target_error} at t={tmax} with nts={nts}")
                break
        nts = [2**k for k in range(1, 20)]  # powers of 2
        for nts in nts:
            t_euler, x_euler, v_euler = SHO_solver_Euler(x0, v0, tmin, tmax, nts, SHO_deriv)
            relative_error_euler_at_tmax = relative_error(x_euler[-1], analytical_solution(x0, v0, t_euler[-1]))
            if relative_error_euler_at_tmax < target_error:
                print(f"Euler method achieves relative error < {target_error} at t={tmax} with nts={nts}")
                break
        
    elif system_choice == 2:
        analytical_solution = analytical_damped_SHO
        nts = [2**k for k in range(1, 20)]  # powers of 2
        for nts in nts:
            t_verlet, x_verlet, v_verlet = verlet_solver(x0, v0, tmin, tmax, nts, A_verlet_damped)
            relative_error_verlet_at_tmax = relative_error(x_verlet[-1], analytical_solution(x0, v0, t_verlet[-1]))
            if relative_error_verlet_at_tmax < target_error:
                print(f"Verlet method achieves relative error < {target_error} at t={tmax} with nts={nts}")
                break

        nts = [2**k for k in range(1, 20)]  # powers of 2
        for nts in nts:
            t_rk2, x_rk2, v_rk2 = SHO_solver_RK2(x0, v0, tmin, tmax, nts, SHO_deriv_damped)
            relative_error_rk2_at_tmax = relative_error(x_rk2[-1], analytical_solution(x0, v0, t_rk2[-1]))
            if relative_error_rk2_at_tmax < target_error:
                print(f"RK2 method achieves relative error < {target_error} at t={tmax} with nts={nts}")
                break

        nts = [2**k for k in range(1, 20)]  # powers of 2
        for nts in nts:
            t_euler, x_euler, v_euler = SHO_solver_Euler(x0, v0, tmin, tmax, nts, SHO_deriv_damped)
            relative_error_euler_at_tmax = relative_error(x_euler[-1], analytical_solution(x0, v0, t_euler[-1]))
            if relative_error_euler_at_tmax < target_error:
                print(f"Euler method achieves relative error < {target_error} at t={tmax} with nts={nts}")
                break

######################################################################33
# loglog plots

elif plot_choice == 5:
    nts_list = [2**k for k in range(4, 18)]  # powers of 2

    relative_errors = np.zeros(len(nts_list))
    if system_choice == 1:
        if method_choice == 1:
            for nts in nts_list:
                t_euler, x_euler, v_euler = SHO_solver_Euler(x0, v0, tmin, tmax, nts, SHO_deriv)
                err = relative_error(x_euler[-1], analytical_SHO(x0, v0, t_euler[-1]))
                relative_errors[nts_list.index(nts)] = err

            plt.loglog(nts_list, relative_errors, 'o-', label="Euler Method")
            plt.xlabel('Number of Time Steps')
            plt.ylabel('Relative Error at tmax')
            plt.title('Relative Error at tmax vs Number of Time Steps for Euler Method')
            plt.legend()
            plt.show()

            # print slope of line
            log_nts = np.log(nts_list)
            log_errors = np.log(relative_errors)
            slope = (log_errors[-1] - log_errors[0]) / (log_nts[-1] - log_nts[0])
            print(f"Slope of log-log plot: {slope:.2f}")

        elif method_choice == 2:
            for nts in nts_list:
                t_rk2, x_rk2, v_rk2 = SHO_solver_RK2(x0, v0, tmin, tmax, nts, SHO_deriv)
                err = relative_error(x_rk2[-1], analytical_SHO(x0, v0, t_rk2[-1]))
                relative_errors[nts_list.index(nts)] = err

            plt.loglog(nts_list, relative_errors, 'o-', label="RK2 Method")
            plt.xlabel('Number of Time Steps')
            plt.ylabel('Relative Error at tmax')
            plt.title('Relative Error at tmax vs Number of Time Steps for RK2 Method')
            plt.legend()
            plt.show()

            # print slope of line
            log_nts = np.log(nts_list)
            log_errors = np.log(relative_errors)
            slope = (log_errors[-1] - log_errors[0]) / (log_nts[-1] - log_nts[0])
            print(f"Slope of log-log plot: {slope:.2f}")

        elif method_choice == 3:
            print(f"loglog plot unavailable for RK4 since Scipy chooses nts")
        
        elif method_choice==4: 

            for nts in nts_list:
                t_verlet, x_verlet, v_verlet = verlet_solver(x0, v0, tmin, tmax, nts, A_verlet_SHO)
                err = relative_error(x_verlet[-1], analytical_SHO(x0, v0, t_verlet[-1]))
                relative_errors[nts_list.index(nts)] = err

            plt.loglog(nts_list, relative_errors, 'o-', label="Verlet Method")
            plt.xlabel('Number of Time Steps')
            plt.ylabel('Relative Error at tmax')
            plt.title('Relative Error at tmax vs Number of Time Steps for Verlet Method')
            plt.legend()
            plt.show()

            # print slope of line
            log_nts = np.log(nts_list)
            log_errors = np.log(relative_errors)
            slope = (log_errors[-1] - log_errors[0]) / (log_nts[-1] - log_nts[0])
            print(f"Slope of log-log plot: {slope:.2f}")

        elif method_choice==5: 
            print(f"loglog plot unavailable for ODEINT")
    elif system_choice == 2:
        if method_choice == 1:
            for nts in nts_list:
                t_euler, x_euler, v_euler = Euler_solver(x0, v0, tmin, tmax, nts, SHO_deriv_damped)
                err = relative_error(x_euler[-1], analytical_damped_SHO(x0, v0, t_euler[-1]))
                relative_errors[nts_list.index(nts)] = err

            plt.loglog(nts_list, relative_errors, 'o-', label="Euler Method")
            plt.xlabel('Number of Time Steps')
            plt.ylabel('Relative Error at tmax')
            plt.title('Relative Error at tmax vs Number of Time Steps for Euler Method')
            plt.legend()
            plt.show()

            # print slope of line
            log_nts = np.log(nts_list)
            log_errors = np.log(relative_errors)
            slope = (log_errors[-1] - log_errors[0]) / (log_nts[-1] - log_nts[0])
            print(f"Slope of log-log plot: {slope:.2f}")
        elif method_choice == 2:
            for nts in nts_list:
                t_rk2, x_rk2, v_rk2 = SHO_solver_RK2(x0, v0, tmin, tmax, nts, SHO_deriv_damped)
                err = relative_error(x_rk2[-1], analytical_damped_SHO(x0, v0, t_rk2[-1]))
                relative_errors[nts_list.index(nts)] = err

            plt.loglog(nts_list, relative_errors, 'o-', label="RK2 Method")
            plt.xlabel('Number of Time Steps')
            plt.ylabel('Relative Error at tmax')
            plt.title('Relative Error at tmax vs Number of Time Steps for RK2 Method')
            plt.legend()
            plt.show()

            # print slope of line
            log_nts = np.log(nts_list)
            log_errors = np.log(relative_errors)
            slope = (log_errors[-1] - log_errors[0]) / (log_nts[-1] - log_nts[0])
            print(f"Slope of log-log plot: {slope:.2f}")
        elif method_choice == 3:
            print(f"loglog plot unavailable for RK4 since Scipy chooses nts")
        elif method_choice==4:
            for nts in nts_list:
                t_verlet, x_verlet, v_verlet = verlet_solver(x0, v0, tmin, tmax, nts, A_verlet_damped)
                err = relative_error(x_verlet[-1], analytical_damped_SHO(x0, v0, t_verlet[-1]))
                relative_errors[nts_list.index(nts)] = err

            plt.loglog(nts_list, relative_errors, 'o-', label="Verlet Method")
            plt.xlabel('Number of Time Steps')
            plt.ylabel('Relative Error at tmax')
            plt.title('Relative Error at tmax vs Number of Time Steps for Verlet Method')
            plt.legend()
            plt.show()

            # print slope of line
            log_nts = np.log(nts_list)
            log_errors = np.log(relative_errors)
            slope = (log_errors[-1] - log_errors[0]) / (log_nts[-1] - log_nts[0])
            print(f"Slope of log-log plot: {slope:.2f}")
        elif method_choice==5:
            print(f"loglog plot unavailable for ODEINT")


