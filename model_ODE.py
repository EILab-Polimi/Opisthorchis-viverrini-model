import numpy as np
from scipy.integrate import solve_ivp

def model_ODE(Time, par, setup, y0):
    

    def odemodel(p, nNodes, H, S, F, A, W, chi, xi, epsilon, theta, tspan, y0_mat):
        
        y0 = y0_mat.flatten(order='F')  # equivalente a MATLAB column-major vectorization
       
        def eqs(t, y):
           
            dy = np.zeros(4 * nNodes)

            

            
            with np.errstate(divide='ignore', invalid='ignore'):
                P4 = W @ (F * y[3::4]) / F - np.sum(W, axis=1) * y[3::4]
            P4[np.isnan(P4)] = 0 
           

            beta_E = p['beta_E']
            theta_C = theta

            
            dy[0::4] = chi * F / H * epsilon * y[3::4] - (p['mu_W'] + p['mu_H']) * y[0::4]

            # dy(2:4:end)
            dy[1::4] = (
                xi * p['rho_E'] * H * y[0::4] / (p['alpha'] + y[0::4]) / A
                - (p['mu_E'] + beta_E * S / A) * y[1::4]
            )

            # dy(3:4:end)
            dy[2::4] = beta_E * (1 - y[2::4]) * y[1::4] - p['mu_S'] * y[2::4]

            # dy(4:4:end)
            dy[3::4] = theta_C * y[2::4] + P4 - (p['mu_F'] + chi) * y[3::4]

            dy[np.isnan(dy)] = 0
            return dy

        
        sol = solve_ivp(eqs, (tspan[0], tspan[-1]), y0, t_eval=tspan, method='RK45')
        y = sol.y.T  
        return y

    
    y = odemodel(
        par,
        setup['nNodes'],
        setup['H'],
        setup['S'],
        setup['F'],
        setup['A'],
        setup['W'],
        setup['chi'],
        setup['par']['xi'],
        setup['par']['epsilon'],
        setup['par']['theta'],
        Time,
        y0,
    )

    # MATLAB: y(:,2:4:end) = [];
    n_cols = y.shape[1]
    print(n_cols)
    keep_mask = np.ones(n_cols, dtype=bool)
    keep_mask[1::4] = False
    y = y[:, keep_mask]

    return y
