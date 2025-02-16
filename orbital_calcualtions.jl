println("Starting simulation...")

#Givens
mu_G = 7.802e-5
mu_E = 2.523e-5

delta_t = 0.0001
v_x = 0.0
v_y = 0.0
x = 0.9992
y = 0.0

#Calculated through L1 Calculations
x_E = 1.020449653771479 

#Partial derivative of the potential with respect to x
function omega_derivative_x(x, y, mu)
    return x-(((1-mu)*(x+mu))/((mu+x)^2+y^2)^(3/2))-((mu*(x-mu-1))/(y^2+(-1-mu+x)^2)^(3/2))
end

#Partial derivative of the potential with respect to y
function omega_derivative_y(x, y, mu)
    return y-((y*(1-mu))/((mu+x)^2+y^2)^(3/2))-((mu*y)/(y^2+(-1-mu+x)^2)^(3/2))
end

#Runge-Kutta method to find velocity and position over time. 
function runge_kutta(v_x, v_y, x, y, mu, delta_t)
    k1x = omega_derivative_x(x, y, mu) + 2 * v_y
    k1y = omega_derivative_y(x, y, mu) - 2 * v_x

    k2x = omega_derivative_x(x + 0.5 * delta_t * v_x, y + 0.5 * delta_t * v_y, mu) + 2 * (v_y + 0.5 * delta_t * k1y)
    k2y = omega_derivative_y(x + 0.5 * delta_t * v_x, y + 0.5 * delta_t * v_y, mu) - 2 * (v_x + 0.5 * delta_t * k1x)

    k3x = omega_derivative_x(x + 0.5 * delta_t * (v_x + 0.5 * delta_t * k2x), y + 0.5 * delta_t * (v_y + 0.5 * delta_t * k2y), mu) + 2 * (v_y + 0.5 * delta_t * k2y)
    k3y = omega_derivative_y(x + 0.5 * delta_t * (v_x + 0.5 * delta_t * k2x), y + 0.5 * delta_t * (v_y + 0.5 * delta_t * k2y), mu) - 2 * (v_x + 0.5 * delta_t * k2x)

    k4x = omega_derivative_x(x + delta_t * (v_x + delta_t * k3x), y + delta_t * (v_y + delta_t * k3y), mu) + 2 * (v_y + delta_t * k3y)
    k4y = omega_derivative_y(x + delta_t * (v_x + delta_t * k3x), y + delta_t * (v_y + delta_t * k3y), mu) - 2 * (v_x + delta_t * k3x)

    v_x += (k1x + 2 * k2x + 2 * k3x + k4x) * delta_t / 6
    v_y += (k1y + 2 * k2y + 2 * k3y + k4y) * delta_t / 6

    x += v_x * delta_t
    y += v_y * delta_t

    return v_x, v_y, x, y
end

for i in 1:300
    global v_x, v_y, x, y = runge_kutta(v_x, v_y, x, y, mu_G, delta_t)
    println("$i, $v_x, $v_y, $x, $y")
    if abs(x - x_E) < 1e-5
        println("Crossed Lagrange Point x_E at iteration $i")
        break
    end
end
