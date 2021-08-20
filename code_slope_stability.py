#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

def solve_quadratic(a, b, c, xm):
    discriminant = b ** 2 - 4 * a * c
    r1 = (-b + np.sqrt(discriminant))/ (2 * a) 
    r2 = (-b - np.sqrt(discriminant))/ (2 * a)
    if xm == 0:
        return r2
    if r1 < xm:
        return r1
    elif r2 < xm:
        return r2

pi = np.pi
cohesion = float(input('cohesion (in kN/m2): '))
friction_ang = int(input('friction angle (deg.): '))
friction_ang *= pi/ 180
unit_wt = float(input('unit weight (kN/m3): '))

H = float(input('height of slope (metres): '))
D = float(input('depth factor for firm base: ')) # depth factor for firm base
beta = float(input('inclination of slope (deg.): ')) # inclination of slope (deg.)
beta *= pi/ 180
num_points = int(input('number of chords: ')) # number of points on the slope-top to be considered for analysis
                                              # for each of these points, slip surfaces of varying curvature will be analysed
#factor = (D * np.tan(beta) + 1) ** (1/ (num_points - 1)) #float(input('value slightly greater than 1: '))
factor = 1.2
data = [[], [], []]

plt.plot([0, H/ np.tan(beta)], [0, H], 'k')
plt.plot([H/ np.tan(beta), factor ** num_points * H/ np.tan(beta)], 2 * [H], 'k')
plt.plot([0, factor ** num_points * H/ np.tan(beta)], 2 * [H * (1 - D)], 'k')

num_slip_surfaces = int(input('number of slip surfaces: ')) 
num_slices = int(input('number of slices: '))
for i in range(num_points):
    xi = factor ** i * H/ np.tan(beta)
    yi = H
    xi_arr = np.linspace(0, xi, 1000)
    
    if xi <= H/ np.tan(beta) + D * H:
        R_min = (xi ** 2 + yi ** 2)/ (2 * xi) # true if {xi - H/ tan(beta)} <= D * H
    else:
        R_min = solve_quadratic(yi ** 2, (xi ** 2 + yi ** 2) * (yi - 2 * D * H), (0.5 * (xi ** 2 + yi ** 2) - yi * D * H) ** 2 + (xi * D * H) ** 2, 0)
    # distribute the slices on both sides proportionally
    num_slices1 = int((H/ np.tan(beta))/ xi * num_slices)
    num_slices2 = num_slices - num_slices1
    
    for j in range(num_slip_surfaces):
        R = factor ** j * R_min
        # having known radius and end points of the arc, arc's centre can be found
        xc = solve_quadratic(1, -xi, yi ** 2/ 4 * (1 + (xi/ yi) ** 2) - R ** 2/ (1 + (xi/ yi) ** 2), xi/ 2)
        yc = -xi/ yi * (xc - xi/ 2) + yi/ 2

        plt.plot(xi_arr, yc - np.sqrt(abs(R ** 2 - (xi_arr - xc) ** 2)), 'k-', alpha = 0.5)
        plt.plot(xc, yc, 'k.')
        plt.plot([xc, 0], [yc, 0], 'k--', alpha = 0.5)
        plt.plot([xc, xi], [yc, yi], 'k--', alpha = 0.5)
        # equation of circle is now known
        resisting_moments = 0
        driving_moments = 0
        
        del_x = (H/ np.tan(beta))/ num_slices1
        xL = 0
        for k in range(num_slices1):
            xR = xL + del_x
            y_tL = np.tan(beta) * xL
            y_tR = np.tan(beta) * xR
            # if x-coordinate is given, y-coordinate can be computed from the equation of circle 
            y_bL = yc - np.sqrt(abs(R ** 2 - (xL - xc) ** 2)) # the term inside sqrt may be negative (but close to zero), hence abs is used
            y_bR = yc - np.sqrt(abs(R ** 2 - (xR - xc) ** 2))
            #plt.plot([xL, xR, xR, xL, xL], [y_bL, y_bR, y_tR, y_tL, y_bL], 'k-', alpha = 0.5) # visualise slices
            
            # each slice can be approximated as a trapezium - composed of 2 triangles and 1 rectangle
            # its centre of mass can be determined by x_cm = sum(Ai * xi)/ sum(Ai) and y_cm = sum(Ai * yi)/ sum(Ai) 
            x_cm = []
            area = []
            x_cm.append(xL + 2/ 3 * del_x) # x-coordinate of upper triangle of the slice
            area.append(0.5 * del_x * (y_tR - y_tL))
            x_cm.append((xL + xR)/ 2)
            if y_bL >= y_bR:
                area.append(del_x * (y_tL - y_bL)) 
                x_cm.append(xL + 2/ 3 * del_x)
                area.append(0.5 * del_x * (y_bL - y_bR))
            else:
                area.append(del_x * (y_tL - y_bR))
                x_cm.append(xL + 1/ 3 * del_x)
                area.append(0.5 * del_x * (y_bR - y_bL))
            
            sum_area = 0
            sum_area_x = 0
            for l in range(3):
                sum_area += area[l]
                sum_area_x += area[l] * x_cm[l]
            X = sum_area_x/ sum_area
            base_length = np.sqrt(del_x ** 2 + (y_bL - y_bR) ** 2)
            alpha = np.arctan((y_bR - y_bL)/ del_x)
            
            resisting_moments += (cohesion * base_length + unit_wt * sum_area * np.cos(alpha) * np.tan(friction_ang)) * R
            driving_moments += unit_wt * sum_area * (X - xc)
            xL += del_x
        if xi != H/ np.tan(beta):
            del_x = (xi - H/ np.tan(beta))/ num_slices2
            xL = H/ np.tan(beta)
            for k in range(num_slices2):
                xR = xL + del_x
                yT = H
                y_bL = yc - np.sqrt(abs(R ** 2 - (xL - xc) ** 2))
                y_bR = yc - np.sqrt(abs(R ** 2 - (xR - xc) ** 2))
                #plt.plot([xL, xR, xR, xL, xL], [y_bL, y_bR, yT, yT, y_bL], 'k-', alpha = 0.5) # visualize slices

                # each slice can be approximated as a trapezium - composed of 2 triangles and 1 rectangle
                # its centre of mass can be determined by x_cm = sum(Ai * xi)/ sum(Ai) and y_cm = sum(Ai * yi)/ sum(Ai) 
                x_cm = []
                area = []
                x_cm.append((xL + xR)/ 2)
                if y_bL >= y_bR:
                    area.append(del_x * (yT - y_bL)) 
                    x_cm.append(xL + 2/ 3 * del_x)
                    area.append(0.5 * del_x * (y_bL - y_bR))
                else:
                    area.append(del_x * (yT - y_bR))
                    x_cm.append(xL + 1/ 3 * del_x)
                    area.append(0.5 * del_x * (y_bR - y_bL))

                sum_area = 0
                sum_area_x = 0
                for l in range(2):
                    sum_area += area[l]
                    sum_area_x += area[l] * x_cm[l]
                X = sum_area_x/ sum_area
                
                base_length = np.sqrt(del_x ** 2 + (y_bL - y_bR) ** 2)
                alpha = np.arctan((y_bR - y_bL)/ del_x)

                resisting_moments += (cohesion * base_length + unit_wt * sum_area * np.cos(alpha) * np.tan(friction_ang)) * R
                driving_moments += unit_wt * sum_area * (X - xc)
                xL += del_x
        
        FS = resisting_moments/ driving_moments
        data[0].append(xi)
        data[1].append(R)
        data[2].append(FS)

print('Min FS:', min(data[2]))
plt.axis('equal')
#plt.savefig('slope.png', dpi = 200)
plt.show()

fig = go.Figure(data = [go.Scatter3d(x = data[0], y = data[1], z = data[2], mode = 'markers')])
fig.show()


# 
