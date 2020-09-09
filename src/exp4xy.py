import parallelstrategies as ps
import numpy as np
import matplotlib.pyplot as plt

def dist(w, p):
    min = 999
    for spot_xy in p.spots_xy:
        d = np.sqrt((w.x-spot_xy[0])**2 + (w.y-spot_xy[1])**2)
        if d < min:
            min = d
    return min

def runDual(duration_mins):
    w = ps.Worm(dt,duration_steps,pir_window)
    w.x = 0.0
    p = ps.Plate()
    p.setGrid()
    p.noise_sd = 0.008
    x_hist = np.zeros(duration_steps)
    y_hist = np.zeros(duration_steps)
    for step in range(duration_steps):
        w.step_rc()
        #w.step_wv(p.conc(w.sx[0],w.sy[0])-p.conc(w.sx[1],w.sy[1]))
        #w.step_pr(step,p.conc(w.x,w.y))
        w.step()
        p.step(dt)
        if np.sqrt(w.x**2 + w.y**2) > p.petri_radius:
            w.bounce(p.petri_radius)
        x_hist[step] = w.x
        y_hist[step] = w.y
    return x_hist,y_hist

reps = 20
pir_window = 0.6
duration_mins = 20
dt = 0.1
duration_steps = int((duration_mins * 60) / dt)
x = np.zeros((reps,duration_steps))
y = np.zeros((reps,duration_steps))
ci = np.zeros(reps)
for i in range(reps):
    x[i],y[i] = runDual(duration_mins)
np.savetxt('exp4xy_null_x.csv', x, delimiter=',')
np.savetxt('exp4xy_null_y.csv', y, delimiter=',')
