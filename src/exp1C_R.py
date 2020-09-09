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
    p = ps.Plate()
    d_hist = np.zeros(duration_steps)
    inside = 0
    for step in range(duration_steps):
        w.step_rc()
        w.step_wv(p.conc(w.sx[0],w.sy[0])-p.conc(w.sx[1],w.sy[1]))
        w.step_pr(step,p.conc(w.x,w.y))
        w.step()
        p.step(dt)
        if np.sqrt(w.x**2 + w.y**2) > p.petri_radius:
            w.bounce(p.petri_radius)
        d_hist[step] = dist(w,p)
        if d_hist[step] < p.CI_radial_radius:
            inside += 1
    outside = duration_steps - inside
    ci = (inside - outside)/duration_steps
    return d_hist,ci

reps = 5000
pir_window = 0.6
duration_mins = 20
dt = 0.6
duration_steps = int((duration_mins * 60) / dt)
d = np.zeros((reps,duration_steps))
ci = np.zeros(reps)
for i in range(reps):
    d[i],ci[i] = runDual(duration_mins)
np.savetxt('exp1A2_R_d.csv', d, delimiter=',')
np.savetxt('exp1A2_R_ci.csv', ci, delimiter=',')
