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

def runDual(steepness,duration_mins):
    w = ps.Worm(dt,duration_steps,pir_window)
    p = ps.Plate()
    p.setGrid()
    p.scale = steepness
    #p.noise_sd = noise
    d_hist = np.zeros(duration_steps)
    inside = 0
    for step in range(duration_steps):
        w.step_rc()
        #w.step_wv(p.conc(w.sx[0],w.sy[0])-p.conc(w.sx[1],w.sy[1]))
        w.step_pr(step,p.conc(w.x,w.y)) # v1.0
        #w.step_pr_v2(step,p.conc(w.x,w.y)) # v2.0
        w.step()
        p.step(dt)
        if np.sqrt(w.x**2 + w.y**2) > p.petri_radius:
            w.bounce(p.petri_radius)
        d_hist[step] = dist(w,p)
        if d_hist[step] < p.CI_radius:
            inside += 1
    outside = duration_steps - inside
    ci = (inside - outside)/duration_steps
    n=int(duration_mins/2)
    return np.mean(d_hist[-n:]),np.min(d_hist),ci

reps = 5000
conds = 11
steepness = np.linspace(0.1,1.9,conds)
pir_window = 0.6
noise = 0.008
duration_mins = 20
dt = 0.1
duration_steps = int((duration_mins * 60) / dt)
d = np.zeros((reps,conds))
m = np.zeros((reps,conds))
ci = np.zeros((reps,conds))
for i in range(reps):
    print("Rep:",i)
    for j in range(conds):
        #print("Steepness:",steepness[j])
        d[i][j],m[i][j],ci[i][j] = runDual(steepness[j],duration_mins)
np.savetxt('exp3b_pr_d.csv', d, delimiter=',')
np.savetxt('exp3b_pr_m.csv', m, delimiter=',')
np.savetxt('exp3b_pr_ci.csv', ci, delimiter=',')
