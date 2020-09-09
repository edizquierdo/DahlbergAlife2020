import numpy as np

class Plate():

    def __init__(self):
        self.spots_xy = [(0,0), (0,0)]  #[(0,22), (0,22)]            # list of coordinates of NaCl spots as tuples: (x,y) (mm)
        self.volumes = np.array([5, 5])             # Volume of spots (uL)
        self.concentrations = np.array([0.5, 0.5])  # NaCl concentration of spotting solution (M)
        self.spots_No = self.volumes * self.concentrations    # initial amount of NaCl spotted - 2.5 for the radial assay, or 0.2 for the grid assay (um)
        self.spots_lt = np.array([19.0, 4.0]) * 60 * 60   # list of spot lead times (hours)
        self.scale = 1.0
        self.D = 0.000015    #diffusion constant of NaCl (cm^2/s)
        self.d = 0.18    #thickness of the plate (cm)
        self.petri_radius = 45    # mm
        self.CI_radius = np.sqrt(2/np.pi)*10  #chemotaxis index radius according to Iino and Yoshida 2009 with correction
        self.CI_radial_radius = np.sqrt((45**2)/2)
        self.noise_sd = 0.0    #noise in detected c, as a %. 0% (off) by default.
        self.t = 0.0

    def setGrid(self):
        self.spots_xy = [(-30,10),(-30,-10),(-10,30),(-10,10),(-10,-10),(-10,-30),(10,30),(10,10),(10,-10),(10,-30),(30,10),(30,-10)]
        self.volumes = np.full(12, 1)
        self.concentrations = np.full(12, 0.200)
        self.spots_No = self.volumes * self.concentrations
        self.spots_lt = np.full(12,1.0) * 60 * 60

    def step(self, dt):
        self.t += dt

    def conc(self, x, y):
        c = 0
        for i in range(len(self.spots_xy)):
            t = self.t + self.spots_lt[i]     #age of spot is simulation time (sec) plus lead time ( converting hours to sec)
            dist = np.sqrt((x-self.spots_xy[i][0])**2 + (y-self.spots_xy[i][1])**2)
            No = self.spots_No[i]
            exponent = - dist**2 / (400 * self.D * t)    #calculate the exponent of e
            dividend = No * np.exp(exponent)    #given exponent, calculate the dividend of the equation (including No)
            divisor = 4 * np.pi * self.d * self.D * t    #calculate divisor
            c += (dividend / divisor)     #latter part is noise, with sd given as a %
        c = (c * self.scale) + np.random.normal(0.0, self.noise_sd)
        if c < 0.0:
            c = 0.0
        return c

class Worm():

    def __init__(self, dt, duration_steps, pt_window_duration):
        self.dt = dt
        self.v = 0.12                   # Velocity (mm/s)
        self.wv_window = 0.01           # Weathervane sampling window width (mm)
        self.rcr_cycle_duration = 12    # Duration of random curving rate cycle (sec)
        self.rcr_cycle_steps = int(np.rint(self.rcr_cycle_duration/self.dt))    # Number of steps in random curving cycle
        self.pf_a = 0.023               # Pirouette frequency sigmoid parameters
        self.pf_b = 0.40
        self.pf_c = 140
        self.pf_d = 0.0033
        self.wv_scale = 12.7            # Weathervane scale/strength (unitless)
        self.pt_window_duration = pt_window_duration   # Length of pirouette temporal window (sec)
        self.pt_window_steps = int(np.rint(self.pt_window_duration/self.dt))  # Lnegth of pirouette window (steps) #ZZZ
        self.bpf = 0.0197       # Baseline pirouette frequency (events/sec)
        self.pf = self.bpf * self.dt    # Initialized as baseline pirouette frequency, updated by self.update_pf(). Input as event/sec (events/dt)
        self.rcsd_degrees = 32.3        # Random curving standard deviation (deg/mm)
        self.rcsd = np.deg2rad(self.rcsd_degrees) * self.v * self.dt    # Random curving std dev, entered as deg/mm and converted to deg/dt
        self.x = 22.5 #45/np.sqrt(2) #22.5 #0                      # X-position (mm)
        self.y = 0.0                      # Y-position (mm)
        self.h = 0                      # Heading (rad)
        self.sh = np.array([np.pi/2,-np.pi/2]) + self.h         # heading of sensors from body (abs = rel + body) (rad)
        self.sx = (np.cos(self.sh) * self.wv_window/2) + self.x    #x-coordinate of the sensors (mm)
        self.sy = (np.sin(self.sh) * self.wv_window/2) + self.y    #y-coordinate of the sensors (mm)
        self.md = 0                     # Distance to nearest spot (mm)
        self.wv = 0                     # Current weathervane curving rate bias (rad/dt)
        self.rcr = 0                    # Current random curving rate (rad/dt)
        self.drcr = 0                   # Stepwise change in random curving rate (rad/dt^2)
        self.rcr_clock = self.rcr_cycle_steps    # Random curving cycle clock
        self.c_hist = np.zeros(duration_steps)

    def step_rc(self):               # RANDOM CURVING
        self.rcr += self.drcr    #modify the curving rate by the stepwise curving rate change: (orcr-nrcr) / (cycle length)
        if self.rcr_clock >= self.rcr_cycle_steps:    # Check if cr needs updating
            orcr = self.rcr    #store old cr "guidepost'
            nrcr = np.random.normal(0, self.rcsd)    #find new cr "guidepost"
            self.drcr = (nrcr - orcr) / self.rcr_cycle_steps    #create new cycle's cr
            self.rcr_clock = 0    #reset clock for new cycle
        self.h += self.rcr    #add curving rate to heading
        self.rcr_clock += 1    #update curving rate clock

    def step_wv(self, cdiff):               # WEATHERVANE
         self.wv = np.deg2rad( self.wv_scale * (cdiff) * 10/self.wv_window) * self.dt * self.v    #10/self.wv_window is the scale factor needed to yield 1 cm
         self.h += self.wv    #add wvcrb to heading

    def step_pr(self, step, c):               # PIROUETTE ZZZZ
        self.c_hist[step] = c    #add current concentration to history
        start = max([0, step-self.pt_window_steps])  # Floor the start value to 0 so that we don't go negative and index from the end of the list
        end = max([2, step+1]) # Floor the end value to 1 so that in the beginning of the run we don't get an empty list (and self.pf=NaN)
        ch = self.c_hist[start:end]  # Get the recent several c from self.c_hist
        dc = np.diff(ch)
        adc = np.mean(dc)  # Getting average diff, reason being so we can use IY2009 sensorimotor transformation
        dcdt = adc/self.dt
        self.pf = (self.pf_a/(self.pf_b + np.exp(self.pf_c * dcdt)) + self.pf_d) * self.dt  #based on IY2009 (sigmoidal fit on empirical data) (event/dt)

    def step_pr_v2(self, step, c):               # PIROUETTE ZZZZ
        self.c_hist[step] = c    #add current concentration to history
        start = max([0, step-self.pt_window_steps])  # Floor the start value to 0 so that we don't go negative and index from the end of the list
        past_c = self.c_hist[start]  # Get the recent several c from self.c_hist
        dcdt = (c - past_c)/(self.dt*self.pt_window_steps)
        self.pf = (self.pf_a/(self.pf_b + np.exp(self.pf_c * dcdt)) + self.pf_d) * self.dt  #based on IY2009 (sigmoidal fit on empirical data) (event/dt)

    def step(self):
        if np.random.random() < self.pf:    #determine if a pirouette event happens
            self.h = np.random.random() * 2  * np.pi    #do a pirouette
        self.x += np.cos(self.h) * self.v * self.dt
        self.y += np.sin(self.h) * self.v * self.dt
        self.sh = np.array([np.pi/2,-np.pi/2]) + self.h         # heading of sensors from body (abs = rel + body) (rad)
        self.sx = (np.cos(self.sh) * self.wv_window/2) + self.x    #x-coordinate of the sensors (mm)
        self.sy = (np.sin(self.sh) * self.wv_window/2) + self.y    #y-coordinate of the sensors (mm)

    def bounce(self, distance):
        bearing = np.arctan2(self.y, self.x)
        self.x = np.cos(bearing) * distance
        self.y = np.sin(bearing) * distance
        self.h = np.random.random() * 2  * np.pi
