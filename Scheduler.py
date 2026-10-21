# Dependencies,
import numpy as np

class SatelliteSystem:
    """Class to represent a satellite system with multiple satellites and ground stations."""
    def __init__(self, planetary_radius, planetary_mass, planetary_angular_velocity, G=6.67430e-11):
        self.sats = []
        self.ground_stations = []
        self.planetary_radius = planetary_radius
        self.planetary_mass = planetary_mass
        self.planetary_angular_velocity = planetary_angular_velocity
        self.G = G

    def add_sat(self, sat):
        sat.planetary_radius = self.planetary_radius
        sat.planetary_mass = self.planetary_mass
        sat.planetary_angular_velocity = self.planetary_angular_velocity
        sat.G = self.G
        self.sats.append(sat)

    def add_station(self, station):
        self.ground_stations.append(station)

    def compute_window(self, sat, station, interval, delta_t=1000, format=False):
        """Compute the visibility of a satellite from a ground station at time t."""

        # PARAMETERS
        rE = self.planetary_radius
        h = sat.altitude
        omega_0 = sat._compute_omega()
        s = station.position
        phi = sat.phase_angle

        # Pre-calculations,
        t = np.linspace(start=interval[0], stop=interval[1], num=delta_t) 
        R = rE/(rE + h)

        if sat.orbit_prograde:
            omega_rel = omega_0 - self.planetary_angular_velocity
        else:
            omega_rel = omega_0 + self.planetary_angular_velocity

        # Visibility calculations,
        delta = np.mod(omega_rel*t - s + phi, 2*np.pi)
        beta  = np.minimum(delta, 2*np.pi - delta)
        epsilon = np.arcsin((np.cos(beta)-R)/np.sqrt(1+R**2-2*R*np.cos(beta)))
        visible = epsilon >= 0

        # Detecting windows,
        changes = np.diff(visible.astype(int))
        starts = np.where(changes == 1)[0] + 1
        ends = np.where(changes == -1)[0]

        if visible[0]:
            starts = np.insert(starts, 0, 0)

        if visible[-1]:
            ends = np.append(ends, len(visible) - 1)

        windows = [(
                    sat.name,
                    t[s].item(), 
                    t[e].item(), 
                    t[e].item()-t[s].item(),
                    (t[e].item()-t[s].item())*sat.transfer_rate) 
                    for s, e in zip(starts, ends)]

        # Converting to HH:MM format,
        if format:
            windows = [(
                        sat_name,
                        self._sec_to_hms(start), 
                        self._sec_to_hms(end), 
                        self._sec_to_hms(duration),
                        round(capacity)
                        ) 
                        for sat_name, start, end, duration, capacity in windows]
        return windows

    def compute_windows(self, station, interval, delta_t=1000, format=False, order=True):
        """Compute the visibility of all satellites from a ground station at in the given time interval."""
        all_windows = []
        for sat in self.sats:
            windows = self.compute_window(sat, station, interval, delta_t=delta_t, format=format)
            all_windows.append(windows)
        if order:
            all_windows = [window for sublist in all_windows for window in sublist]
            all_windows.sort(key=lambda x: x[1])  # Sorting by start time.
        return all_windows

    def _sec_to_hms(self, t):
        """Helper function to convert seconds into the HH:MM:SS format."""
        t = int(round(t))
        h = t // 3600
        m = (t % 3600) // 60
        s = t % 60
        return f"{h:02d}:{m:02d}:{s:02d}"

class Satellite:
    """Class to represent a satellite with its properties such as name, altitude, and phase angle."""
    def __init__(self, name, altitude, phase_angle, transfer_rate, orbit_prograde=True):
        self.name = name
        self.altitude = altitude
        self.phase_angle = phase_angle
        self.transfer_rate = transfer_rate
        self.orbit_prograde = orbit_prograde
        self.planetary_radius = None
        self.planetary_mass = None
        self.planetary_angular_velocity = None
        self.G = None

    def get_position(self, t):
        """Calculate the satellite's position at time t. Returns theta which is the angle taken from the direction of 
         the positive x-axis in radians. Note that the reference frame is centered at the center of the planet."""
        omega_0 = self._compute_omega()
        if not self.orbit_prograde:
            theta = (omega_0 + self.planetary_angular_velocity) * t + self.phase_angle
        else:
            theta = (omega_0 - self.planetary_angular_velocity) * t + self.phase_angle
        return theta.item()

    def _compute_omega(self):
        """Compute the angular velocity of the satellite based on its altitude. The units of omega are in radians per second 
        (rad/s); the SI unit of angular velocity."""
        # Assuming circular orbit: omega_0 = sqrt(G*M/r^3)
        G = self.G
        M = self.planetary_mass
        r = (self.planetary_radius + self.altitude)*1e3  # Convert km to m
        return np.sqrt(G*M/ r**3)

    def _compute_orbital_period(self):
        """Compute the orbital period of the satellite based on its altitude. The units of T are in seconds (s); the SI unit of time."""
        omega = self._compute_omega()
        T = 2*np.pi / omega
        return T.item()

class Station:
    """Class to represent a ground station with its properties such as name and position."""
    def __init__(self, name, position):
        self.name = name
        self.position = position