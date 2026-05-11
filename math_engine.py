# pyrefly: ignore [missing-import]
import numpy as np

class BeamSolver:
    def __init__(self, spans, loads, end_conditions=('pin', 'pin')):
        """
        spans: list of span lengths [L1, L2, ..., Ln]
        loads: dict {span_idx: [{'type': 'uniform', 'w': w}, {'type': 'point', 'P': P, 'a': a}]}
               w: uniform load intensity (positive downwards)
               P: point load magnitude (positive downwards)
               a: distance of point load from the left support of the span
        end_conditions: tuple of strings for left and right ends (currently assumed 'pin' for simple supports)
        """
        self.spans = spans
        self.loads = loads
        self.n_spans = len(spans)
        self.n_supports = self.n_spans + 1
        self.end_conditions = end_conditions
        self.moments = np.zeros(self.n_supports)
        self.reactions = np.zeros(self.n_supports)
        
    def solve_moments(self):
        """Solves for internal bending moments at supports using the Three-Moment Equation."""
        n_eq = self.n_supports - 2
        if n_eq <= 0:
            return self.moments # Statically determinate or single span with pinned ends
            
        A = np.zeros((n_eq, n_eq))
        B = np.zeros(n_eq)
        
        for i in range(1, self.n_supports - 1):
            L_left = self.spans[i-1]
            L_right = self.spans[i]
            
            # Fill Matrix A
            if i > 1:
                A[i-1, i-2] = L_left
            A[i-1, i-1] = 2 * (L_left + L_right)
            if i < self.n_supports - 2:
                A[i-1, i] = L_right
                
            # Fill Vector B
            rhs = 0
            
            # Span to the left (span i-1)
            if (i-1) in self.loads:
                for load in self.loads[i-1]:
                    if load['type'] == 'uniform':
                        w = load['w']
                        rhs -= w * L_left**3 / 4
                    elif load['type'] == 'point':
                        P = load['P']
                        a = load['a'] # Distance from left support
                        rhs -= P * a / L_left * (L_left**2 - a**2)
                        
            # Span to the right (span i)
            if i in self.loads:
                for load in self.loads[i]:
                    if load['type'] == 'uniform':
                        w = load['w']
                        rhs -= w * L_right**3 / 4
                    elif load['type'] == 'point':
                        P = load['P']
                        a = load['a']
                        b = L_right - a # Distance from right support
                        rhs -= P * b / L_right * (L_right**2 - b**2)
                        
            B[i-1] = rhs
            
        # Solve Ax = B
        x = np.linalg.solve(A, B)
        self.moments[1:-1] = x
        return self.moments

    def compute_internal_forces(self, num_points=100):
        """Calculates arrays of position, shear force, and bending moment for visualization."""
        x_total = []
        V_total = []
        M_total = []
        
        current_x = 0
        
        # Calculate reactions as well
        self.reactions = np.zeros(self.n_supports)
        
        for i in range(self.n_spans):
            L = self.spans[i]
            x_span = np.linspace(0, L, num_points)
            
            # Base shear and moment for simply supported condition
            V_span = np.zeros_like(x_span)
            M_span = np.zeros_like(x_span)
            
            R_left_0 = 0
            R_right_0 = 0
            
            if i in self.loads:
                for load in self.loads[i]:
                    if load['type'] == 'uniform':
                        w = load['w']
                        R_left_0 += w * L / 2
                        R_right_0 += w * L / 2
                        
                        V_span += w * (L/2 - x_span)
                        M_span += w * x_span / 2 * (L - x_span)
                        
                    elif load['type'] == 'point':
                        P = load['P']
                        a = load['a']
                        b = L - a
                        
                        R_left_0 += P * b / L
                        R_right_0 += P * a / L
                        
                        V_span += np.where(x_span < a, P * b / L, -P * a / L)
                        M_span += np.where(x_span < a, P * b / L * x_span, P * a / L * (L - x_span))
            
            # Adjust for continuous beam end moments
            M_L = self.moments[i]
            M_R = self.moments[i+1]
            
            V_m = (M_L - M_R) / L
            V_span += V_m
            
            M_m = M_L + (M_R - M_L) * (x_span / L)
            M_span += M_m
            
            # Update reactions
            self.reactions[i] += R_left_0 + V_m
            self.reactions[i+1] += R_right_0 - V_m
            
            # Append to total lists
            x_total.extend(x_span + current_x)
            V_total.extend(V_span)
            M_total.extend(M_span)
            
            current_x += L
            
        return np.array(x_total), np.array(V_total), np.array(M_total)
