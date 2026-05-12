import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


class Span:
    def __init__(self, length):
        self.length = length
        self.loads = []

    def add_load(self, load_type, magnitude, distance=0):
        self.loads.append({'type': load_type, 'mag': magnitude, 'dist': distance})

    def rhs_left(self):
        total = 0
        for ld in self.loads:
            L = self.length
            if ld['type'] == 'uniform':
                total += ld['mag'] * L**3 / 4
            elif ld['type'] == 'point':
                a = ld['dist']
                total += ld['mag'] * a / L * (L**2 - a**2)
        return total

    def rhs_right(self):
        total = 0
        for ld in self.loads:
            L = self.length
            if ld['type'] == 'uniform':
                total += ld['mag'] * L**3 / 4
            elif ld['type'] == 'point':
                a, b = ld['dist'], L - ld['dist']
                total += ld['mag'] * b / L * (L**2 - b**2)
        return total


def solve_moments(spans, left_end='Pin', right_end='Pin'):
    n = len(spans) + 1
    moments = np.zeros(n)
    # Internal supports are always unknowns; end supports depend on type
    unknowns = list(range(1, n - 1))
    if left_end == 'Fixed':
        unknowns = [0] + unknowns
    if right_end == 'Fixed':
        unknowns = unknowns + [n - 1]
    n_eq = len(unknowns)
    if n_eq == 0:
        return moments
    idx_map = {s: j for j, s in enumerate(unknowns)}
    A = np.zeros((n_eq, n_eq))
    B = np.zeros(n_eq)
    for j, i in enumerate(unknowns):
        if i == 0:  # fixed left end — imaginary zero-length span to the left
            Ll, Lr = 0, spans[0].length
            rhs = -spans[0].rhs_right()
        elif i == n - 1:  # fixed right end — imaginary zero-length span to the right
            Ll, Lr = spans[-1].length, 0
            rhs = -spans[-1].rhs_left()
        else:  # internal support
            Ll, Lr = spans[i-1].length, spans[i].length
            rhs = -(spans[i-1].rhs_left() + spans[i].rhs_right())
        if i - 1 in idx_map:
            A[j, idx_map[i-1]] = Ll
        A[j, j] = 2 * (Ll + Lr)
        if i + 1 in idx_map:
            A[j, idx_map[i+1]] = Lr
        B[j] = rhs
    x = np.linalg.solve(A, B)
    for j, i in enumerate(unknowns):
        moments[i] = x[j]
    return moments


def compute_forces(spans, moments, pts=100):
    x_all, V_all, M_all = [], [], []
    reactions = np.zeros(len(spans) + 1)
    cx = 0
    for i, sp in enumerate(spans):
        L = sp.length
        x = np.linspace(0, L, pts)
        V, M = np.zeros(pts), np.zeros(pts)
        rl, rr = 0, 0
        for ld in sp.loads:
            if ld['type'] == 'uniform':
                w = ld['mag']
                rl += w * L / 2; rr += w * L / 2
                V += w * (L/2 - x); M += w * x / 2 * (L - x)
            elif ld['type'] == 'point':
                P, a, b = ld['mag'], ld['dist'], L - ld['dist']
                rl += P * b / L; rr += P * a / L
                V += np.where(x < a, P*b/L, -P*a/L)
                M += np.where(x < a, P*b/L*x, P*a/L*(L-x))
        vm = (moments[i] - moments[i+1]) / L
        V += vm; M += moments[i] + (moments[i+1] - moments[i]) * x / L
        reactions[i] += rl + vm; reactions[i+1] += rr - vm
        x_all.extend(x + cx); V_all.extend(V); M_all.extend(M)
        cx += L
    return np.array(x_all), np.array(V_all), np.array(M_all), reactions


class BeamApp:
    def __init__(self, root):
        self.root = root
        self.root.title("CE 257 Beam Solver (Three-Moment Equation)")
        self.root.geometry("1050x720")
        self.spans = []
        self.setup_ui()
        self.root.after(100, self.draw_beam)

    def setup_ui(self):
        left = tk.Frame(self.root, width=340)
        left.pack(side=tk.LEFT, fill=tk.Y, padx=8, pady=8)
        left.pack_propagate(False)

        # Span entry
        sf = ttk.LabelFrame(left, text="Beam Geometry", padding=6)
        sf.pack(fill=tk.X, pady=(0, 6))
        ttk.Label(sf, text="Span Length (m):").grid(row=0, column=0, padx=4, pady=3)
        self.ent_len = tk.Entry(sf, width=10)
        self.ent_len.grid(row=0, column=1, padx=4, pady=3)
        ttk.Button(sf, text="Add Span", command=self.add_span).grid(row=0, column=2, padx=4)
        ttk.Label(sf, text="Defined Spans (click to select):").grid(row=1, column=0, columnspan=3, sticky="w", padx=4, pady=(6,2))
        self.span_list = tk.Listbox(sf, height=5, exportselection=False)
        self.span_list.grid(row=2, column=0, columnspan=3, sticky="we", padx=4, pady=2)

        # Support conditions
        sup = ttk.LabelFrame(left, text="Support Conditions", padding=6)
        sup.pack(fill=tk.X, pady=(0, 6))
        types = ["Pin", "Roller", "Fixed", "Free End"]
        ttk.Label(sup, text="Left End:").grid(row=0, column=0, padx=4, pady=3, sticky="w")
        self.left_sup = tk.StringVar(value="Pin")
        ttk.Combobox(sup, textvariable=self.left_sup, values=types, width=10,
                     state="readonly").grid(row=0, column=1, padx=4, pady=3)
        ttk.Label(sup, text="Right End:").grid(row=1, column=0, padx=4, pady=3, sticky="w")
        self.right_sup = tk.StringVar(value="Pin")
        ttk.Combobox(sup, textvariable=self.right_sup, values=types, width=10,
                     state="readonly").grid(row=1, column=1, padx=4, pady=3)

        # Load entry
        lf = ttk.LabelFrame(left, text="Applied Loading", padding=6)
        lf.pack(fill=tk.X, pady=(0, 6))
        ttk.Label(lf, text="Load Type:").grid(row=0, column=0, padx=4, pady=3, sticky="w")
        self.load_type = tk.StringVar(value="uniform")
        ttk.Combobox(lf, textvariable=self.load_type, values=["uniform", "point"],
                     width=10, state="readonly").grid(row=0, column=1, padx=4, pady=3)
        ttk.Label(lf, text="Magnitude (kN or kN/m):").grid(row=1, column=0, padx=4, pady=3, sticky="w")
        self.ent_mag = tk.Entry(lf, width=10)
        self.ent_mag.grid(row=1, column=1, padx=4, pady=3)
        ttk.Label(lf, text="Distance 'a' (m) [point only]:").grid(row=2, column=0, padx=4, pady=3, sticky="w")
        self.ent_dist = tk.Entry(lf, width=10)
        self.ent_dist.grid(row=2, column=1, padx=4, pady=3)
        ttk.Button(lf, text="Add Load to Selected Span", command=self.add_load).grid(
            row=3, column=0, columnspan=2, pady=6, sticky="we", padx=4)

        # Action buttons
        af = ttk.Frame(left, padding=4)
        af.pack(fill=tk.X, pady=(0, 6))
        ttk.Button(af, text="Clear All", command=self.clear_all).pack(side=tk.LEFT, padx=4)
        tk.Button(af, text="Solve & Plot BMD/SFD", bg="#27ae60", fg="white",
                  font=("Arial", 10, "bold"), command=self.solve).pack(side=tk.RIGHT, padx=4)

        # Status log
        logf = ttk.LabelFrame(left, text="Status Log", padding=4)
        logf.pack(fill=tk.BOTH, expand=True)
        self.log_text = tk.Text(logf, height=6, width=38, wrap=tk.WORD)
        self.log_text.pack(fill=tk.BOTH, expand=True)
        self.log_text.insert(tk.END, "Status: Ready\n")
        self.log_text.config(state=tk.DISABLED)

        # Right panel
        right = tk.Frame(self.root)
        right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=8, pady=8)
        cf = ttk.LabelFrame(right, text="Beam Visualization", padding=4)
        cf.pack(fill=tk.X, pady=(0, 4))
        self.canvas = tk.Canvas(cf, bg="white", height=170)
        self.canvas.pack(fill=tk.X)
        self.plot_frame = ttk.LabelFrame(right, text="Analysis Results", padding=4)
        self.plot_frame.pack(fill=tk.BOTH, expand=True)

    def log(self, msg):
        self.log_text.config(state=tk.NORMAL)
        self.log_text.insert(tk.END, msg + "\n")
        self.log_text.see(tk.END)
        self.log_text.config(state=tk.DISABLED)

    def add_span(self):
        try:
            L = float(self.ent_len.get())
            if L <= 0: raise ValueError
            self.spans.append(Span(L))
            self.span_list.insert(tk.END, f"Span {len(self.spans)}: L = {L} m")
            self.log(f"Added span {L} m. Total spans: {len(self.spans)}")
            self.draw_beam(); self.ent_len.delete(0, tk.END)
        except ValueError:
            messagebox.showerror("Error", "Enter a valid positive span length.")

    def add_load(self):
        sel = self.span_list.curselection()
        if not sel:
            messagebox.showwarning("Warning", "Select a span from the list first.")
            return
        idx = sel[0]; sp = self.spans[idx]
        try:
            mag = float(self.ent_mag.get())
            lt = self.load_type.get()
            if lt == 'point':
                d = float(self.ent_dist.get())
                if d < 0 or d > sp.length:
                    messagebox.showerror("Error", f"Distance must be 0–{sp.length} m."); return
                sp.add_load('point', mag, d)
                self.log(f"Point load {mag} kN at {d} m on Span {idx+1}.")
            else:
                sp.add_load('uniform', mag)
                self.log(f"UDL {mag} kN/m on Span {idx+1}.")
            self.draw_beam(); self.ent_mag.delete(0, tk.END); self.ent_dist.delete(0, tk.END)
        except ValueError:
            messagebox.showerror("Error", "Enter valid numbers for load parameters.")

    def clear_all(self):
        self.spans = []; self.span_list.delete(0, tk.END)
        self.left_sup.set("Pin"); self.right_sup.set("Pin")
        self.log("Cleared all."); self.draw_beam()
        for w in self.plot_frame.winfo_children(): w.destroy()

    # ── Support Drawing ──────────────────────────────────────────────────
    def _draw_support(self, x, y, stype):
        if stype == "Pin":
            self.canvas.create_polygon(x, y, x-10, y+18, x+10, y+18, fill="#2ecc71", outline="black")
            self.canvas.create_line(x-14, y+18, x+14, y+18, width=3)
        elif stype == "Roller":
            r = 7
            self.canvas.create_oval(x-r, y+4, x+r, y+4+2*r, fill="#f1c40f", outline="black", width=2)
            self.canvas.create_line(x-14, y+4+2*r+2, x+14, y+4+2*r+2, width=3)
        elif stype == "Fixed":
            self.canvas.create_rectangle(x-4, y-20, x+4, y+20, fill="#7f8c8d", outline="black", width=2)
            for dy in range(-18, 20, 8):
                self.canvas.create_line(x-4, y+dy, x-14, y+dy+8, width=1.5, fill="#555")
        elif stype == "Free End":
            self.canvas.create_oval(x-4, y-4, x+4, y+4, fill="#e74c3c", outline="black", width=2)

    def draw_beam(self):
        self.canvas.delete("all"); self.canvas.update()
        cw, ch = self.canvas.winfo_width(), self.canvas.winfo_height()
        if cw <= 1: cw, ch = 680, 170
        if not self.spans:
            self.canvas.create_text(cw/2, ch/2, text="Add spans to begin.", fill="gray", font=("Arial", 12))
            return
        total = sum(s.length for s in self.spans)
        scale, sx, y = (cw - 100) / total, 50, ch / 2
        self.canvas.create_line(sx, y, sx + total*scale, y, width=5, fill="#3498db")
        n_sup = len(self.spans) + 1
        cx = sx
        # Draw left-end support
        self._draw_support(cx, y, self.left_sup.get())
        for i, sp in enumerate(self.spans):
            ex = cx + sp.length * scale
            # Intermediate supports are Pin; last support uses right-end type
            sup_type = self.right_sup.get() if i == len(self.spans) - 1 else "Pin"
            self._draw_support(ex, y, sup_type)
            self.canvas.create_text((cx+ex)/2, y+30, text=f"{sp.length} m", fill="#2c3e50", font=("Arial", 9, "bold"))
            for ld in sp.loads:
                if ld['type'] == 'uniform':
                    pw = sp.length * scale
                    for j in range(max(5, int(pw/25)) + 1):
                        lx = cx + j / max(5, int(pw/25)) * pw
                        self.canvas.create_line(lx, y-35, lx, y, arrow=tk.LAST, fill="#e74c3c", width=2)
                    self.canvas.create_text(cx + pw/2, y-45, text=f"{ld['mag']} kN/m", fill="#c0392b", font=("Arial", 9, "bold"))
                elif ld['type'] == 'point':
                    lx = cx + ld['dist'] * scale
                    self.canvas.create_line(lx, y-55, lx, y, arrow=tk.LAST, fill="#e74c3c", width=3)
                    self.canvas.create_text(lx, y-63, text=f"{ld['mag']} kN", fill="#c0392b", font=("Arial", 9, "bold"))
            cx = ex

    def solve(self):
        if not self.spans:
            messagebox.showwarning("Warning", "Add at least one span first."); return
        try:
            self.log("Solving...")
            moms = solve_moments(self.spans, self.left_sup.get(), self.right_sup.get())
            for i, m in enumerate(moms): self.log(f"  M{i} = {m:.2f} kNm")
            x, V, M, rcts = compute_forces(self.spans, moms)
            for i, r in enumerate(rcts): self.log(f"  R{i} = {r:.2f} kN")
            self.plot(x, V, M); self.log("Done.")
        except Exception as e:
            messagebox.showerror("Error", f"Solution failed:\n{e}")

    def plot(self, x, V, M):
        for w in self.plot_frame.winfo_children(): w.destroy()
        fig, (a1, a2) = plt.subplots(2, 1, figsize=(6, 4.5), sharex=True)
        a1.plot(x, V, color='#2980b9', lw=2); a1.fill_between(x, V, 0, color='#2980b9', alpha=0.15)
        a1.set_title('Shear Force Diagram (SFD)'); a1.set_ylabel('V (kN)'); a1.grid(True, ls='--', alpha=0.5); a1.axhline(0, color='k', lw=0.8)
        a2.plot(x, M, color='#c0392b', lw=2); a2.fill_between(x, M, 0, color='#c0392b', alpha=0.15)
        a2.set_title('Bending Moment Diagram (BMD)'); a2.set_xlabel('Position (m)'); a2.set_ylabel('M (kNm)')
        a2.grid(True, ls='--', alpha=0.5); a2.axhline(0, color='k', lw=0.8); a2.invert_yaxis()
        sx = 0
        for sp in self.spans:
            a1.axvline(sx, color='gray', ls=':', lw=1); a2.axvline(sx, color='gray', ls=':', lw=1); sx += sp.length
        a1.axvline(sx, color='gray', ls=':', lw=1); a2.axvline(sx, color='gray', ls=':', lw=1)
        fig.tight_layout()
        c = FigureCanvasTkAgg(fig, master=self.plot_frame); c.draw(); c.get_tk_widget().pack(fill=tk.BOTH, expand=True)


root = tk.Tk()
app = BeamApp(root)
root.mainloop()
