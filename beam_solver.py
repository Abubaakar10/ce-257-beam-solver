# =============================================================================
# CE 257 Beam Solver — Three-Moment Equation Analysis Tool
# All-in-one single-file application for continuous beam analysis
# =============================================================================

import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


# ── Data Classes ─────────────────────────────────────────────────────────────

class Load:
    """Represents a single load applied on a span."""
    def __init__(self, load_type, magnitude, distance=0):
        self.type = load_type        # 'uniform' or 'point'
        self.magnitude = magnitude   # kN/m for UDL, kN for point
        self.distance = distance     # distance 'a' from left support (point loads only)

    def rhs_left(self, L):
        """RHS contribution when this load's span is to the LEFT of the support."""
        if self.type == 'uniform':
            return self.magnitude * L**3 / 4
        elif self.type == 'point':
            a = self.distance
            return self.magnitude * a / L * (L**2 - a**2)
        return 0

    def rhs_right(self, L):
        """RHS contribution when this load's span is to the RIGHT of the support."""
        if self.type == 'uniform':
            return self.magnitude * L**3 / 4
        elif self.type == 'point':
            a = self.distance
            b = L - a
            return self.magnitude * b / L * (L**2 - b**2)
        return 0


class Span:
    """Represents a single span of the continuous beam."""
    def __init__(self, length):
        self.length = length
        self.loads = []

    def add_load(self, load):
        self.loads.append(load)

    def total_rhs_left(self):
        """Sum of RHS terms when this span is on the left side of a support."""
        return sum(ld.rhs_left(self.length) for ld in self.loads)

    def total_rhs_right(self):
        """Sum of RHS terms when this span is on the right side of a support."""
        return sum(ld.rhs_right(self.length) for ld in self.loads)


# ── Solver ───────────────────────────────────────────────────────────────────

def solve_three_moment(spans):
    """
    Solves support moments for a continuous beam using the Three-Moment Equation.
    Assumes simple supports (M=0) at both ends.
    Returns a list of moments at each support [M0, M1, ..., Mn].
    """
    n_spans = len(spans)
    n_supports = n_spans + 1
    moments = np.zeros(n_supports)

    n_eq = n_supports - 2  # number of internal supports
    if n_eq <= 0:
        return moments  # single span — simply supported

    A = np.zeros((n_eq, n_eq))
    B = np.zeros(n_eq)

    for i in range(1, n_supports - 1):
        L_left = spans[i - 1].length
        L_right = spans[i].length

        # Coefficient matrix
        if i > 1:
            A[i - 1, i - 2] = L_left
        A[i - 1, i - 1] = 2 * (L_left + L_right)
        if i < n_supports - 2:
            A[i - 1, i] = L_right

        # Right-hand side
        rhs = 0
        rhs -= spans[i - 1].total_rhs_left()
        rhs -= spans[i].total_rhs_right()
        B[i - 1] = rhs

    x = np.linalg.solve(A, B)
    moments[1:-1] = x
    return moments


def compute_internal_forces(spans, moments, num_points=100):
    """
    Computes arrays of position, shear force, and bending moment along the beam.
    Returns (x_array, V_array, M_array, reactions).
    """
    x_total = []
    V_total = []
    M_total = []
    n_supports = len(spans) + 1
    reactions = np.zeros(n_supports)
    current_x = 0

    for i, span in enumerate(spans):
        L = span.length
        x_span = np.linspace(0, L, num_points)

        V_span = np.zeros_like(x_span)
        M_span = np.zeros_like(x_span)

        R_left_0 = 0
        R_right_0 = 0

        for ld in span.loads:
            if ld.type == 'uniform':
                w = ld.magnitude
                R_left_0 += w * L / 2
                R_right_0 += w * L / 2
                V_span += w * (L / 2 - x_span)
                M_span += w * x_span / 2 * (L - x_span)

            elif ld.type == 'point':
                P = ld.magnitude
                a = ld.distance
                b = L - a
                R_left_0 += P * b / L
                R_right_0 += P * a / L
                V_span += np.where(x_span < a, P * b / L, -P * a / L)
                M_span += np.where(x_span < a,
                                   P * b / L * x_span,
                                   P * a / L * (L - x_span))

        # Adjustments from continuous beam end moments
        M_L = moments[i]
        M_R = moments[i + 1]

        V_m = (M_L - M_R) / L
        V_span += V_m

        M_m = M_L + (M_R - M_L) * (x_span / L)
        M_span += M_m

        reactions[i] += R_left_0 + V_m
        reactions[i + 1] += R_right_0 - V_m

        x_total.extend(x_span + current_x)
        V_total.extend(V_span)
        M_total.extend(M_span)

        current_x += L

    return np.array(x_total), np.array(V_total), np.array(M_total), reactions


# ── GUI Application ──────────────────────────────────────────────────────────

class BeamApp:
    def __init__(self, root):
        self.root = root
        self.root.title("CE 257 Beam Solver (Three-Moment Equation)")
        self.root.geometry("1100x750")

        self.spans = []           # list of Span objects
        self.support_moments = [] # solved moments at each support

        self.setup_ui()
        # Give the canvas time to render before first draw
        self.root.after(100, self.draw_beam_canvas)

    # ── UI Layout ────────────────────────────────────────────────────────

    def setup_ui(self):
        # ── Left Panel: Controls ─────────────────────────────────────────
        left_panel = tk.Frame(self.root, width=360)
        left_panel.pack(side=tk.LEFT, fill=tk.Y, padx=8, pady=8)
        left_panel.pack_propagate(False)

        # --- Section 1: Add Span ---
        span_frame = ttk.LabelFrame(left_panel, text="Beam Geometry", padding=8)
        span_frame.pack(fill=tk.X, pady=(0, 8))

        ttk.Label(span_frame, text="Span Length (m):").grid(row=0, column=0,
                                                            padx=5, pady=4, sticky="w")
        self.span_length_var = tk.StringVar()
        ttk.Entry(span_frame, textvariable=self.span_length_var,
                  width=12).grid(row=0, column=1, padx=5, pady=4)
        ttk.Button(span_frame, text="Add Span",
                   command=self.add_span).grid(row=0, column=2, padx=5, pady=4)

        # Span Listbox (used for both display AND selection when adding loads)
        ttk.Label(span_frame, text="Defined Spans:").grid(row=1, column=0,
                                                          columnspan=3, sticky="w",
                                                          padx=5, pady=(8, 2))
        self.span_listbox = tk.Listbox(span_frame, height=5,
                                       exportselection=False,
                                       selectmode=tk.SINGLE)
        self.span_listbox.grid(row=2, column=0, columnspan=3,
                               sticky="we", padx=5, pady=2)

        # --- Section 2: Add Load ---
        load_frame = ttk.LabelFrame(left_panel, text="Applied Loading", padding=8)
        load_frame.pack(fill=tk.X, pady=(0, 8))

        ttk.Label(load_frame, text="Load Type:").grid(row=0, column=0,
                                                      padx=5, pady=4, sticky="w")
        self.load_type_var = tk.StringVar(value="uniform")
        type_combo = ttk.Combobox(load_frame, textvariable=self.load_type_var,
                                  values=["uniform", "point"], width=10,
                                  state="readonly")
        type_combo.grid(row=0, column=1, padx=5, pady=4)
        type_combo.bind("<<ComboboxSelected>>", self.on_load_type_change)

        ttk.Label(load_frame, text="Magnitude (kN or kN/m):").grid(row=1, column=0,
                                                                    padx=5, pady=4,
                                                                    sticky="w")
        self.load_mag_var = tk.StringVar()
        ttk.Entry(load_frame, textvariable=self.load_mag_var,
                  width=12).grid(row=1, column=1, padx=5, pady=4)

        self.dist_label = ttk.Label(load_frame, text="Distance 'a' (m):")
        self.dist_label.grid(row=2, column=0, padx=5, pady=4, sticky="w")
        self.load_dist_var = tk.StringVar()
        self.dist_entry = ttk.Entry(load_frame, textvariable=self.load_dist_var,
                                    width=12)
        self.dist_entry.grid(row=2, column=1, padx=5, pady=4)

        # Hide distance field initially (uniform load selected by default)
        self.dist_label.grid_remove()
        self.dist_entry.grid_remove()

        ttk.Label(load_frame,
                  text="↑ Select a span from the list above,\n   then fill in load details.",
                  foreground="gray").grid(row=3, column=0, columnspan=2,
                                          padx=5, pady=4, sticky="w")

        ttk.Button(load_frame, text="Add Load to Selected Span",
                   command=self.add_load).grid(row=4, column=0, columnspan=2,
                                               pady=8, sticky="we", padx=5)

        # --- Section 3: Action Buttons ---
        action_frame = ttk.Frame(left_panel, padding=4)
        action_frame.pack(fill=tk.X, pady=(0, 8))

        ttk.Button(action_frame, text="Clear All",
                   command=self.clear_all).pack(side=tk.LEFT, padx=5)
        solve_btn = tk.Button(action_frame, text="Solve & Plot BMD/SFD",
                              bg="#27ae60", fg="white",
                              font=("Arial", 10, "bold"),
                              command=self.solve_and_plot)
        solve_btn.pack(side=tk.RIGHT, padx=5)

        # --- Section 4: Status Log ---
        log_frame = ttk.LabelFrame(left_panel, text="Status Log", padding=4)
        log_frame.pack(fill=tk.BOTH, expand=True)

        self.info_text = tk.Text(log_frame, height=8, width=40, wrap=tk.WORD)
        self.info_text.pack(fill=tk.BOTH, expand=True)
        self.info_text.insert(tk.END, "Status: Ready\n")
        self.info_text.config(state=tk.DISABLED)

        # ── Right Panel: Beam Canvas + Embedded Plots ────────────────────
        right_panel = tk.Frame(self.root)
        right_panel.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=8, pady=8)

        # Beam visualization canvas
        canvas_frame = ttk.LabelFrame(right_panel, text="Beam Visualization",
                                      padding=4)
        canvas_frame.pack(fill=tk.X, pady=(0, 4))

        self.canvas = tk.Canvas(canvas_frame, bg="white", height=180)
        self.canvas.pack(fill=tk.X)

        # Plot area (SFD + BMD embedded in the same window)
        self.plot_frame = ttk.LabelFrame(right_panel, text="Analysis Results",
                                         padding=4)
        self.plot_frame.pack(fill=tk.BOTH, expand=True)

        # Placeholder label when no results yet
        self.placeholder_label = ttk.Label(
            self.plot_frame,
            text="Add spans and loads, then click 'Solve & Plot BMD/SFD'\n"
                 "to display the analysis results here.",
            foreground="gray",
            anchor="center",
            justify="center"
        )
        self.placeholder_label.pack(expand=True)

    # ── Event Handlers ───────────────────────────────────────────────────

    def on_load_type_change(self, event=None):
        """Show/hide the distance field depending on load type."""
        if self.load_type_var.get() == "point":
            self.dist_label.grid()
            self.dist_entry.grid()
        else:
            self.dist_label.grid_remove()
            self.dist_entry.grid_remove()

    def log(self, message):
        """Append a message to the status log."""
        self.info_text.config(state=tk.NORMAL)
        self.info_text.insert(tk.END, message + "\n")
        self.info_text.see(tk.END)
        self.info_text.config(state=tk.DISABLED)

    # ── Add / Clear Actions ──────────────────────────────────────────────

    def add_span(self):
        try:
            length = float(self.span_length_var.get())
            if length <= 0:
                raise ValueError
            new_span = Span(length)
            self.spans.append(new_span)
            self.span_listbox.insert(tk.END,
                                     f"Span {len(self.spans)}: L = {length} m")
            self.log(f"Added span of length {length} m. "
                     f"Total spans: {len(self.spans)}")
            self.draw_beam_canvas()
            self.span_length_var.set("")
        except ValueError:
            messagebox.showerror("Error",
                                 "Please enter a valid positive number for span length.")

    def add_load(self):
        # Check that a span is selected in the listbox
        selection = self.span_listbox.curselection()
        if not selection:
            messagebox.showwarning("Warning",
                                   "Please select a span from the list before adding a load.")
            return

        span_idx = selection[0]
        span = self.spans[span_idx]

        try:
            mag = float(self.load_mag_var.get())
            l_type = self.load_type_var.get()

            if l_type == 'point':
                dist_str = self.load_dist_var.get()
                if not dist_str:
                    messagebox.showerror("Error",
                                         "Please enter distance 'a' for the point load.")
                    return
                dist = float(dist_str)
                if dist < 0 or dist > span.length:
                    messagebox.showerror(
                        "Error",
                        f"Distance 'a' must be between 0 and {span.length} m "
                        f"(the length of the selected span).")
                    return
                new_load = Load('point', mag, dist)
                span.add_load(new_load)
                self.log(f"Added point load {mag} kN at {dist} m on Span {span_idx + 1}.")
            else:
                new_load = Load('uniform', mag)
                span.add_load(new_load)
                self.log(f"Added uniform load {mag} kN/m on Span {span_idx + 1}.")

            self.draw_beam_canvas()
            self.load_mag_var.set("")
            self.load_dist_var.set("")

        except ValueError:
            messagebox.showerror("Error",
                                 "Please enter valid numbers for load parameters.")

    def clear_all(self):
        self.spans = []
        self.span_listbox.delete(0, tk.END)
        self.log("Cleared all spans and loads.")
        self.draw_beam_canvas()
        # Reset plot area
        for widget in self.plot_frame.winfo_children():
            widget.destroy()
        self.placeholder_label = ttk.Label(
            self.plot_frame,
            text="Add spans and loads, then click 'Solve & Plot BMD/SFD'\n"
                 "to display the analysis results here.",
            foreground="gray", anchor="center", justify="center"
        )
        self.placeholder_label.pack(expand=True)

    # ── Beam Canvas Drawing ──────────────────────────────────────────────

    def draw_beam_canvas(self):
        """Draws the beam, supports, and loads on the canvas."""
        self.canvas.delete("all")
        self.canvas.update()
        w = self.canvas.winfo_width()
        h = self.canvas.winfo_height()

        if w <= 1:
            w, h = 700, 180

        # Background grid
        for gx in range(0, w, 20):
            self.canvas.create_line(gx, 0, gx, h, fill="#f0f0f0")
        for gy in range(0, h, 20):
            self.canvas.create_line(0, gy, w, gy, fill="#f0f0f0")

        if not self.spans:
            self.canvas.create_text(w / 2, h / 2,
                                    text="Enter a span length to begin.",
                                    fill="gray", font=("Arial", 13))
            return

        total_length = sum(s.length for s in self.spans)
        if total_length == 0:
            return

        margin = 60
        scale = (w - 2 * margin) / total_length
        start_x = margin
        y = h / 2

        # Draw beam line
        end_x = start_x + total_length * scale
        self.canvas.create_line(start_x, y, end_x, y, width=6, fill="#3498db")

        # Draw supports and span labels
        current_x = start_x
        self.draw_support(current_x, y)

        for i, span in enumerate(self.spans):
            span_end_x = current_x + span.length * scale
            self.draw_support(span_end_x, y)

            # Span length label
            mid_x = (current_x + span_end_x) / 2
            self.canvas.create_text(mid_x, y + 35,
                                    text=f"{span.length} m",
                                    fill="#2c3e50",
                                    font=("Arial", 9, "bold"))

            # Draw loads for this span
            for ld in span.loads:
                if ld.type == 'uniform':
                    span_len_px = span.length * scale
                    num_arrows = max(5, int(span_len_px / 25))
                    for j in range(num_arrows + 1):
                        lx = current_x + (j / num_arrows) * span_len_px
                        self.canvas.create_line(lx, y - 40, lx, y,
                                                arrow=tk.LAST, fill="#e74c3c",
                                                width=2)
                    # Label
                    self.canvas.create_text(
                        current_x + span_len_px / 2, y - 50,
                        text=f"{ld.magnitude} kN/m",
                        fill="#c0392b", font=("Arial", 9, "bold"))

                elif ld.type == 'point':
                    lx = current_x + ld.distance * scale
                    self.canvas.create_line(lx, y - 60, lx, y,
                                            arrow=tk.LAST, fill="#e74c3c",
                                            width=3)
                    self.canvas.create_text(lx, y - 68,
                                            text=f"{ld.magnitude} kN",
                                            fill="#c0392b",
                                            font=("Arial", 9, "bold"))

            current_x = span_end_x

    def draw_support(self, x, y):
        """Draws a triangular support at position (x, y)."""
        self.canvas.create_polygon(x, y, x - 10, y + 20, x + 10, y + 20,
                                   fill="#2ecc71", outline="black")
        self.canvas.create_line(x - 15, y + 20, x + 15, y + 20,
                                width=3, fill="black")

    # ── Solve & Plot ─────────────────────────────────────────────────────

    def solve_and_plot(self):
        if not self.spans:
            messagebox.showwarning("Warning",
                                   "Please add at least one span before solving.")
            return

        try:
            self.log("Solving using Three-Moment Equation...")

            # Solve moments
            self.support_moments = solve_three_moment(self.spans)

            self.log("Moments at supports:")
            for i, m in enumerate(self.support_moments):
                self.log(f"  Support {i} (M{i}): {m:.2f} kNm")

            # Compute internal forces
            x, V, M, reactions = compute_internal_forces(
                self.spans, self.support_moments)

            self.log("Reactions at supports:")
            for i, r in enumerate(reactions):
                self.log(f"  Support {i} (R{i}): {r:.2f} kN")

            self.log("Plotting BMD and SFD...")
            self.plot_results(x, V, M)
            self.log("Done — results displayed.")

        except Exception as e:
            messagebox.showerror("Error",
                                 f"An error occurred during solution:\n{str(e)}")

    def plot_results(self, x, V, M):
        """Embeds the SFD and BMD matplotlib plots inside the Tkinter window."""
        # Clear previous plots
        for widget in self.plot_frame.winfo_children():
            widget.destroy()

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 5), sharex=True)
        fig.patch.set_facecolor("#fafafa")

        # ── Shear Force Diagram ──
        ax1.plot(x, V, color='#2980b9', linewidth=2)
        ax1.fill_between(x, V, 0, color='#2980b9', alpha=0.15)
        ax1.set_title('Shear Force Diagram (SFD)', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Shear Force (kN)', fontsize=10)
        ax1.grid(True, linestyle='--', alpha=0.5)
        ax1.axhline(0, color='black', linewidth=0.8)

        # ── Bending Moment Diagram ──
        ax2.plot(x, M, color='#c0392b', linewidth=2)
        ax2.fill_between(x, M, 0, color='#c0392b', alpha=0.15)
        ax2.set_title('Bending Moment Diagram (BMD)', fontsize=12, fontweight='bold')
        ax2.set_xlabel('Position along beam (m)', fontsize=10)
        ax2.set_ylabel('Bending Moment (kNm)', fontsize=10)
        ax2.grid(True, linestyle='--', alpha=0.5)
        ax2.axhline(0, color='black', linewidth=0.8)
        ax2.invert_yaxis()

        # Mark support locations
        support_x = 0
        for sp in self.spans:
            ax1.axvline(support_x, color='gray', linestyle=':', linewidth=1.2)
            ax2.axvline(support_x, color='gray', linestyle=':', linewidth=1.2)
            support_x += sp.length
        ax1.axvline(support_x, color='gray', linestyle=':', linewidth=1.2)
        ax2.axvline(support_x, color='gray', linestyle=':', linewidth=1.2)

        fig.tight_layout()

        # Embed into Tkinter
        canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)


# ── Main Execution ───────────────────────────────────────────────────────────

root = tk.Tk()
app = BeamApp(root)
root.mainloop()
