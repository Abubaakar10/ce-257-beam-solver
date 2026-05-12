# CE 257 Beam Solver — Three-Moment Equation Analysis Tool

A Python desktop application for analyzing **continuous beams** using the **Three-Moment Equation (Clapeyron's Theorem)**. Built with Tkinter and Matplotlib, designed for civil/structural engineering coursework.

![Python](https://img.shields.io/badge/Python-3.8%2B-blue)
![License](https://img.shields.io/badge/License-MIT-green)

---

## Features

- **Three-Moment Equation Solver** — Automatically constructs and solves the system of equations for statically indeterminate continuous beams
- **Point Loads & UDL Support** — Apply point loads at any position or uniformly distributed loads on any span
- **Interactive Span Selection** — Select spans from a listbox to apply loads (no manual index typing)
- **Live Beam Visualization** — See your beam, supports, and loads drawn in real-time as you build
- **Embedded SFD & BMD Plots** — Shear Force and Bending Moment diagrams displayed directly inside the app window
- **Support Reactions & Moments** — All computed values logged in the status panel
- **Single-File Design** — The entire application is contained in one Python file (`beam_solver.py`)

---

## Screenshots

<!-- Add screenshots of your app here after running it -->
<!-- ![Beam Solver Screenshot](screenshots/app.png) -->

---

## Getting Started

### Prerequisites

- **Python 3.8+**
- **Spyder IDE** (recommended) or any Python environment

### Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/Abubaakar10/ce-257-beam-solver.git
   cd ce-257-beam-solver
   ```

2. **Create a virtual environment (optional but recommended):**
   ```bash
   python -m venv .venv
   .venv\Scripts\activate    # Windows
   source .venv/bin/activate  # macOS/Linux
   ```

3. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

### Running the App

**Option A — Spyder (Recommended):**
1. Open `beam_solver.py` in Spyder
2. Press **F5** (or click Run)

**Option B — Command Line:**
```bash
python beam_solver.py
```

---

## How to Use

1. **Add Spans** — Enter a span length in meters and click "Add Span". Repeat for each span of your continuous beam.
2. **Select a Span** — Click on a span in the listbox to select it.
3. **Add Loads** — Choose load type (uniform or point), enter the magnitude, and distance (for point loads only). Click "Add Load to Selected Span".
4. **Solve** — Click the green **"Solve & Plot BMD/SFD"** button to run the analysis.
5. **View Results** — The SFD and BMD appear in the right panel. Support moments and reactions are printed in the status log.
6. **Clear All** — Reset everything and start a new beam.

---

## Project Structure

```
ce-257-beam-solver/
├── beam_solver.py      # All-in-one application (solver + GUI + plotting)
├── report.md           # Technical implementation report
├── requirements.txt    # Python dependencies (numpy, matplotlib)
├── .gitignore
└── README.md
```

---

## Theory: Three-Moment Equation

For a continuous beam with spans $L_i$ and $L_{i+1}$, the Three-Moment Equation relates the bending moments at three consecutive supports:

$$M_{i-1}L_{i} + 2M_{i}(L_{i} + L_{i+1}) + M_{i+1}L_{i+1} = -6 \left( \frac{A_i \bar{x}_i}{L_i} + \frac{A_{i+1} \bar{x}_{i+1}}{L_{i+1}} \right)$$

The application assembles this into a matrix equation **Ax = B** and solves it using `numpy.linalg.solve`.

For more details, see [`report.md`](report.md).

---

## Dependencies

| Package | Purpose |
|---------|---------|
| `numpy` | Matrix algebra & linear system solver |
| `matplotlib` | SFD/BMD plotting, embedded in Tkinter |
| `tkinter` | GUI framework (included with Python) |

---

## Author

**Abubakar** — CE 257 Capstone Project

---

## License

This project is open source and available under the [MIT License](LICENSE).
