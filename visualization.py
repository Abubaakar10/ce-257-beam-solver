# pyrefly: ignore [missing-import]
import matplotlib.pyplot as plt

def plot_diagrams(x, V, M, spans):
    """
    Plots the Shear Force Diagram (SFD) and Bending Moment Diagram (BMD).
    
    x: array of positions along the beam
    V: array of shear force values
    M: array of bending moment values
    spans: list of span lengths to mark supports
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Plot Shear Force Diagram (SFD)
    ax1.plot(x, V, color='blue', linewidth=2)
    ax1.fill_between(x, V, 0, color='blue', alpha=0.2)
    ax1.set_title('Shear Force Diagram (SFD)', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Shear Force (kN)', fontsize=12)
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.axhline(0, color='black', linewidth=1)
    
    # Plot Bending Moment Diagram (BMD)
    # Note: By standard civil engineering conventions, positive bending moment (sagging) 
    # is often plotted on the tension side (below the axis). We invert the y-axis for this effect.
    ax2.plot(x, M, color='red', linewidth=2)
    ax2.fill_between(x, M, 0, color='red', alpha=0.2)
    ax2.set_title('Bending Moment Diagram (BMD)', fontsize=14, fontweight='bold')
    ax2.set_xlabel('Position along beam (m)', fontsize=12)
    ax2.set_ylabel('Bending Moment (kNm)', fontsize=12)
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.axhline(0, color='black', linewidth=1)
    
    # Invert y-axis for BMD
    ax2.invert_yaxis()
    
    # Mark support locations with vertical dashed lines
    support_x = 0
    ax1.axvline(support_x, color='gray', linestyle=':', linewidth=1.5)
    ax2.axvline(support_x, color='gray', linestyle=':', linewidth=1.5)
    
    for L in spans:
        support_x += L
        ax1.axvline(support_x, color='gray', linestyle=':', linewidth=1.5)
        ax2.axvline(support_x, color='gray', linestyle=':', linewidth=1.5)
        
    plt.tight_layout()
    plt.show()
