"""
Memristor Characterization with PySpice (Pickett/Simmons tunnel-barrier model)
----------------------------------------------------------------------------

Produces three standard experiments:
  1) I–V curve under a triangular stimulus
  2) Potentiation/Depression via write–read pulse trains
  3) Dynamic Route Map (DRM): |dw/dt| vs w for multiple DC voltages

This script is designed to be numerically robust in ngspice by:
  - Using a "safe-math" subcircuit (no ln(0), sqrt(neg), divide-by-zero)
  - Forcing strictly increasing PWL time points
  - Using Gear integration and a capped maximum timestep

Files expected in the same folder:
  - mem_pickett_simmons_safe.subckt
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from PySpice.Spice.Netlist import Circuit


HERE = Path(__file__).resolve().parent
SUBCKT_PATH = HERE / "mem_pickett_simmons_safe.subckt"


# -----------------------------
# Utility: strict PWL builder
# -----------------------------
def pwl_from_points(points, eps=1e-12) -> str:
    """
    Build a SPICE PWL(...) string from (t, v) points ensuring strictly increasing time.

    ngspice warns/behaves poorly if any time points are equal or decreasing.
    We enforce monotonicity by nudging non-increasing times by a tiny epsilon.
    """
    out = []
    last_t = -np.inf
    for t, v in points:
        t = float(t)
        if t <= last_t:
            t = last_t + eps
        out.append((t, float(v)))
        last_t = t
    return "PWL(" + " ".join(f"{t:.12e} {v:.6g}" for t, v in out) + ")"


def add_spice_numerics(circuit: Circuit, max_timestep: float | None = None) -> None:
    """
    Add ngspice options helpful for stiff behavioral models.
    """
    opts = [
        ".options method=gear maxord=2",
        ".options reltol=1e-4 abstol=1e-12 vntol=1e-8 trtol=7",
    ]
    if max_timestep is not None:
        # Transient timestep limiting is typically done via .tran ... <tmax>
        # PySpice exposes this as simulator.transient(..., max_time=...)
        # but we also keep this here in case users run the generated netlist directly.
        opts.append(f"* recommended max timestep: {max_timestep:g} s")
    circuit.raw_spice += "\n" + "\n".join(opts) + "\n"


# -----------------------------
# Core simulator wrapper
# -----------------------------
def run_simulation(
    voltage_source: str,
    x0: float,
    end_time: float,
    step_time: float,
    max_timestep: float | None = None,
    input_node: str = "te",
    source_name: str = "in",
):
    circuit = Circuit("PickettSimmonsMemristor")
    circuit.raw_spice += SUBCKT_PATH.read_text()
    add_spice_numerics(circuit, max_timestep=max_timestep)

    circuit.V(source_name, input_node, circuit.gnd, voltage_source)
    circuit.X("m1", "MEM", input_node, circuit.gnd, "x", "dxdt")

    sim = circuit.simulator()
    sim.initial_condition(x=x0)

    # In PySpice/ngspice, max_time limits internal timestep.
    # If None, ngspice chooses adaptively (can miss fast edges).
    if max_timestep is None:
        analysis = sim.transient(step_time=step_time, end_time=end_time)
    else:
        analysis = sim.transient(step_time=step_time, end_time=end_time, max_time=max_timestep)

    return analysis


# -----------------------------
# Experiment 1: I–V
# -----------------------------
def run_iv(x0: float, amplitude: float, period: float, step_time: float, max_timestep: float):
    times = [0, period / 4, period / 2, 3 * period / 4, period]
    volts = [0, amplitude, 0, -amplitude, 0]
    pwl = pwl_from_points(list(zip(times, volts)))

    analysis = run_simulation(
        voltage_source=pwl,
        x0=x0,
        end_time=period,
        step_time=step_time,
        max_timestep=max_timestep,
        input_node="te",
        source_name="in",
    )

    v = np.array(analysis["te"])
    i = -np.array(analysis.branches["vin"])  # current through V(in)
    x = np.array(analysis.nodes["x"])
    dxdt = np.array(analysis.nodes["dxdt"])
    t = np.array(analysis.time)

    return t, v, i, x, dxdt


# -----------------------------
# Experiment 2: Potentiation/Depression (write–read)
# -----------------------------
def build_write_read_train(
    V_write: float,
    T_write: float,
    V_read: float,
    T_read: float,
    t_slope: float,
    num_pulses: int,
):
    """
    Build a single continuous write–read pulse train:
    [ramp->Vw][hold][ramp->0][idle][ramp->Vr][hold][ramp->0][idle]...
    """
    pts = []
    t = 0.0
    pts.append((t, 0.0))

    def ramp(to_v):
        nonlocal t
        pts.append((t, pts[-1][1]))
        t += t_slope
        pts.append((t, to_v))

    def hold(v, dt):
        nonlocal t
        pts.append((t, v))
        t += dt
        pts.append((t, v))

    def idle(dt):
        nonlocal t
        pts.append((t, 0.0))
        t += dt
        pts.append((t, 0.0))

    read_sample_times = []

    gap = t_slope  # small idle between phases; keep simple/consistent

    for _ in range(num_pulses):
        # WRITE
        ramp(V_write)
        hold(V_write, T_write)
        ramp(0.0)
        idle(gap)

        # READ
        ramp(V_read)
        hold(V_read, T_read)
        read_sample_times.append(t)  # sample at end of read hold
        ramp(0.0)
        idle(gap)

    return pwl_from_points(pts), t, np.array(read_sample_times)


def run_pot_dep(
    xoff: float,
    xon: float,
    V_pot: float,
    T_pot: float,
    V_dep: float,
    T_dep: float,
    V_read: float,
    T_read: float,
    t_slope: float,
    num_pulses: int,
    step_time: float,
    max_timestep: float,
):
    # Potentiation train (start from OFF-ish)
    pwl_pot, t_end_pot, t_read_pot = build_write_read_train(V_pot, T_pot, V_read, T_read, t_slope, num_pulses)
    a_pot = run_simulation(pwl_pot, x0=xoff, end_time=t_end_pot, step_time=step_time, max_timestep=max_timestep, source_name="in")

    tp = np.array(a_pot.time)
    idx_p = np.searchsorted(tp, t_read_pot)
    idx_p = np.clip(idx_p, 0, len(tp) - 1)

    i_read_p = -np.array(a_pot.branches["vin"])[idx_p]
    v_read_p = np.array(a_pot["te"])[idx_p]
    g_pot = i_read_p / np.maximum(v_read_p, 1e-12)

    # Depression train (start from ON-ish)
    pwl_dep, t_end_dep, t_read_dep = build_write_read_train(V_dep, T_dep, V_read, T_read, t_slope, num_pulses)
    a_dep = run_simulation(pwl_dep, x0=xon, end_time=t_end_dep, step_time=step_time, max_timestep=max_timestep, source_name="in")

    td = np.array(a_dep.time)
    idx_d = np.searchsorted(td, t_read_dep)
    idx_d = np.clip(idx_d, 0, len(td) - 1)

    i_read_d = -np.array(a_dep.branches["vin"])[idx_d]
    v_read_d = np.array(a_dep["te"])[idx_d]
    g_dep = i_read_d / np.maximum(v_read_d, 1e-12)

    return g_pot, g_dep


# -----------------------------
# Experiment 3: Dynamic Route Map
# -----------------------------
def run_drm(voltages, t_write: float, step_time: float, max_timestep: float, xon: float, xoff: float):
    results = {}
    for V in voltages:
        x0 = xoff if V < 0 else xon
        a = run_simulation(f"{V}", x0=x0, end_time=t_write, step_time=step_time, max_timestep=max_timestep, source_name="in")
        x = np.array(a.nodes["x"])
        dxdt = np.array(a.nodes["dxdt"])
        results[float(V)] = (x, dxdt)
    return results


# -----------------------------
# Main
# -----------------------------
def main():
    # Recommended state bounds from the Pickett-style fit (meters)
    xon = 1.3e-9   # near a_off
    xoff = 1.7e-9  # near a_on

    # Global numeric settings (important)
    max_timestep = 1e-6   # cap internal timestep (adjust if you use faster edges)
    step_time_iv = 1e-4
    step_time_pulse = 1e-6
    step_time_drm = 1e-3

    # ---- Plot 1: I–V ----
    t, v, i, x, dxdt = run_iv(x0=xon, amplitude=1.45, period=1.0, step_time=step_time_iv, max_timestep=max_timestep)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))

    ax1.set_title("I–V (triangular)")
    ax1.plot(v, i, ".-", linewidth=1)
    ax1.set_xlabel("V (V)")
    ax1.set_ylabel("I (A)")
    ax1.grid(True)

    # ---- Plot 2: State vs time ----
    ax2.set_title("State x(t) = w(t)")
    ax2.plot(t, x * 1e9, "-", linewidth=1)
    ax2.set_xlabel("t (s)")
    ax2.set_ylabel("w (nm)")
    ax2.grid(True)

    # ---- Plot 3: Pot/Dep ----
    V_pot, T_pot = -1.15, 10e-3
    V_dep, T_dep = +1.2, 100e-3
    V_read, T_read = 0.1, 10e-6
    t_slope = 1e-6
    num_pulses = 40

    g_pot, g_dep = run_pot_dep(
        xoff=xoff,
        xon=xon,
        V_pot=V_pot,
        T_pot=T_pot,
        V_dep=V_dep,
        T_dep=T_dep,
        V_read=V_read,
        T_read=T_read,
        t_slope=t_slope,
        num_pulses=num_pulses,
        step_time=step_time_pulse,
        max_timestep=max_timestep,
    )

    ax3.set_title("Potentiation / Depression (read conductance)")
    ax3.plot(np.arange(len(g_pot)), g_pot * 1e6, ".-")
    ax3.plot(len(g_pot) + np.arange(len(g_dep)), g_dep * 1e6, ".-")
    ax3.set_xlabel("Pulse index")
    ax3.set_ylabel("G_read (µS)")
    ax3.grid(True)

    plt.tight_layout()

    # ---- DRM figure ----
    fig2, ax = plt.subplots(1, 1, figsize=(6, 4))
    drm = run_drm(np.arange(-1.2, 1.3, 0.4), t_write=2.0, step_time=step_time_drm, max_timestep=max_timestep, xon=xon, xoff=xoff)
    for V, (xv, dxv) in drm.items():
        ls = "-" if V > 0 else "--"
        ax.plot(xv * 1e9, np.abs(dxv), ls, label=f"V={V:.1f}V")
    ax.set_yscale("log")
    ax.set_xlabel("w (nm)")
    ax.set_ylabel("|dw/dt| (m/s)")
    ax.set_title("Dynamic Route Map (DRM)")
    ax.grid(True, which="both")
    ax.legend()

    # Save artifacts
    out = HERE / "results"
    out.mkdir(exist_ok=True)
    fig.savefig(out / "fig_main.png", dpi=200, bbox_inches="tight")
    fig2.savefig(out / "fig_drm.png", dpi=200, bbox_inches="tight")

    np.savetxt(out / "iv_trace.csv", np.column_stack([t, v, i, x, dxdt]),
               delimiter=",", header="t,V,I,w,dwdt", comments="")
    np.savetxt(out / "pot_dep.csv", np.column_stack([np.arange(len(g_pot)), g_pot, g_dep]),
               delimiter=",", header="k,G_pot,G_dep", comments="")

    plt.show()


if __name__ == "__main__":
    main()
