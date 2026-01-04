from pathlib import Path
import subprocess
import numpy as np
import matplotlib.pyplot as plt

# Path to ngspice executable (edit if needed on another machine)
NGSPICE_EXE = "ngspice"
# NOTE: If ngspice is not in PATH, replace with full path to ngspice executable

HERE = Path(__file__).resolve().parent
SUBCKT = HERE / "mem_simmons_pickett.subckt"


def run_ngspice(netlist_text: str, workname: str):
    cir = HERE / f"{workname}.cir"
    log = HERE / f"{workname}.log"

    cir.write_text(netlist_text, encoding="utf-8")

    cmd = [NGSPICE_EXE, "-b", "-o", str(log), str(cir)]
    print("Running:", " ".join(cmd))

    p = subprocess.run(cmd, cwd=str(HERE))
    if p.returncode != 0:
        print("\n--- NGSPICE FAILED. Showing log ---\n")
        if log.exists():
            print(log.read_text(errors="ignore"))
        raise SystemExit(p.returncode)


def load_wrdata(path: Path):
    return np.loadtxt(path)


def make_iv():
    out = HERE / "iv.csv"
    T = 10e-3
    Vmax = 1.2

    netlist = f"""
* IV test (triangular stimulus)
.include "{SUBCKT.as_posix()}"

* instance (override params explicitly)
Xmem te 0 w dwdt MEM xinit=0.2 Ron=200 Roff=200k k=2e3 p=2

Vin te 0 PWL( 0 0  {T/4} {Vmax}  {3*T/4} {-Vmax}  {T} 0 )

.tran 10u {T}

.control
set filetype=ascii
run
* IMPORTANT: write to local file name (Windows ngspice hates long paths)
wrdata iv.csv time v(te) i(vin) v(w) v(dwdt)
quit
.endc

.end
"""
    run_ngspice(netlist, "iv_run")

    if not out.exists():
        print("ERROR: ngspice did not create iv.csv. Check iv_run.log above.")
        raise SystemExit(1)

    d = load_wrdata(out)
    time = d[:, 0]
    v = d[:, 1]
    i = d[:, 2]
    w = d[:, 3]

    plt.figure()
    plt.plot(v, -i)
    plt.xlabel("V (V)")
    plt.ylabel("I (A)")
    plt.title("I–V (triangular stimulus)")
    plt.grid(True)

    plt.figure()
    plt.plot(time, w)
    plt.xlabel("t (s)")
    plt.ylabel("state x (stored as V(w))")
    plt.title("State x(t)")
    plt.grid(True)


def make_pulses():
    out = HERE / "pulses.csv"

    Vp = 1.5
    pw = 200e-6
    gap = 200e-6
    n_pairs = 20

    t = 0.0
    pts = [(t, 0.0)]
    for _ in range(n_pairs):
        # + pulse
        pts += [(t, 0.0), (t + 1e-9, Vp), (t + pw, Vp), (t + pw + 1e-9, 0.0)]
        t += pw + gap
        # - pulse
        pts += [(t, 0.0), (t + 1e-9, -Vp), (t + pw, -Vp), (t + pw + 1e-9, 0.0)]
        t += pw + gap
    pts.append((t, 0.0))

    pwl = " ".join([f"{tt} {vv}" for tt, vv in pts])

    netlist = f"""
* Pulse test
.include "{SUBCKT.as_posix()}"

Xmem te 0 w dwdt MEM xinit=0.2 Ron=200 Roff=200k k=2e3 p=2

Vin te 0 PWL( {pwl} )

.tran 5u {t}

.control
set filetype=ascii
run
wrdata pulses.csv time v(te) i(vin) v(w) v(dwdt)
quit
.endc

.end
"""
    run_ngspice(netlist, "pulse_run")

    if not out.exists():
        print("ERROR: ngspice did not create pulses.csv. Check pulse_run.log above.")
        raise SystemExit(1)

    d = load_wrdata(out)
    time = d[:, 0]
    w = d[:, 3]
    dwdt = d[:, 4]

    # Sample end-of-pulses -> effective resistance trend
    sample_times = []
    tt = 0.0
    for _ in range(n_pairs):
        tt += pw
        sample_times.append(tt)
        tt += gap + pw
        sample_times.append(tt)
        tt += gap

    sample_times = np.array(sample_times)
    idx = np.searchsorted(time, sample_times, side="left")
    idx = np.clip(idx, 0, len(time) - 1)
    w_s = w[idx]

    # same linear mapping as subckt: R = Ron*x + Roff*(1-x)
    Ron, Roff = 200.0, 200e3
    xcl = np.minimum(np.maximum(w_s, 0.0), 1.0)
    R_s = Ron * xcl + Roff * (1.0 - xcl)

    plt.figure()
    plt.plot(np.arange(len(R_s)), R_s, marker="o")
    plt.xlabel("Pulse index")
    plt.ylabel("R_eff (Ohm)")
    plt.title("Potentiation / Depression (effective R vs pulse)")
    plt.grid(True)

    plt.figure()
    plt.plot(w, dwdt, ".", markersize=2)
    plt.xlabel("x (stored as V(w))")
    plt.ylabel("dx/dt (V/s)")
    plt.title("Dynamic route map: dx/dt vs x")
    plt.grid(True)


if __name__ == "__main__":
    make_iv()
    make_pulses()

    # save all figures
    outdir = HERE / "results"
    outdir.mkdir(exist_ok=True)
    for i, num in enumerate(plt.get_fignums(), start=1):
        plt.figure(num)
        plt.savefig(outdir / f"fig_{i}.png", dpi=200, bbox_inches="tight")

    plt.show()
