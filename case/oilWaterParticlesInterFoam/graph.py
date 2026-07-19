import os
import numpy as np
import matplotlib.pyplot as plt

CLOUD_NAME = "kinematicCloud"

def get_time_dirs(case_dir="."):
    """Return sorted list of numeric time directories."""
    time_dirs = []
    for name in os.listdir(case_dir):
        path = os.path.join(case_dir, name)
        if os.path.isdir(path):
            try:
                float(name)
                time_dirs.append(name)
            except ValueError:
                continue
    return sorted(time_dirs, key=float)

def read_we_values(we_path):
    """Read scalar values from an OpenFOAM lagrangian field file We."""
    values = []
    with open(we_path, "r") as f:
        in_list = False
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s == "(":
                in_list = True
                continue
            if s == ")":
                break
            if in_list:
                try:
                    values.append(float(s))
                except ValueError:
                    pass
    return values

def main():
    times = []
    we_mean = []

    for t in get_time_dirs("."):
        we_file = os.path.join(t, "lagrangian", CLOUD_NAME, "We")
        if not os.path.exists(we_file):
            continue

        values = read_we_values(we_file)
        if not values:
            continue

        times.append(float(t))
        we_mean.append(float(np.mean(values)))

    if not times:
        print("No We files found. Check that */lagrangian/kinematicCloud/We exists.")
        return

    plt.figure(figsize=(10, 6))
    plt.plot(times, we_mean, marker="o", linestyle="-", color="b", markersize=4)
    plt.xlabel("Time")
    plt.ylabel("Average Weber number")
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
