import matplotlib.pyplot as plt

with open("metrics.log", "r") as f:
    lines = f.readlines()

# Filter and parse lines into a list of metric dictionaries
parsed_lines = [line.strip().split() for line in lines if "±" in line]

# --- Simplification Start ---

# 1. Group metrics by iteration
metric_names = [line[0] for line in parsed_lines]
num_unique_metrics = len(set(metric_names))

iterations_data = []
# Assume metrics appear in a consistent order for each iteration
for i in range(0, len(parsed_lines), num_unique_metrics):
    iteration_chunk = parsed_lines[i:i + num_unique_metrics]
    # Convert the chunk of metrics into a dictionary for that iteration
    iteration_dict = {
        parts[0]: {"mean": float(parts[1]), "std": float(parts[3])}
        for parts in iteration_chunk
    }
    iterations_data.append(iteration_dict)

# 2. Sort the list of iterations by accuracy_score
sorted_iterations = sorted(
    iterations_data,
    key=lambda x: x["mutual_info_proba_classif"]["mean"],
    reverse=True
)

# 3. Reorganize the sorted data for plotting
# The plot function expects a dictionary of lists, not a list of dictionaries
metrics_for_plot = {key: {"mean": [], "std": []} for key in metric_names}
for iteration in sorted_iterations:
    for metric_name, values in iteration.items():
        metrics_for_plot[metric_name]["mean"].append(values["mean"])
        metrics_for_plot[metric_name]["std"].append(values["std"])

# --- Simplification End ---

def plot_metrics(metrics):
    fig, ax = plt.subplots(figsize=(10, 6))
    x_values = range(len(next(iter(metrics.values()))["mean"]))
    
    for metric, values in metrics.items():
        ax.errorbar(
            x_values,
            values["mean"],
            yerr=values["std"],
            label=metric,
            capsize=5,
            marker='o'
        )
    
    ax.set_title("Metrics Comparison")
    ax.set_xlabel("Iterations (Sorted by ROC AUC Score)")
    ax.set_ylabel("Metric Value")
    ax.set_xticks(x_values) # Optional: keeps integer ticks
    ax.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    plt.show()

plot_metrics(metrics_for_plot)