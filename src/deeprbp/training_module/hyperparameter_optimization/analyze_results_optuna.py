# /src/deeprbp/training_module/hyperparameter_optimization/analyze_results_optuna.py

import os
import argparse
import seaborn as sns
import os
import optuna
import optuna.visualization.matplotlib as optuna_plt
from optuna.storages import RDBStorage
import matplotlib.pyplot as plt

def save_plot(name, output_dir, use_tight_layout=True, remove_legend=False, tight_layout_kwargs=None):
    '''Function to save plot in PNG and PDF'''
    path_png = os.path.join(output_dir, f"{name}.png")
    path_pdf = os.path.join(output_dir, f"{name}.pdf")
    if remove_legend:
        legend = plt.gca().get_legend()
        if legend:
            legend.remove()
    if use_tight_layout:
        try:
            if tight_layout_kwargs is None:
                plt.tight_layout()
            else:
                plt.tight_layout(**tight_layout_kwargs)
        except Exception as e:
            print(f"[Warning] tight_layout() failed: {e}")
    plt.savefig(path_png, dpi=300)
    plt.savefig(path_pdf, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Analyze Optuna study results.")
    parser.add_argument("--storage_path", type=str, required=True, help="Path to Optuna SQLite database")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save output plots")
    args = parser.parse_args()
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Load the study
    study = optuna.load_study(study_name="deeprbp_gridsearch_optuna", storage=RDBStorage(url=f"sqlite:///{args.storage_path}"))

    # Save results to CSV
    study.trials_dataframe().to_csv(os.path.join(args.output_dir, "results_hyp.csv"), index=False)

    # Plotting results
    sns.set_theme(style="white", context="paper", font="DejaVu Sans")
    plt.rcParams.update({
        "axes.edgecolor": "black",
        "axes.linewidth": 0.8,
        "figure.facecolor": "white",   # Fondo blanco fuera del área del plot
        "axes.facecolor": "white",     # Fondo blanco dentro del área del plot
        "savefig.facecolor": "white"   # Fondo blanco en imágenes guardadas
    })

    width_cm = 18
    height_cm = 12
    figsize_inch = (width_cm / 2.54, height_cm / 2.54)

    # 📊 Optimization History
    plt.figure(figsize=figsize_inch)
    optuna_plt.plot_optimization_history(study)
    plt.title("Optimization History", fontsize=14)
    plt.xlabel("Trial", fontsize=12)
    plt.ylabel("Objective Value", fontsize=12)
    save_plot("optimization_history", args.output_dir)

    # 📈 Parameter Importances
    plt.figure(figsize=figsize_inch)
    optuna_plt.plot_param_importances(study)
    plt.xlabel("Importance", fontsize=12)
    save_plot("param_importances", args.output_dir)

    # 📊 Empirical Distribution Function (EDF): Cumulative distribution of objective values
    plt.figure(figsize=figsize_inch)
    optuna_plt.plot_edf(study)
    plt.title("Empirical Distribution Function")
    save_plot("edf", args.output_dir)

    # ⏱️ Timeline Plot: Shows when each trial started and how long it took
    plt.figure(figsize=figsize_inch)
    optuna_plt.plot_timeline(study)
    plt.title("Timeline")
    save_plot("timeline", args.output_dir)

    # 📉 Intermediate Values: Epoch-wise objective values (needs `report_intermediate_values`)
    plt.figure(figsize=figsize_inch)
    optuna_plt.plot_intermediate_values(study)
    plt.title("Intermediate Values")
    save_plot("intermediate_values", output_dir=args.output_dir, remove_legend=True)
    print("✅ Analysis completed! 💾 Results saved to: {args.output_dir}")

if __name__ == "__main__":
    main()
