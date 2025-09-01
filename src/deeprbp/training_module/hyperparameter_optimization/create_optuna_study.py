# create_study.py
import os
import argparse
import optuna
from optuna.storages import RDBStorage
from optuna.samplers import TPESampler
from optuna.pruners import MedianPruner

# def get_sampler_and_pruner(): # used for first 80 trials (fase 1 broad)
#     sampler = TPESampler(
#                     n_startup_trials=75, # Exploración inicial
#                     multivariate=True,
#                     group=True # por nuestro espacio jerárquico
#     )
#     pruner = MedianPruner(
#         n_startup_trials=25, # Espera 25 trials antes de podar
#         n_warmup_steps=470, # Espera alrededor de 10 epochs antes de evaluar
#         interval_steps=47   # Revisa cada 1 epoch promedio
#     )
#     print(f"🔍 Sampler configuration: n_startup_trials=75, multivariate=True, group=True")
#     print(f"🔍 Pruner configuration: n_startup_trials=20, n_warmup_steps=470, interval_steps=47")
#     return sampler, pruner

def get_sampler_and_pruner(seed: int = 42):
    sampler = TPESampler(
        n_startup_trials=10,       # pocas aleatorias antes de TPE “pleno”
        multivariate=True,
        group=True,
        n_ei_candidates=64,
        seed=seed,
    )
    pruner = MedianPruner(
        n_startup_trials=0,        # permitir podas desde el principio de fase 2
        n_warmup_steps=3,          # epochs
        interval_steps=1,
    )
    print("🔍 Sampler (phase 2): n_startup_trials=10, multivariate=True, group=True, n_ei_candidates=64")
    print("🔍 Pruner  (phase 2): n_startup_trials=0, n_warmup_steps=3, interval_steps=1 (epoch-based)")
    return sampler, pruner

def main():
    parser = argparse.ArgumentParser(description="Setup step for Optuna optimization: create an Optuna study using a median pruner and TPE sampler.")
    parser.add_argument("--output_dir",type=str, required=True, help="Directory where the Optuna storage (optuna.db) will be saved.")
    args = parser.parse_args()
    
    # 👉 Option 1: SQLite (local file)
    os.makedirs(args.output_dir, exist_ok=True)
    db_path = f"sqlite:///{os.path.join(args.output_dir, 'optuna.db')}"
    storage = RDBStorage(url=db_path)
    
    # 👉 Opción 2: PostgreSQL (descomenta si usas PostgreSQL)
    # # storage = RDBStorage(url="postgresql://optuna_user:supersecurepassword@your-db-host:5432/optuna_db")
    sampler, pruner = get_sampler_and_pruner()

    print(f"🧪 Creating the study: 'deeprbp_gridsearch_optuna' at '{args.output_dir}'...")
    optuna.create_study(
        study_name="deeprbp_gridsearch_optuna",
        direction="minimize",
        pruner=pruner,
        sampler=sampler,
        storage=storage
    )

    print("✅ Optuna study created successfully.")

if __name__ == "__main__":
    main()