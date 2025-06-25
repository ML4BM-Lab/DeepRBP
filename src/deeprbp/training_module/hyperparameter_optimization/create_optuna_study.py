# create_study.py
import os
import argparse
import optuna
from optuna.storages import RDBStorage
from optuna.samplers import TPESampler
from optuna.pruners import MedianPruner

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
    sampler = TPESampler(seed=42)
    pruner = MedianPruner(n_startup_trials=5, n_warmup_steps=30, interval_steps=10) #parameter values used by optuna in their tutorials
    
    print(f"🧪 Creating the study: 'deeprbp_gridsearch_optuna' at '{args.output_dir}'...")
    optuna.create_study(
        study_name="deeprbp_gridsearch_optuna",
        direction="minimize",
        pruner=pruner,
        sampler=sampler,
        storage=storage)
    
    print("✅ Optuna study created successfully.")

