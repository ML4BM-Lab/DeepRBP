import sqlite3
import os
import datetime
import time

class BackupManager:
    def __init__(self, db_path, backup_dir, backup_interval_seconds=18000):
        self.db_path = db_path
        self.backup_dir = backup_dir
        self.backup_interval_seconds = backup_interval_seconds
        self.last_backup_time = 0

    def backup_database(self):
        os.makedirs(self.backup_dir, exist_ok=True)
        conn = sqlite3.connect(self.db_path)
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        backup_path = os.path.join(self.backup_dir, f"optuna_backup_{timestamp}.db")
        backup_conn = sqlite3.connect(backup_path)
        with backup_conn:
            conn.backup(backup_conn)
        print(f"Backup creado en: {backup_path}")
        backup_conn.close()
        conn.close()

    def should_backup(self):
        current_time = time.time()
        return (current_time - self.last_backup_time) >= self.backup_interval_seconds

    def update_backup_time(self):
        self.last_backup_time = time.time()
