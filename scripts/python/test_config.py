from pathlib import Path
import sqlite3
import hashlib
import json
import subprocess

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DB_FILE = PROJECT_ROOT / "data/database/screening_database.sqlite"

# Current 1IEP docking configuration
config = {
    "target": "1IEP_chainA",
    "center_x": 15.190,
    "center_y": 53.903,
    "center_z": 16.917,
    "size_x": 19.0,
    "size_y": 27.0,
    "size_z": 24.0,
    "exhaustiveness": 8,
    "num_modes": 10
}

# Get the installed Vina version
result = subprocess.run(
    ["vina", "--version"],
    capture_output=True,
    text=True
)

vina_version = result.stdout.strip() or result.stderr.strip()

config["vina_version"] = vina_version

# Create a stable representation of the configuration
config_text = json.dumps(config, sort_keys=True)

config_hash = hashlib.sha256(
    config_text.encode("utf-8")
).hexdigest()

connection = sqlite3.connect(DB_FILE)
cursor = connection.cursor()

# Target ID 1 was created in the previous step
cursor.execute("""
INSERT INTO docking_configs (
    target_id,
    center_x,
    center_y,
    center_z,
    size_x,
    size_y,
    size_z,
    exhaustiveness,
    num_modes,
    vina_version,
    config_hash
)
VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
""", (
    1,
    config["center_x"],
    config["center_y"],
    config["center_z"],
    config["size_x"],
    config["size_y"],
    config["size_z"],
    config["exhaustiveness"],
    config["num_modes"],
    vina_version,
    config_hash
))

config_id = cursor.lastrowid
connection.commit()

print(f"Configuration inserted: ID={config_id}")
print(f"Vina version: {vina_version}")
print(f"Configuration hash: {config_hash}")

print("\nDocking configurations:")
for row in cursor.execute("""
    SELECT
        config_id,
        target_id,
        center_x,
        center_y,
        center_z,
        size_x,
        size_y,
        size_z,
        exhaustiveness,
        num_modes,
        vina_version
    FROM docking_configs
"""):
    print(row)

connection.close()
