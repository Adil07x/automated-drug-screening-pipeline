import sqlite3

DB_PATH = "data/database/screening_database.sqlite"

NEW_COLUMNS = {
    "tpsa": "REAL",
    "rotatable_bonds": "INTEGER",
    "pains_alert_count": "INTEGER",
    "pains_alerts": "TEXT",
    "brenk_alert_count": "INTEGER",
    "brenk_alerts": "TEXT",
    "sa_score": "REAL",
    "sc_score": "REAL",
}

conn = sqlite3.connect(DB_PATH)
cursor = conn.cursor()

existing_columns = {
    row[1]
    for row in cursor.execute("PRAGMA table_info(compound_properties)")
}

for column, data_type in NEW_COLUMNS.items():
    if column not in existing_columns:
        cursor.execute(
            f"ALTER TABLE compound_properties ADD COLUMN {column} {data_type}"
        )
        print(f"Added column: {column}")
    else:
        print(f"Already exists: {column}")

conn.commit()
conn.close()

print("\nDatabase schema update completed.")
