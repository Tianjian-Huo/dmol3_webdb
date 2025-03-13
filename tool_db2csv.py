import sqlite3
import pandas as pd
import os

def select_db_file():
    """
    让用户选择要转换的数据库文件。
    """
    db_files = [f for f in os.listdir() if f.endswith(".db")]

    if not db_files:
        print("当前目录下没有找到任何 .db 文件！")
        return None

    print("可用的数据库文件:")
    for idx, db_file in enumerate(db_files, start=1):
        print(f"{idx}. {db_file}")

    while True:
        try:
            choice = int(input("\n请输入要转换的数据库文件编号: "))
            if 1 <= choice <= len(db_files):
                return db_files[choice - 1]
            else:
                print("请输入有效的编号！")
        except ValueError:
            print("请输入数字编号！")

def convert_db_to_single_csv(db_file):
    """
    读取 SQLite 数据库文件，并将所有表的数据合并到一个 CSV 文件。
    """
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    # 获取所有表名
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = cursor.fetchall()

    if not tables:
        print(f"数据库 {db_file} 中没有可用的表。")
        conn.close()
        return

    db_name = os.path.splitext(db_file)[0]  # 去掉 .db 后缀
    csv_filename = f"{db_name}_all_tables.csv"

    all_data = []  # 用于存储所有表的数据

    for table in tables:
        table_name = table[0]
        df = pd.read_sql_query(f"SELECT * FROM {table_name}", conn)
        df.insert(0, "table_name", table_name)  # 在第一列添加表名
        all_data.append(df)

    conn.close()

    # **合并所有表的数据**
    final_df = pd.concat(all_data, ignore_index=True)

    # **保存为 CSV**
    final_df.to_csv(csv_filename, index=False, encoding="utf-8-sig")  # 适用于中文
    print(f"\n所有表数据已成功转换为单个 CSV 文件：{csv_filename}")

if __name__ == "__main__":
    db_file = select_db_file()
    if db_file:
        convert_db_to_single_csv(db_file)
