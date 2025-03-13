import sqlite3
import os

def last_id(db_path: str) -> int:
    from ase.db import connect
    db = connect(db_path)
    for row in db.select():
        last_id = row.id
    return last_id

current_path = os.path.dirname(os.path.abspath(__file__))
file_1 = os.path.join(current_path, 'DATABASE.db')  # 目标数据库
file_2 = os.path.join(current_path, "DMOL_RESULTS.db")  # 被合并的数据库

conn1 = sqlite3.connect(file_1)
conn2 = sqlite3.connect(file_2)
c1 = conn1.cursor()
c2 = conn2.cursor()

c1.execute('SELECT id from systems')
row_num1 = len(c1.fetchall())  
db1_last_id = last_id(file_1)  

c2.execute('SELECT id from systems')
row_num2 = len(c2.fetchall())

print(f'初始 DATABASE.db 行数: {row_num1}, 被合并 DMOL_RESULTS.db 行数: {row_num2}')

for i in range(1, row_num2+1):
    c2.execute("UPDATE systems set ID=? where ID=?", (10000 + i, i))
    c2.execute("UPDATE keys set id=? where id=?", (10000 + i, i))
    c2.execute("UPDATE species set id=? where id=?", (10000 + i, i))
    c2.execute("UPDATE text_key_values set id=? where id=?", (10000 + i, i))
    c2.execute("UPDATE number_key_values set id=? where id=?", (10000 + i, i))
conn2.commit()

for j in range(1, row_num2+1):
    new_id = db1_last_id + j
    c2.execute("UPDATE systems set ID=? where ID=?", (new_id, 10000 + j))
    c2.execute("UPDATE keys set id=? where id=?", (new_id, 10000 + j))
    c2.execute("UPDATE species set id=? where id=?", (new_id, 10000 + j))
    c2.execute("UPDATE text_key_values set id=? where id=?", (new_id, 10000 + j))
    c2.execute("UPDATE number_key_values set id=? where id=?", (new_id, 10000 + j))
conn2.commit()
conn2.close()

c1.execute(f"attach '{file_2}' as SecondaryDB")
c1.execute('insert into systems select * from SecondaryDB.systems')
c1.execute('insert into species select * from SecondaryDB.species')
c1.execute('insert into keys select * from SecondaryDB.keys')
c1.execute('insert into text_key_values select * from SecondaryDB.text_key_values')
c1.execute('insert into number_key_values select * from SecondaryDB.number_key_values')

conn1.commit()
c1.execute('detach SecondaryDB')
conn1.close()
print("合并 DMOL 结果到 DATABASE.db 完成！")