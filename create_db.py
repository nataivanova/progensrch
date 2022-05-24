#!/usr/bin/env python3

import pickle
import numpy as np
import sqlite3


cur.execute("drop table if exists tracks")
cur.execute('''CREATE TABLE tracks (mdonor real, maccretor real, mtrate real, period real, teff real, age real, radius real, dt real, srcfile text)''')
cur.execute("create index if not exists idx_mdonor ON tracks(mdonor)")
cur.execute("create index if not exists idx_maccretor ON tracks(maccretor)")
cur.execute("create index if not exists idx_mtrate ON tracks(mtrate)")
cur.execute("create index if not exists idx_period ON tracks(period)")
cur.execute("create index if not exists idx_teff ON tracks(teff)")
cur.execute("create index if not exists idx_age ON tracks(age)")
cur.execute("create index if not exists idx_radius ON tracks(radius)")
cur.execute("create index if not exists idx_dt ON tracks(dt)")
cur.execute("create index if not exists idx_srcfile ON tracks(srcfile)")
for fname in allfiles:
    with open(fname.strip(),'rb') as _:
        v = pickle.load(_)
        u = np.transpose(v)
    for _ in u:
        _ = list(_)
        _.append(fname.strip())
        cur.execute("insert into tracks values (?,?,?,?,?,?,?,?,?)", _)
        con.commit()
