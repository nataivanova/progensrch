#!/usr/bin/env python3

import pickle
import numpy as np
import sqlite3 as sql

con = sql.connect('/var/tmp/progen_tool/db.sqlite')

cur = con.cursor()

with open('/tmp/alldbfiles','r') as f:
    allfiles = f.readlines()

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

cur.execute(
    "CREATE TABLE meta AS
    SELECT CASE WHEN srcfile LIKE '%/ns_data/%' THEN FALSE ELSE TRUE END as isbh
    , min(mdonor) as mdonor_min
    , max(mdonor) as mdonor_max
    , min(maccretor) as maccretor_min
    , max(maccretor) as maccretor_max
    , srcfile FROM tracks group by srcfile"
)
cur.execute("create index if not exists idx_meta_srcfile       ON meta(srcfile)")
cur.execute("create index if not exists idx_meta_maccretor_max ON meta(maccretor_max)")
cur.execute("create index if not exists idx_meta_maccretor_min ON meta(maccretor_min)")
cur.execute("create index if not exists idx_meta_donor_max     ON meta(donor_max)")
cur.execute("create index if not exists idx_meta_donor_min     ON meta(donor_min)")
