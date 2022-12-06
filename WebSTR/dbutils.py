import pandas as pd
import sqlite3

def run_query_withparms(sql, BasePath, glochrom):
    conn = sqlite3.connect(BasePath + "dbSTR" + str(glochrom) + ".db")
    df   = pd.read_sql_query(sql, conn)
    return df

def get_db(BasePathM):
    """
    Opens a new database connection if there is none yet for the
    current application context.
    """
    db_conn = connect_db(BasePathM)
    return db_conn

def connect_db(BasePathM):
    """
    Connects to the specific database.
    """
    tfile =  BasePathM + "dbSTR.db"
    print(tfile)
    conn = sqlite3.connect(tfile)
    return conn
