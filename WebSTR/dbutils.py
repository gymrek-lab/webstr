import sqlite3

def connect_db(BasePathM):
    """
    Connects to the specific database.
    """
    tfile =  BasePathM + "dbSTR.db"
    conn = sqlite3.connect(tfile)
    return conn
