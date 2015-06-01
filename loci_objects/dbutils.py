# from contextlib import contextmanager
# from sqlalchemy.orm import sessionmaker,session
from sqlalchemy.ext.declarative import declarative_base

'''At the moment this is a stub of a module, which provides the declarative base for all the others.
In the future, it should hold the class which manages the sessions.
'''

dbBase=declarative_base()