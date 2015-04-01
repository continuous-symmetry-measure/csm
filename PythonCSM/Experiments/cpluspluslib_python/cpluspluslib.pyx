__author__ = 'zmbq'

cimport library


def HelloWorld():
    library.HelloWorld()

def Add(a, b):
    return library.Add(a * 2, b) / 2

def AddList(lst):
    return library.AddList(lst)

def Print(msg):
    library.Print(msg.encode('UTF8'))

def GetName(first, last):
    result = library.GetName(first.encode('UTF8'), last.encode('UTF8'))
    return result.decode('UTF8')

def Concat(strings):
    encoded = [s.encode('UTF8') for s in strings]
    result = library.Concat(encoded)
    return result.decode('UTF8')