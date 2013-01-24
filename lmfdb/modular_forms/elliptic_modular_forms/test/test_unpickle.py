import pymongo, gridfs
from sage.all import loads
dbport = 37010
C = pymongo.connection.Connection(port=int(37010))

def test_unpickle():
    r"""
        We test whether the ModularSymbols spaces
        stored in the database unpickle and that
        that :trac:`13998` is fixed
        (patch needs to be applied).
    """
    modforms=C['modularforms']
    modular_symbols = modforms['Modular_symbols']
    files = modular_symbols.files
    finds = files.find()
    n=finds.count()
    print n, ' records found.'
    i=1
    for rec in finds:
        fid = rec['_id']
        fs = gridfs.GridFS(C['modularforms'], 'Modular_symbols')
        f = fs.get(fid)
        M = loads(f.read())
        A = M.ambient_hecke_module()
        if hasattr(A,"_HeckeModule_free_module__decomposition"):
            S = M.cuspidal_submodule().decomposition()[0].free_module().basis()
            try:
                print hash(S)
                print 'Successfully tested hash ', i , '/', n
                i=i+1
            except:
                'Hashing failed for' + M, S
                return M
            
