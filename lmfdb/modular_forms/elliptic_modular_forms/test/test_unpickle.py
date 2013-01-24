import pymongo
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
    files = C[db_name][collection].files
    finds = files.find()
    n=finds.count()
    i=1
    for rec in finds:
        fid = rec['_id']
        fs = gridfs.GridFS(C[db_name], collection)
        f = fs.get(fid)
        M = loads(f.read())
        S = M.cuspidal_submodule().decomposition()[0].free_module().basis()
        try:
            print hash(S)
            print 'Successfully tested hash ' + i + '/' + n
            i=i+1
        except:
            'Hashing failed for' + M, S
            return M
            
