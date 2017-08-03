from sfqm.fqm.genus_symbol import GenusSymbol
from psage.modules.finite_quadratic_module import FiniteQuadraticModule

def assign_genus_symbols(data, limit=None, even=True):
    #only positive definite for now
    if limit == None:
        limit = len(data)
    genus_symbols = []
    for j in range(limit):
        A = data[j]
        try:
            F = FiniteQuadraticModule(A)
            s = F.jordan_decomposition().genus_symbol()
            s = GenusSymbol(s)
            #now determine the 2-type
            tp = twotype(A,s)
            #print [A,tp,str(s)]
            genus_symbols.append([A,tp,str(s)])
        except:
            if even:
                continue
            d = A.det()
            print A, d
            a = crt([2^-1 for p in d.prime_factors() if not p == 2]+[1], [p^valuation(d,p) for p in d.prime_factors() if not p == 2]+[2^(valuation(d,2)+1)])
            print d, a
            A = 2*a*A
            F = FiniteQuadraticModule(A)
            s = F.jordan_decomposition().genus_symbol()
            t = GenusSymbol(s)
            if t._symbol_dict.has_key(2):
                for tt in t._symbol_dict[2]:
                    tt[0] = tt[0]-1
                t._reduce()
                s = str(t)
            tp = twotype(A,s)
            genus_symbols.append([A,tp,str(s)])
    return genus_symbols

def twotype(A,s):
    if s.level() % 2 == 0:
        tp = "I" if s.is_odd() else "II"
    else:
        Q = QuadraticForm(A)
        l = Q.local_genus_symbol(2).symbol_tuple_list()
        tp = "I" if l[0][0][3] == 1 else "II"
    return tp
            

def list_all_genus_symbols(signature=(1,0), det=1, level=1, even=True, test_global=True, unique=True):
    det = Integer(det)
    level = Integer(level)
    #list all genus symbols of given signature for determinant and level up to given bound
    if not even:
        raise NotImplementedError("So far we only support even lattices")
    m,n = signature
    d = m+n #dimension
    asym = []
    for N in range(1,level+1):
        asym = asym + anisotropic_symbols(N, (m-n) % 8)
    sym = filter(lambda x: x.order() <= det, asym)
    #print sym
    print "Starting with {} anisotropic symbols".format(len(asym))
    for p in prime_range(ceil(sqrt(det))):
        print p
        added = True
        while added:
            sym_work = copy(sym)
            added = False
            for s in sym_work:
                if s.order()*p**2 > det:
                    print "Skipping {}".format(s)
                    continue
                C = s.C(p, unique=False)
                for t in C:
                    if t.level() <= level and t.order() <= det:
                        if not t in sym:
                            sym.append(t)
                            added = True
    #print sym
    if test_global:
        sym = filter(lambda x: x.is_global(m,n), sym)
    if unique:
        CU = []
        for s in sym:
            cont = False
            for t in CU:
                if t.defines_isomorphic_module(s):
                    if s<t:
                        CU.remove(t)
                        print "Choosing {} over {}".format(s,t)
                    else:
                        cont = True
                        StopIteration()
            if not cont:
                CU.append(s)
        sym = CU
    return sym

def sort_genus_symbols(symbols):
    ssymbols = {}
    for s in symbols:
        if not ssymbols.has_key(s.order()):
            ssymbols[s.order()] = []
        ssymbols[s.order()].append(s)
    for d in ssymbols.keys():
        ssymbols[d] = sorted(ssymbols[d])
    return ssymbols

def normalize_genus_symbols(data, genus_symbols):
    for d in data:
        s = GenusSymbol(d[2])
        N = s.level()
        D = s.order()
        if not s in genus_symbols[D]:
            print "Need to change symbol for {} (level = {}, det = {})".format(s,s.level(),s.order())
            changed = False
            for t in genus_symbols[D]:
                if t.defines_isomorphic_module(s):
                    d[2] = str(t)
                    print "Changed to {}".format(t)
                    changed = True
                    StopIteration()
            if not changed:
                print "ERROR: No equivalent symbol for {} found!".format(s)
    return data

def number_of_genera(genus_symbols,det,level):
    l = genus_symbols[det]
    return len(s for s in l if l.level() == level)
