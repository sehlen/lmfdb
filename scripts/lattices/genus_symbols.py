from sfqm.fqm.genus_symbol import GenusSymbol
from psage.modules.finite_quadratic_module import FiniteQuadraticModule
from itertools import product
import re

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
            

def list_genus_symbols(signature=(1,0), det=1, level=1, even=True, test_global=True, unique=True):
    det = Integer(det)
    level = Integer(level)
    #list all genus symbols of given signature for determinant and level as given
    if not even:
        raise NotImplementedError("So far we only support even lattices")
    m,n = signature
    d = m+n #dimension
    sym = []
    psyms = {}
    for p in det.prime_factors():
        psyms[p] = []
        v = det.valuation(p)
        l = Partitions(v)
#        print p, v
        from collections import Counter
        psyms[p] = []
        for a in l:
            c = Counter(a)
            if p != 2:
                s = ".".join([str(p**e) + "^" + "{}" + str(r) for e,r in c.iteritems()])
                for eps in product(["+","-"], repeat=len(c)):
                    psyms[p].append(s.format(*eps))
            if p == 2:
                s = ".".join([str(p**e) + "{}" + "^{}" + str(r) for e,r in c.iteritems()])
                fillings = product(product(['','_0','_1','_2','_3','_4','_5','_6','_7'],["+","-"]), repeat=len(c))
                fillings = map(lambda x: flatten(x), fillings)
#                print fillings
                for teps in fillings:
                        ss = s.format(*teps)
#                        print ss
                        psyms[p].append(ss)
#        print psyms
    for s in product(*psyms.values()):
        #print '.'.join(s)
        try:
 #           print s
            t = GenusSymbol('.'.join(s))
            if t.level() == level:
                sym.append(t)
#                print "HERE", t
#                print sym
            if t.p_rank(2) < m+n:
                t1 = deepcopy(t)
                if not t1._symbol_dict.has_key(2):
                    t1._symbol_dict[2] = []
                a = prod(x[2] for x in t1._symbol_dict[2])
                b = kronecker(odd_part(det),2)
                n2 = m+n-t.p_rank(2)
                o = (m-n-t.signature()) % 8
 #               print a,b,n2,o
                if n2 == 2 :
                    if a*b == -1:
                        if is_odd(o) or o == 0:
                            continue
                    if a*b == 1:
                        if is_odd(o) or o == 4:
                            continue
                t1._symbol_dict[2].insert(0, [0,n2,a*b,1,o])
                if t1.level() == level:
                    sym.append(t1)
        except ValueError:
            continue

    if test_global:
        sym = filter(lambda x: x.is_global(m,n), sym)
    if unique:
        CU = []
        for s in sym:
            cont = False
            for t in CU:
                if t.defines_isomorphic_module(s, ignore_unimodular_2_component=False):
                    if s<t:
                        CU.remove(t)
                        print "Choosing {} over {}".format(s,t)
                    else:
                        cont = True
                        StopIteration()
            if not cont:
                CU.append(s)
        sym = uniq(CU)
    return sym

@parallel
def list_genus_symbols_range(signature=(1,0), det_range=range(1,10)):
    if det_range in ZZ:
        det_range=range(1,det_range)
    ls = {}
    for d in det_range:
        ls[d] = []
        for N in Integer(2*d).divisors():
            ldn=list_genus_symbols(signature,d,N)
            #print ldn
            if len(ldn)>0:
                ls[d] = ls[d] + ldn
    return ls

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
        s = re.findall('\(([^\)]*)\)',d[2])
        if len(s) == 0:
            s = d[2]
        else:
            s = s[0]
        s = GenusSymbol(s)
        N = s.level()
        D = s.order()
        if not genus_symbols.has_key(D):
            print "No symbols of oder {} in list".format(D)
            continue
        if not s in genus_symbols[D]:
            print "Need to change symbol for {} (level = {}, det = {})".format(s,s.level(),s.order())
            changed = False
            for t in genus_symbols[D]:
                if t.defines_isomorphic_module(s, ignore_unimodular_2_component=False):
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
