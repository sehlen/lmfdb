
from lmfdb.base import LmfdbTest, getDBConnection

from flask import request
import unittest2

from views.emf_main import *
from . import emf_logger
emf_logger.setLevel(100)

def anum(a):
    return sum([26**i*(ord(a[len(a)-i-1])-ord('a')+1) for i in range(len(a))])

def check_orbit_list(a):
    return set([anum(c) for c in a]) == set(range(1,len(a)+1))

class EmfTest(LmfdbTest):

    def runTest():
        pass

    def test_gamma0_pages(self):
        version=1.3
        wmax=40
        Nmax=24
        errors = []
        spaces = getDBConnection().modularforms2.webmodformspace
        forms = getDBConnection().modularforms2.webnewforms
        data = spaces.find({'weight':{'$ge':int(2)},'weight':{'$lt':int(wmax+1)},'level':{'$lt':int(Nmax+1)},'character':int(1),'version':float(version)})
        print "Checking %d spaces with trivial character of weight w <= %d and level N <= %d"%(data.count(),wmax,Nmax)
        for s in data:
            if s['space_label'] != "%d.%d.%d"%(s['level'],s['weight'],s['character']):
                 print "Label %s does not match level=%d, weight=%d, character=%d"%(s['space_label'],s['level'],s['weight'],s['character'])
                 errors.append(s['space_label'])
            if not check_orbit_list(s['hecke_orbits']):
                print "Space %s has a bad list of Hecke orbits: %s" % (s['space_label'], s['hecke_orbits'])
                errors.append(s['space_label'])
            orbits = forms.find({'version':float(version),'parent':s['space_label']})
            olabels = [r['label'] for r in orbits]
            if len(olabels) != len(s['hecke_orbits']) or set(olabels) != set(s['hecke_orbits']):
                print "Hecke orbit data in webnewforms for space %s is incomplete or inconsistent" % label
                print "    %s versus %s" % (olabels,stab[label][0])
                errors.append(s['space_label'])
            orbits = orbits.rewind()
            odims = [r['dimension'] for r in orbits] 
            if sum(odims) != s['dimension_new_cusp_forms']:
                print "Hecke orbit dimensions %s do not sum to %d for space %s" % (odims, s['dimension_new_cusp_forms'], s['space_label'])
                errors.append(s['space_label'])
            l = s['space_label'].split('.')
            url = "ModularForm/GL2/Q/holomorphic/%s/%s/%s/"%(l[0],l[1],l[2])
            print "Checking %s (%d Hecke orbits, total dimension %d)"%(url,len(s['hecke_orbits']),s['dimension_new_cusp_forms'])
            page = self.tc.get(url, follow_redirects=True)
            if s['dimension_new_cusp_forms'] == 0:
                if not "no newforms of this weight, level and character" in page.data:
                    print "Failed on", url
                    errors.append(s['space_label'])
            else:
                if not ("Decomposition of" in page.data and "irreducible Hecke orbits" in page.data):
                    print "Failed on", url
                    errors.append(s['space_label'])
                orbits.rewind()
                for r in orbits:
                    if not r['hecke_orbit_label'] in page.data:
                        print "Hecke orbit label %s does not appear on page %s"%(r['hecke_orbit_label'],url)
                        errors.append(r['hecke_orbit_label'])
                orbits.rewind()
                for r in orbits:
                    if not 'hecke_orbit_label' in r:
                        print "no hecke_orbit_label in ", r
                        errors.append(r['hecke_orbit_label'])
                    l = r['hecke_orbit_label'].split('.')
                    if len(l) != 4:
                        print 'bad hecke_orbit_label', l
                        errors.append(r['hecke_orbit_label'])
                    url = "ModularForm/GL2/Q/holomorphic/%s/%s/%s/%s/"%(l[0],l[1],l[2],l[3])
                    print "Checking %s"%url
                    page = self.tc.get(url, follow_redirects=True)
                    if not "Fourier coefficients" in page.data and r['hecke_orbit_label'] in page.data:
                        print 'Failed on', url
                        errors.append(r['hecke_orbit_label'])
        if errors:
            print "Errors occurred for the following labels: ", errors
        assert not errors
