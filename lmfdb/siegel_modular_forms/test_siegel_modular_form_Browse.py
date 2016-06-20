# -*- coding: utf8 -*-
from lmfdb.base import LmfdbTest
import math
import unittest2


class HomePageTest(LmfdbTest):

    def check(self,path,text):
        data = self.tc.get(path).data
        assert text in data

    # All tests should pass: these are all the links in the browse page 
    def test_Siegel_links_browse_page(self):
        r"""
        Check that the links work.
        """
        self.check("/ModularForm/GSp/Q/Sp4Z_j/10/0/", 'M_{10,0}')
        self.check("/ModularForm/GSp/Q/Sp4Z_j/",  'Upsilon')
        self.check("/ModularForm/GSp/Q/Kp/",  'in level 277, the')
        self.check("/ModularForm/GSp/Q/Sp6Z/",  'Miyawaki (1)')
        self.check("/ModularForm/GSp/Q/Sp8Z/",  'Other_II (2)')
        self.check("/ModularForm/GSp/Q/Gamma0_2/",  'Gamma_0(2)')
        self.check("/ModularForm/GSp/Q/Gamma1_2/",  'Gamma_1(2)')
        self.check("/ModularForm/GSp/Q/Gamma_2/",  'Gamma(2)')
        self.check("/ModularForm/GSp/Q/Gamma0_3/",  'Gamma_0(3)')
        self.check("/ModularForm/GSp/Q/Gamma0_3_psi_3/",  'T.Ibukiyama:')
        self.check("/ModularForm/GSp/Q/Gamma0_4/",  'Gamma_0(4)')
        self.check("/ModularForm/GSp/Q/Gamma0_4_psi_4/",  'psi_4')
        self.check("/ModularForm/GSp/Q/Gamma0_4_half/",  'k-1/2')


