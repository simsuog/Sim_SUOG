#
#  Copyright (c) 2023 - The SUOG Project.
#  @Author: Mirna El Ghosh
#

from owlready2 import *


class SUOGOnto(object):

    """this class implements the SUOG ontology"""

    def __init__(self, src):
        self._onto = owlready2.get_ontology(src).load()
        self._classes = [s for s in self._onto.classes()]
        self._rootiri = 'http://www.suog.org/ontology#Root'
        self._listclass = list(self._onto.classes())
        self._instances = list(self._onto.individuals())
        self._root = self._onto.search_one(iri=self._rootiri)

    def is_a(self, c):
        return (c.is_a)

    def getInstances(self, c):
        list = c.instances()
        return list

    def search_one(self, c):
        return self._onto.search_one(iri='*#'+c)

    def search_one_hpo(self, c):
        return self._onto.search_one(iri='*' + c)

    def search(self, c):
        return self._onto.search_one(iri=c)

    def search_label(self, c):
        return self._onto.search_one(prefLabel='*'+c)

    def ancestorsOf(self, c):
        return c.ancestors()

    def subClassesOf(self, c):
        return list(c.subclasses())

    def descendantOf(self, c):
        return list(c.descendants())

    def get_Label(self, irid):
        entity = self._onto.search_one(iri='*#'+irid)
        return entity.prefLabel.en[0]

    def get_Label_hpo(self, irid):
        entity = self._onto.search_one(iri='*'+irid)
        return entity.label[0]

    def get_Label2(self, irid):
        return irid.altLabel

    def getDisorders(self):
        disorder_top = self._onto.search(iri='*#FM0004*')
        return disorder_top[0].descendants()

    def getFindings(self):
        finding_top = self._onto.search(iri='*#FM0001*')
        return finding_top[0].descendants()

    def getPhenotypes(self, iri):
        finding_top = self._onto.search(iri='*'+iri)
        print(finding_top[0].descendants())
        return finding_top[0].descendants()

    def descendants(self, c):
        return c.descendants()

    def labelClassOf(self, c):
        return c.label

    def prefLabelClassOf(self, c):
        return c.prefLabel