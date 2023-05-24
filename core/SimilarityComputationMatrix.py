#
#  Copyright (c) 2023 - The SUOG Project.
#  @Author: Mirna El Ghosh
#

#
#
import itertools
import numpy as np
import csv
import time
from datetime import date
import pandas as pd

from core.OntologyGraph import SUOGOntologyGraph, Taxonomy
from core.SUOGOntology import SUOGOnto


class SimilarityComputation(object):

    def __init__(self, src, file, measure):
        self._onto = SUOGOnto(src)
        self.og = SUOGOntologyGraph(src)
        self._tax = Taxonomy(self.og)
        self._length = len(src[src.find('/'):src.find('.owl')])
        self._version = src[src.find('v')+1:src.find('v')+5]
        self._ontologyfile = src[src.find('/')+1:src.find('.owl')]
        self.imagesfile = 'csv/input/'+self._version+'/images_similarity.csv'
        self.input = file
        self.output = 'csv/output/'+self._version+'/'
        self.measure = measure
        self._alpha = 0.4
        self._beta = 0.3
        self._gamma = 0.3
        self._alphap = 0.1
        self._betap = 0.3
        self._gammap = 0.6
        self.afindings = []
        self.pfindings = []
        self.disorders = []
        self.routes = []
        self.modes = []
        self.views = []

    def select_annotations(self):
        images=[]
        with open(self.input, newline='') as f:
            reader = csv.reader(f)
            for r in reader:
                images.append(r[0])
        print("list---------------------", images)
        #self.outputfile = 'images_similarity.csv'

        with open(self.imagesfile, 'w', newline='', encoding="utf-8") as file:
            writer = csv.writer(file, delimiter=',')
            for image in images:
                i = self._onto.search_one(image)
                label = i.FMr0020

                disorders = []
                findings = []
                views = []
                modes = []
                routes = []

                findings_i = list(i.FMr0022)
                disorders_i = list(i.FMr0021)
                print('findings: ', findings_i)
                modes_i = list(i.FMr0026)
                routes_i = list(i.FMr0027)
                views_i = list(i.FMr0025)

                if findings_i:
                    for f in findings_i:
                        print(self._onto.is_a(f))
                        findingsl = self._onto.is_a(f)
                        for fl in findingsl:
                            findings.append(fl)
                        print(findings)
                if disorders_i:
                    for d in disorders_i:
                        print('disorder', d)
                        disordersl = d.is_a
                        for dl in disordersl:
                            disorders.append(dl)
                if modes_i:
                    for m in modes_i:
                        modesl = m.is_a
                        for ml in modesl:
                            modes.append(ml)
                if routes_i:
                    for r in routes_i:
                        routesl = r.is_a
                        for rl in routesl:
                            routes.append(rl)
                if views_i:
                    for v in views_i:
                        viewsl = v.is_a
                        for vl in viewsl:
                            views.append(vl)
                print('annotations', label, i, findings, disorders, modes, routes, views)
                annotations = [label, i, findings, disorders, modes, routes, views]
                #writer.writerow((label, i, findings, disorders, modes, routes, views))
                #annotations_ = (','.join(annotations))
                writer.writerow(annotations)
           # for a in annotations:
                #writer.writerow(a)

    def get_annotations(self):
        listimages = []
        imagesdict2 = {}
        imagesid = []
        with open(self.imagesfile, mode='r') as infile:
            reader = csv.reader(infile)
            for r in reader:
                print(r)
                listimages.append(r)

        for i in listimages:
            #print('----->>>>', i[1])
            image=i[1]
            imageid = str(image)[self._length:]
            #print('image id:', imageid)
            imagesid.append(imageid)
            imagesdict2[imageid] = {}
            imagesdict2[imageid]['test'] = (i[0]).split(',')
            imagesdict2[imageid]['findings'] = [x.strip().lstrip('[').strip(']') for x in i[2].split(',')]
            imagesdict2[imageid]['disorders'] = [x.strip().lstrip('[').strip(']') for x in i[3].split(',')]
            imagesdict2[imageid]['modes'] = [x.strip().lstrip('[').strip(']') for x in i[4].split(',')]
            imagesdict2[imageid]['routes'] = [x.strip().lstrip('[').strip(']') for x in i[5].split(',')]
            imagesdict2[imageid]['views'] = [x.strip().lstrip('[').strip(']') for x in i[6].split(',')]
           # for i in imagesdict2.items():
            #    print((i[1]['findings']))
        return imagesdict2

    def splitAnnotations(self):
        listimages = self.get_annotations()
        for i in listimages.items():
            print(i)

            dlist = i[1]['disorders']
            if dlist != ['']:
                for d in dlist:
                    if d[self._length:] not in self.disorders:
                        self.disorders.append(d[self._length:])
                        print('add disorder: ', d)
            mlist = i[1]['modes']
            if mlist != ['']:
                for m in mlist:
                    if m[self._length:] not in self.modes:
                        self.modes.append(m[self._length:])
                        print('add mode: ', m)

            rlist = i[1]['routes']
            if rlist != ['']:
                for r in rlist:
                    if r[self._length:] not in self.routes:
                        self.routes.append(r[self._length:])
                        print('add route: ', r)

            vlist = i[1]['views']
            if vlist != ['']:
                for v in vlist:
                    if v[self._length:] not in self.views:
                        self.views.append(v[self._length:])
                        print('add view: ', v)

            flist = i[1]['findings']
            if flist != ['']:
                for f in flist:
                    if self.is_anatomical(self.in_suog(f)) and f[self._length:] not in self.afindings:
                        self.afindings.append(f[self._length:])
                        print('add anatomical findings ', f)

                    elif f[self._length:] not in self.pfindings:
                         self.pfindings.append(f[self._length:])
                         print('add pathological finding: ', f)

    def computeSimilarityAllAnnotations(self):
        self.splitAnnotations()
        self.computeSimilarityAnnotations('pfindings', self.pfindings)
        self.computeSimilarityAnnotations('disorders', self.disorders)
        self.computeSimilarityAnnotations('afindings', self.afindings)
        self.computeSimilarityAnnotations('views', self.views)
        self.computeSimilarityAnnotations('routes', self.routes)
        self.computeSimilarityAnnotations('modes', self.modes)

    def computeSimilarityAnnotations(self, category, list):
        m = self.measure
        if m == "sim_suog":
            if category == 'pfindings' or category == 'disorders':
                measure = self._tax.sim_suog
            else:
                measure = self._tax.sim_suog_t
        elif m == "sim_ic":
            measure = self._tax.sim_ic
        elif m == "sim_resnik":
            measure = self._tax.sim_resnik

        with open('csv/output/'+self._version+'/'+category+'.csv','w', newline='', encoding="utf-8") as categoryfile:
            writer = csv.writer(categoryfile, delimiter=',')
            list2 = ['x']
            writer.writerow(list2 + list)
            counter1 = 0
            counter2 = 0
            for a1 in list:
                counter1 += 1
                sima = []
                sima.append(a1)
                for a2 in list:
                    counter2 += 1
                    print(category,': ', a1, ':', counter1, '-', a2, ':', counter2)
                    sima.append(measure(self.in_suog_matrix(a1),self.in_suog_matrix(a2)))
                    print(sima)
                writer.writerow(sima)

    def addAnnotations(self, list1, list2, annot1, annot2, label1, label2):
        if list1 == '':
            label1 += ""
        else:
            for f1 in list1:
                f11 = f1[self._length:]
                label1 += "-" + str(self._onto.get_Label(f11))

        if list2 == '':
            label2 += ""
        else:
            for f2 in list2:
                f22 = f2[self._length:]
                label2 += "-" + self._onto.get_Label(f22)
        annot1.append(label1)
        annot2.append(label2)

    def in_suog_matrix(self, c):
        c1 = self._onto.search_one(c)
        return c1

    def in_hpo(self, c):
        c1 = self._onto.search_one_hpo(c)
        return c1

    def in_suog(self, c):
        c1 = c[self._length:]
        c11 = self._onto.search_one(c1)
        return c11

    def locateSimilarity_list(self, list1, list2, result, category):
        path = 'csv/output/'+self._version+'/' + category + '.csv'
        df = pd.read_csv(path)
        df = df.set_index('x')
        res = set(itertools.product(list1, list2))
        for (e1, e2) in res:
            #print(e1, e2)
            data = df[e1[self._length:]]
            sim = data[e2[self._length:]]
            result.append(sim)
        return result

    def normalAppearance_exist(self, findings1):
       # normal = self._onto.search_one("UF5714")
       # singleton = self._onto.search_one("UF2727")
        normal = 'UF5714'
        singleton = 'UF2727'
        f1='[%s]' % ', '.join(map(str, findings1))
        #print('---------findings: ', f1 , str(normal))
        #print(f1)

        return str(normal) in f1 or str(singleton) in f1

    def similarity_entities_Matrix(self, image1, image2, annot1, annot2, dlabel1, dlabel2, label, index):
        print('----------------------------', label, '-----------------------')
        print(image1[0], image1[index][label], image2[0], image2[index][label])
        e1 = image1[index][label]
        e2 = image2[index][label]
        sim = []

        if e1 == [''] or e2 == ['']:
            simentities = 0
        else:
            self.addAnnotations(e1, e2, annot1, annot2, dlabel1, dlabel2)
            if e1 == e2:
                simentities = 1
                print('similar ', label, ' : ', simentities)
            else:
                if label == 'findings':
                    if self.normalAppearance_exist(e1) and self.normalAppearance_exist(e2):
                            e11 = [e for e in e1 if self.is_in_anatomical(e[self._length:])]
                            e22 = [e for e in e2 if self.is_in_anatomical(e[self._length:])]
                            self.locateSimilarity_list(e11, e22, sim, "afindings")
                    else:
                        if not self.normalAppearance_exist(e1) and not self.normalAppearance_exist(e2) :
                            #e11 = [e for e in e1 if not self.is_anatomical(self.in_suog(e))]
                            #e22 = [e for e in e2 if not self.is_anatomical(self.in_suog(e))]
                            e11 = [e for e in e1 if not self.is_in_anatomical(e[self._length:])]
                            e22 = [e for e in e2 if not self.is_in_anatomical(e[self._length:])]
                            self.locateSimilarity_list(e11, e22, sim, "pfindings")

                    if sim == []:
                        simentities = 0
                    else:
                        simentities = np.average(sim)
                        print('mean sim findings: ', simentities)

                ####################disorders-technical elements#######################
                else:
                    sim = []
                    self.locateSimilarity_list(e1, e2, sim, label)
                    simentities = np.average(sim)
                    print('mean sim - ', label, '- : ', simentities)

        return simentities

    def is_anatomical(self, e):
        anatomicalfindings = self._onto.search_one("UF1220")
        super = self._onto.ancestorsOf(e)
        return anatomicalfindings in super

    def is_in_anatomical(self, e):
        path = 'csv/output/'+self._version+'/afindings.csv'
        df = pd.read_csv(path)
        anatomicalfindings = df.iloc[0]
        return e in anatomicalfindings

    def similarity_phenotypes(self):

        print(self._ontologyfile)
        phenotypes = []
        with open(self.input, 'r') as f:
            reader = csv.reader(f, delimiter=',')
            for row in reader:
                phenotypes.append(row[0])

        m = self.measure
        if m == "sim_suog":
            measure = self._tax.sim_suog
        elif m == "sim_ic":
            measure = self._tax.sim_ic
        elif m == "sim_resnik":
            measure = self._tax.sim_resnik

        with open('csv/phenotypes/results_'+self._ontologyfile+'.csv', 'w', newline='',
                  encoding="utf-8") as phenotypefile:
            writer = csv.writer(phenotypefile, delimiter=',')
            list2 = ['x']
            writer.writerow(list2 + phenotypes)
            counter1 = 0
            counter2 = 0
            for a1 in phenotypes:
                counter1 += 1
                sima = [a1]
                for a2 in phenotypes:
                    counter2 += 1
                    print('phenotypes: ', a1, ':', counter1, '-', a2, ':', counter2)
                    if self._ontologyfile == 'hp':
                        sima.append(measure(self.in_hpo(a1), self.in_hpo(a2)))

                    else:
                        sima.append(measure(self.in_suog_matrix(a1), self.in_suog_matrix(a2)))
                    print(sima)
                writer.writerow(sima)

    def similarity_images(self):

        self.select_annotations()
        self.computeSimilarityAllAnnotations()

        previous = []
        counter1 = 0
        today = date.today()
        d = str(today)
        with open(self.output+'/similarity_images_' + self.measure + '_' + d + '.csv', 'w',
                  newline='', encoding="utf-8") as file:
            start = time.time()
            writer = csv.writer(file, delimiter=',')
            writer.writerow(["Ultrasound Study 1", "Image ID 1", "Annotations 1", "Ultrasound Study 2", "Image ID 2",
                             "Annotations 2", "SimFindings", "SimDisorders", "SimRoutes", "SimModes",
                             "SimViews", "Aggregation","Similarity MEAN", "Ultrasound Image1", "Ultrasound Image 2"])
            annotations = self.get_annotations()
            for image1 in annotations.items():
                imagesim = []
                imageid1 = image1[0]
                test1 = image1[1]['test']
                counter1 += 1
                counter2 = 0
                for image2 in annotations.items():
                    imageid2 = image2[0]
                    test2 = image2[1]['test']
                    if {imageid1, imageid2} not in previous and imageid1 != imageid2:
                        counter2 += 1
                        print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Similarity computation of image1:', counter1,
                              '- image2:', counter2)
                        annot1 = []
                        annot2 = []
                        flabel1 = ""
                        dlabel1 = ""
                        rlabel1 = ""
                        vlabel1 = ""
                        mlabel1 = ""

                        flabel2 = ""
                        dlabel2 = ""
                        rlabel2 = ""
                        vlabel2 = ""
                        mlabel2 = ""

                        simfindings = self.similarity_entities_Matrix(image1, image2, annot1, annot2, flabel1, flabel2,
                                                                      'findings', 1)
                        simdisorders = self.similarity_entities_Matrix(image1, image2, annot1, annot2, dlabel1, dlabel2,
                                                                       'disorders', 1)
                        simmodes = self.similarity_entities_Matrix(image1, image2, annot1, annot2, mlabel1, mlabel2,
                                                                   'modes', 1)
                        simroutes = self.similarity_entities_Matrix(image1, image2, annot1, annot2, rlabel1, rlabel2,
                                                                    'routes', 1)
                        simviews = self.similarity_entities_Matrix(image1, image2, annot1, annot2, vlabel1, vlabel2,
                                                                   'views', 1)
                        simtech = (self._alphap * simroutes + self._betap * simmodes + self._gammap * simviews)
                        agg = np.round(
                                np.abs((self._alpha * simfindings + self._beta * simdisorders + self._gamma * simtech)), 2)
                        print('aggregation: ', agg)
                        similarity = np.abs(1 - np.abs(np.log(agg)))

                        print('----Similarity result----:', test1, imageid1, test2, imageid2, ' - Sim Findings: ',
                              simfindings, ' - Sim Disorders: ',
                              simdisorders, ' - Sim Technical Elements: ', simtech, ' - Similarity: ', similarity)
                        imagesim.append((test1, imageid1, annot1, test2, imageid2, annot2, simfindings, simdisorders,
                                         simroutes, simmodes, simviews, agg, similarity))
                        print('similarity score: ', np.round(similarity, 2))
                        previous.append({imageid1, imageid2})

                        elapsed_time = time.time() - start
                        print('Execution time:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

                timetotal = time.time() - start
                print('Execution time:', time.strftime("%H:%M:%S", time.gmtime(timetotal)))
                for r in imagesim:
                    writer.writerow(r)