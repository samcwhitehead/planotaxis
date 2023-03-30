# -*- coding: utf-8 -*-
"""
Simple example of loading UI template created with Qt Designer.

This example uses uic.loadUiType to parse and load the ui at runtime. It is also
possible to pre-compile the .ui file using pyuic (see VideoSpeedTest and 
ScatterPlotSpeedTest examples; these .ui files have been compiled with the
tools/rebuildUi.py script).
"""
# IMPORTS
import os
import h5py
import numpy as np
import pyqtgraph as pg

# muscle model
import muscle_model as mm

# import functions from extract_muscle_signals
from extract_muscle_signals import fit_to_model, extract_gcamp_signals

# import pyqt tools
from pyqtgraph.Qt import QtCore, QtGui#QStringList,QString

# ------------------------------------------------------------------------
import pyqtgraph.Qt
print 'pyqtgraph.Qt.USE_PYSIDE = ', pyqtgraph.Qt.USE_PYSIDE
print 'pyqtgraph.Qt.QtVersion  = ', pyqtgraph.Qt.QtVersion
# ------------------------------------------------------------------------

this_dir, this_filename = os.path.split(__file__)
model_path = os.path.join(this_dir, "models")

pg.mkQApp()

## Define main window class from template
path = os.path.dirname(os.path.abspath(__file__))
uiFile = os.path.join(path, 'viewer.ui')
WindowTemplate, TemplateBaseClass = pg.Qt.loadUiType(uiFile)

default_rframe_data = {'a1': np.array([ 51.5848967 ,  -5.93928407]),
                       'a2': np.array([ -0.09151179,  88.42505672]),
                       'p': np.array([ 26.66908747,  34.43488385])}


# --------------------------------------------------------------------------------------
class ModelView(object):

    def __init__(self,model):
        import copy
        self.model = model
        self.plot_frame = copy.copy(model.frame)
        self.curves = None
        self.element_list = []

    def plot(self,basis,plotobject):
        if self.curves:
            for pitem in self.curves:
                plotobject.removeItem(pitem)
        lines = self.model.coords_from_frame(basis)
        self.curves = list()
        for element_name, line in lines.items():
            if element_name in self.element_list:
                self.curves.append(plotobject.plot(line[0,:],line[1,:]))

    def update_basis(self,basis):
        lines = self.model.coords_from_frame(basis)
        lines = [l for k,l in lines.items() if k in self.element_list]
        if self.curves:
            for curve,line in zip(self.curves,lines):#lines.values()):
                curve.setData(line[0,:],line[1,:])

    def basis_changed(self,roi):
        pnts = roi.saveState()['points']
        p = np.array(pnts[1])

        a1 = np.array(pnts[0])-p
        a2 = np.array(pnts[2])-p

        self.plot_frame['p'] = p
        self.plot_frame['a1'] = a1
        self.plot_frame['a2'] = a2
        self.update_basis(self.plot_frame)


class RefrenceFrameROI(pg.ROI):
    
    def __init__(self, basis, closed=False, pos=None, **args):
        
        pos = [0,0]
        
        self.closed = closed
        self.segments = []
        pg.ROI.__init__(self, pos, **args)
        
        self.addFreeHandle((basis['p'][0]+basis['a1'][0],basis['p'][1]+basis['a1'][1]))
        self.addFreeHandle((basis['p'][0],basis['p'][1]))
        self.addFreeHandle((basis['p'][0]+basis['a2'][0],basis['p'][1]+basis['a2'][1]))

        for i in range(0, len(self.handles)-1):
            self.addSegment(self.handles[i]['item'], self.handles[i+1]['item'])
            
    def addSegment(self, h1, h2, index=None):
        seg = pg.LineSegmentROI(handles=(h1, h2), pen=self.pen, parent=self, movable=False)
        if index is None:
            self.segments.append(seg)
        else:
            self.segments.insert(index, seg)
        #seg.sigClicked.connect(self.segmentClicked)
        #seg.setAcceptedMouseButtons(QtCore.Qt.LeftButton)
        seg.setZValue(self.zValue()+1)
        for h in seg.handles:
            h['item'].setDeletable(False)
        
    def saveState(self):
        state = pg.ROI.saveState(self)
        state['closed'] = self.closed
        state['points'] = [tuple(h.pos()) for h in self.getHandles()]
        return state

    def setState(self,state):
        print state
        pg.ROI.setState(self,state,update = False)
        #state = pg.ROI.saveState(self)
        for h,p in zip(self.getHandles(),state['points']):
            self.movePoint(h,p)

        self.stateChanged(finish=True)
        return state

class MainWindow(TemplateBaseClass):  

    def __init__(self):
        TemplateBaseClass.__init__(self)
        self.setWindowTitle('muscle imaging browser')
        
        # Create the main window
        self.ui = WindowTemplate()
        #initialize the items created in designer
        self.ui.setupUi(self)
        
        #frame view
        self.plt = pg.PlotItem()
        self.ui.frameView.setCentralItem(self.plt)
        self.frameView = pg.ImageItem()
        self.plt.addItem(self.frameView)

        #gama plot
        self.gammaPlt = pg.PlotItem()
        self.ui.gammaPlot.setCentralItem(self.gammaPlt)
        self.ui.gammaSlider.valueChanged.connect(self.gammaChange)
        
        #default gama
        self.gammaf = lambda x: x**1
        self.gammax = np.linspace(0,2,100)
        self.gammaCurve = self.gammaPlt.plot(self.gammax,self.gammaf(self.gammax))

        #timeSeries
        self.timeSeriesPlt = pg.PlotItem()
        self.ui.timeSeriesPlt.setCentralItem(self.timeSeriesPlt)
        self.tserTrace = self.timeSeriesPlt.plot(np.ones(1000))
        self.tpointLine = pg.InfiniteLine(pos = 0,movable = True)
        self.tpointLine.sigPositionChanged.connect(self.tpointLineMoved)
        self.timeSeriesPlt.addItem(self.tpointLine)

        #load frames button
        self.ui.loadFrames.clicked.connect(self.loadFrames)

        #save data button
        self.ui.saveFit.clicked.connect(self.saveFit)
        self.ui.loadFit.clicked.connect(self.loadFit)

        ##scroll bar
        self.ui.frameScrollBar.valueChanged.connect(self.frameScrollBar_valueChanged)

        # Contrast/color control
        self.hist = pg.HistogramLUTItem()
        self.hist.setImageItem(self.frameView)
        self.ui.frameHist.setCentralItem(self.hist)

        self.componentsModel = QtGui.QStandardItemModel(self.ui.componentsView)
        self.ui.componentsView.setModel(self.componentsModel)
        self.componentsModel.itemChanged.connect(self.componentsChanged)

        #modelSelector
        self.loadedComponents = list()
        self.updateModelList()
        self.ui.modelselectBox.currentIndexChanged.connect(self.modelSelected)
        self.modelSelected(0)

        #profileSelector
        self.updateProfileList()
        self.ui.profileselectBox.currentIndexChanged.connect(self.profileSelected)
        self.profileSelected(0)
        self.ui.saveProfile.clicked.connect(self.saveProfile)

        #load outlines
        self.loadLines()
        self.current_frame = 0
        self.show()
        
        #self.ui.commentBox
        self.ui.frameNumber.setText(str(self.current_frame))
        self.ui.frameNumber.textEdited.connect(self.frameInput)

        #addEpoch
        self.epochPlots = dict()
        self.epoch_dict = dict()
        self.ui.newEpoch.clicked.connect(self.newEpoch)
        self.ui.saveEpoch.clicked.connect(self.saveEpoch)

        self.ui.epochStart.textEdited.connect(self.updateEpochFromText)
        self.ui.epochEnd.textEdited.connect(self.updateEpochFromText)
        
        #muscle demixing
        self.ui.applyDemixing.clicked.connect(self.extract_signals)  
        self.ui.subtractBackground.stateChanged.connect(self.subtractBackgroundChecked)
        self.subtract_background = self.ui.subtractBackground.isChecked()

    def subtractBackgroundChecked(self,i):
        self.subtract_background = self.ui.subtractBackground.isChecked()
        print self.subtract_background

    def profileSelected(self,i):
        import cPickle
        profile = self.ui.profileselectBox.currentText()
        with open(model_path + '/%s/profiles/%s'%(self.cur_model,profile),'rb') as f:
            profile_data = cPickle.load(f)
        for component in self.loadedComponents:
            if component['name'] in profile_data['selected_components']:
                component['checkbox'].setCheckState(True)
            else:
                component['checkbox'].setCheckState(False)
        self.ui.profileName.setText(profile)
        print profile_data


    def updateProfileList(self):
        import os
        profile_list = os.listdir(model_path + '/%s/profiles'%(self.cur_model))
        if len(profile_list) == 0:
            print 'creating default profile'
            import cPickle
            with open(model_path + '/%s/profiles/default.cpkl'%(self.cur_model),'wb') as f:
                cPickle.dump({'selected_components':[]},f)
            self.updateProfileList()
        else:
            for profile in profile_list:
                self.ui.profileselectBox.addItem(profile)
        index = self.ui.profileselectBox.findText('default.cpkl', QtCore.Qt.MatchFixedString)
        if index >= 0:
            self.ui.profileselectBox.setCurrentIndex(index)

    def saveProfile(self):
        profile_dir = model_path + '/%s/profiles/'%(self.cur_model)
        name = str(self.ui.profileName.text())
        print profile_dir + name
        f = open(profile_dir + name,'wb')
        import cPickle
        cPickle.dump({'selected_components':self.thorax_view.element_list},f)
        #self.updateProfileList()

    def componentsChanged(self):
        # If the changed item is not checked, don't bother checking others
        #if not item.checkState():
        #    return
        # Loop through the items until you get None, which
        # means you've passed the end of the list
        i = 0
        item_list = list()
        while self.componentsModel.item(i):
            if self.componentsModel.item(i).checkState():
                item_list.append(i)
                #return
            i += 1
        #skeys = self.signalshelf.keys()
        self.checked_signals = [self.loadedComponents[i]['name'] for i in item_list]
        print self.checked_signals
        self.thorax_view.element_list = self.checked_signals
        self.thorax_view.plot(self.thorax_view.plot_frame,self.plt)
        self.roi.stateChanged()

        #self.thorax_view.plot()
        #self.update_tser_plot()

    def modelSelected(self,i):
        import cPickle
        self.cur_model = os.listdir(model_path)[i]
        #print self.cur_modelcomponentsModel
        with open(model_path + '/%s/outlines.cpkl'%(self.cur_model),'rb') as f:
            self.outlines = cPickle.load(f)
        for key in sorted(self.outlines.keys()):
            print key
            item = QtGui.QStandardItem(key)
            #check = 1 if np.random.randint(0, 1) == 1 else 0
            item.setCheckable(True)
            item.setCheckState(False)
            self.loadedComponents.append({'checkbox':item,'name':key})
            self.componentsModel.appendRow(item)
            #self.color_dict[key] = 'r'
    
    def roiClicked(self,item):
        #print item.mname
        #print 'here'
        color = pg.QtGui.QColorDialog.getColor()
        #self.color_dict[item.mname] = color
        #self.update_tser_plot()

    def updateModelList(self):
        import os
        for mstr in os.listdir(model_path):
            self.ui.modelselectBox.addItem(mstr)

    def newEpoch(self):
        name = str(self.ui.epochName.text())
        print name
        if (not(name in self.epoch_dict.keys()) and not(name == '')):
            epoch_range = [self.current_frame,self.current_frame + 100]
            self.epoch_dict[name] = epoch_range
            self.plotEpoch(name)
            ep_plot = self.epochPlots[name]
            sta,stp = ep_plot.getRegion()
            self.ui.epochStart.setText(str(int(sta)))
            self.ui.epochEnd.setText(str(int(stp)))

    def clearEpochs(self):
        for k in self.epoch_dict.keys():
            self.timeSeriesPlt.removeItem(self.epochPlots[k])
            self.epochPlots.pop(k)
            self.epoch_dict.pop(k)

    def plotEpoch(self,k):
        ep = pg.LinearRegionItem(values= self.epoch_dict[k])
        ep.epoch_name = k
        ep.sigRegionChanged.connect(self.updateEpochPlot)
        self.epochPlots[k] = ep
        self.timeSeriesPlt.addItem(ep)
        self.tpointLine.setZValue(ep.zValue()+1)

    def updateEpochPlot(self,ep):
        self.ui.epochName.setText(ep.epoch_name)
        self.updateCurrentEpochState()

    def updateEpochFromText(self):#negative indexing
        k = str(self.ui.epochName.text())
        ep_plot = self.epochPlots[k]
        sta = int(self.ui.epochStart.text())
        stp = int(self.ui.epochEnd.text())
        ep_plot.setRegion((sta,stp))
        self.epoch_dict[k] = [sta,stp]

    def updateCurrentEpochState(self):
        k = str(self.ui.epochName.text())
        ep = self.epoch_dict[k]
        ep_plot = self.epochPlots[k]
        sta,stp = ep_plot.getRegion()
        self.ui.epochStart.setText(str(int(sta)))
        self.ui.epochEnd.setText(str(int(stp)))
        self.epoch_dict[k] = [int(sta),int(stp)]

    def saveEpoch(self):
        flydir = '%s%s/'%(dba.root_dir,self.current_fly)
        f = open(flydir + 'epoch_data.cpkl','wb')
        import cPickle
        cPickle.dump(self.epoch_dict,f)
        print self.epoch_dict

    def frameInput(self,value):
        self.current_frame = int(value)
        self.showFrame()

    def tpointLineMoved(self):
        self.current_frame = int(self.tpointLine.value())
        self.showFrame()

    def gammaChange(self,value):
        gamma = value/50.0
        self.gammaf = lambda x: x**gamma
        #print gamma
        self.gammaCurve.setData(self.gammax,self.gammaf(self.gammax))
        self.showFrame()

    def loadfileTree(self):
        self.ui.fileTree.setColumnCount(1)
        items = []
        #for key,fly in zip(fly_db.keys(),fly_db.values()):
        for key,fly in sorted(fly_db.items()):#zip(fly_db.keys(),fly_db.values()):
            #print key
            try:
                exp1 = fly['experiments'].values()[0]
                exptype = fly['experiments'].keys()[0]
                if 'tiff_data' in exp1.keys():
                    #item_list.append('fly%s'%key)
                    item = QtGui.QTreeWidgetItem(None,['Fly%04d'%int(key)])
                    for img_key in ['images','refstack']:
                        if img_key in exp1['tiff_data'].keys():
                            #data_ref = exp1['tiff_data'][img_key]
                            child = QtGui.QTreeWidgetItem(None,[img_key])
                            child.setData(0,QtCore.Qt.UserRole,key)
                            item.insertChild(0,child)
                            items.append(item)
                            #print (img_key,np.shape(exp1['tiff_data'][img_key]))
                        else:
                            pass
                else:
                    print exp1.keys()
            except KeyError:
                pass
        self.ui.fileTree.insertTopLevelItems(0,items)

    def loadLines(self):
        import cPickle
        #f = open('model_data.cpkl','rb')
        ###f = open('/media/flyranch/ICRA_2015/model_data.cpkl','rb')
        model_data = self.outlines
        #f.close()

        ########################
        #model_keys = []
        e1 = model_data['e1']
        e2 = model_data['e2']

        muscle_dict = dict()
        for key in model_data.keys():
            if not(key in ['e1','e2']):
                muscle_dict[key] = model_data[key]
        frame = mm.Frame()
        frame['a2'] = e1[1]-e2[0]
        frame['a1'] = e2[1]-e2[0]
        frame['p'] = e2[0]
        thorax = mm.GeometricModel(muscle_dict,frame)
        self.thorax_view = ModelView(thorax)
        self.roi = RefrenceFrameROI(thorax.frame)
        self.roi.sigRegionChanged.connect(self.thorax_view.basis_changed)
        #self.roi.sigRegionChanged.connect(self.affineWarp)

        self.plt.disableAutoRange('xy')
        
        state = self.roi.getState()
        rf = default_rframe_data
        pnts = [(rf['p'][0]+rf['a1'][0],rf['p'][1]+rf['a1'][1]),
                 (rf['p'][0],rf['p'][1]),
                 (rf['p'][0]+rf['a2'][0],rf['p'][1]+rf['a2'][1])]
        state['points'] = pnts
        self.roi.setState(state)
        self.roi.stateChanged()
        self.plt.addItem(self.roi)

        self.thorax_view.plot(self.thorax_view.plot_frame,self.plt)

    def loadFrames(self):
        self.CurrentHDF5FileName = str(QtGui.QFileDialog.getOpenFileName(self, 
                                        'Dialog Title',
                                        '')[0])

        self.image_prefix = self.CurrentHDF5FileName.split('/')[-1].split('.hdf5')[0]
        print self.image_prefix
        #tfile = tifffile.TiffFile(self.CurrentTiffFileName)
        #self.images = np.array(fly_db[fnum]['experiments'].values()[0]['tiff_data']['images'])
        import os
        self.CurrentFlyPath = os.path.split(self.CurrentHDF5FileName)[0]
        print self.CurrentFlyPath
        self.FlyFile = h5py.File(self.CurrentHDF5FileName,'r')
        key = self.FlyFile.keys()[0]

        #tfile = tifffile.TiffFile(self.CurrentTiffFileName)
        self.images = np.array(self.FlyFile[key])
        
        try:
            f = open('basis_fits.cpkl','rb')
            import cPickle
            basis = cPickle.load(f)
            state = self.roi.getState()
            pnts = [(basis['p'][0]+basis['a1'][0],basis['p'][1]+basis['a1'][1]),
                    (basis['p'][0],basis['p'][1]),
                    (basis['p'][0]+basis['a2'][0],basis['p'][1]+basis['a2'][1])]
            state['points'] = pnts
            self.roi.setState(state)
            self.roi.stateChanged()
            self.ui.commentBox.setPlainText(basis['commentBox'])
        except IOError:
            print 'no file'
            self.ui.commentBox.setPlainText('')

        self.clearEpochs()

        try:
            f = open('epoch_data.cpkl','rb')
            import cPickle
            self.epoch_dict = cPickle.load(f)
            for k in self.epoch_dict.keys():
                self.plotEpoch(k)
            self.ui.epochName.setText(self.epoch_dict.keys()[0])
            self.updateCurrentEpochState()
        except IOError:
            print 'no epoch file'
            self.ui.epochName.setText('')
            self.ui.epochStart.setText('')
            self.ui.epochEnd.setText('')

        #self.frameView.setImage(self.images[0,:,:])
        self.current_frame = 0
        self.showFrame()
        #self.transformImage.setImage(self.transform_img.astype(np.float32))
        self.ui.frameScrollBar.setMaximum(np.shape(self.images)[0])
        self.plt.autoRange()
        #set transformImage
        

    def showFrame(self):
        try:
            if self.current_frame > 0:
                img = self.gammaf(np.array(self.images[self.current_frame,:,:]).astype(np.float32))
                # typically need to do vertical flip (I think this has to do with the way ROS saves images)
                img = np.flipud(img)
                # check if we need to rotated image 90 degrees
                if img.shape[1] > img.shape[0]:
                    img = np.transpose(img)

                self.frameView.setImage(img.astype(np.float32))
                self.ui.frameNumber.setText(str(self.current_frame))
                self.ui.frameScrollBar.setValue(self.current_frame)
                self.tpointLine.setValue(self.current_frame)
            else:
                raise ValueError('Negative indexing not allowed')
        except ValueError as er:
            print er


    def affineWarp(self,roi):
        src_f = self.thorax_view.plot_frame
        dst_f = self.thorax_view.model.basis

        dst_p0 = dst_f['a1'] + dst_f['p']
        dst_p1 = dst_f['p']
        dst_p2 = dst_f['a2'] + dst_f['p']

        src_p0 = src_f['a1'] + src_f['p']
        src_p1 = src_f['p']
        src_p2 = src_f['a2'] + src_f['p']

    def frameScrollBar_valueChanged(self,value):
        #self.frameView.setImage(self.images[value,:,:])
        self.current_frame = value
        self.showFrame()
        
    def saveFit(self):
        import cPickle
        savedata = dict(self.thorax_view.plot_frame)
        comment_text = self.ui.commentBox.toPlainText()
        savedata['commentBox'] = comment_text
        
        with open(os.path.join(self.CurrentFlyPath,self.image_prefix + '_rframe_fits.cpkl'),'wb') as f:
            cPickle.dump(savedata,f)

    def loadFit(self):
        import cPickle
        with open(os.path.join(self.CurrentFlyPath, self.image_prefix + '_rframe_fits.cpkl'),'rb') as f:
            loaddata = cPickle.load(f)
        print loaddata        
        state = self.roi.getState()
        rf = loaddata
        pnts = [(rf['p'][0]+rf['a1'][0],rf['p'][1]+rf['a1'][1]),
                 (rf['p'][0],rf['p'][1]),
                 (rf['p'][0]+rf['a2'][0],rf['p'][1]+rf['a2'][1])]
        state['points'] = pnts
        self.roi.setState(state)
        self.roi.stateChanged()
        #self.thorax_view.update_basis(loaddata)
        #print self.ui.fileTree.selectedItems()[0].data(0,QtCore.Qt.UserRole).toPyObject()

    # ------------------------------------------------------------------------------
    # SCW update to use an imported function, thereby removing some redundant code
    def extract_signals(self):
        # read data fields from self
        imgs = np.array(self.images)
        fly_frame = self.thorax_view.plot_frame
        muscles = self.thorax_view.element_list

        # some hard-coded values that should really be made into variables
        fit_mode = 'pinv'
        chunk_sz = 100

        # run signal extraction
        signals_dict = extract_gcamp_signals(imgs, fly_frame, driver=None, model_type='volumetric',
                                             model_name='thorax', model_dir=model_path, fit_mode=fit_mode,
                                             fit_pix_mask=None, baseline=None, chunk_sz=chunk_sz,
                                             muscles=muscles)

        # save extracted signals
        import cPickle
        fname = os.path.join(self.CurrentFlyPath, self.image_prefix + '_model_fits.cpkl')
        with open(fname, 'wb') as f:
            cPickle.dump(signals_dict, f)

        # print a little update
        print muscles
        print 'Done!'

win = MainWindow()

## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
    #fly_db.close()
