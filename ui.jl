using PyCall
include("src/Unveilapp.jl")
@pyimport PyQt6
@pyimport sys
qwid = pyimport("PyQt6.QtWidgets")
qpro = pyimport("PyQt6.QtCore")



@pydef mutable struct MainWindow <: qwid.QMainWindow #class MainWindow(QMainWindow):
    function __init__(self)  #def __init__(self):
        py"""super($MainWindow,$self).__init__()"""
        #qwid.QMainWindow.__init__(self)
        #super(MainWindow,self).__init__()
        self.ui()


    end

    function ui(self)
        backcolid = "#5D5D5D"
        bordsize = 2
        stylesheet = "QGroupBox { border: $(bordsize)px solid black; background: $backcolid; font: bold 20px ; margin-top: 4.5ex;subcontrol-origin:  margin } QGroupBox::title{subcontrol-origin: margin ; subcontrol-position : top middle; padding: 25px 0px -10px 0px;  }"
        self.setWindowTitle("Unveil")
        #self.setStyleSheet("background-color: #e6e6e6;")
        # GLOBAL LAYOUT, CONTAINS EVERYTHING
        layout = qwid.QGridLayout()
        self.setLayout(layout)

        #GROUP ESSENTIAL PARAMETER
        kinggroup = qwid.QGroupBox("Essential parameters")
        form_layout = qwid.QGridLayout()
        kinggroup.setObjectName("ColoredGroupBox")  
        kinggroup.setStyleSheet("QGroupBox { border: $(bordsize)px solid black; background: $backcolid; font: bold 20px ; margin-top: 3ex;subcontrol-origin:  margin } QGroupBox::title{subcontrol-origin: margin ; subcontrol-position : top middle; padding: 10px 0px -10px 0px;  }")  
        # self.fitsourcename = "FITSNAME"
        # self.fitsourcepath = qwid.QLineEdit("FITSOURCE PATH")
        # self.fitsourcebtn = qwid.QPushButton("Browse")
        # form_layout.addWidget(qwid.QLabel("Enter the path to fits source : "),0,0)
        # form_layout.addWidget(self.fitsourcepath,0,1)
        # form_layout.addWidget(self.fitsourcebtn,0,3) 
        # self.namefitsource="FITSNAME.fits"
        # py"""
        # $self.fitsourcebtn.clicked.connect(lambda: $self.open_file($self.fitsourcepath,$self.fitsourcename))
        # """

        self.savepath = qwid.QLineEdit("PATH TO SAVE DATA")
        self.savepathbtn = qwid.QPushButton("Browse")
        form_layout.addWidget(qwid.QLabel("Enter the path to saving directory : "),1,0)
        form_layout.addWidget(self.savepath,1,1)
        form_layout.addWidget(self.savepathbtn,1,3) 
        py"""
        $self.savepathbtn.clicked.connect(lambda: $self.open_dir($self.savepath))
        """

        self.savename = qwid.QLineEdit("NAME FOR SAVED DATA")
        form_layout.addWidget(qwid.QLabel("Enter the generic name to give to saved data : "),2,0)
        form_layout.addWidget(self.savename,2,1)


        self.noisebtn1 = qwid.QSpinBox()
        self.noisebtn1.setMinimum(1)
        self.noisebtn1.setMaximum(2e9)
        form_layout.addWidget(qwid.QLabel("Indices of first channel with noise "),3,0)
        form_layout.addWidget(self.noisebtn1,3,1)
        self.noisebtn2 = qwid.QSpinBox()
        self.noisebtn2.setMaximum(2e9)
        form_layout.addWidget(qwid.QLabel("Indices of first channel with noise "),4,0)
        form_layout.addWidget(self.noisebtn2,4,1)

        self.blankvalbtn = qwid.QDoubleSpinBox()
        py"""$self.blankvalbtn.setMinimum(float("-inf"))
        $self.blankvalbtn.setMaximum(float("inf"))    """
        form_layout.addWidget(qwid.QLabel("Blank value of your data "),5,0)
        form_layout.addWidget(self.blankvalbtn,5,1)


        self.units = qwid.QLineEdit("Units")
        form_layout.addWidget(qwid.QLabel("Units of the velocity axis (see header, should be m/s or km/s) "),6,0)
        form_layout.addWidget(self.units,6,1)
        kinggroup.setLayout(form_layout)





        # PCA GROUP
        pcabox = qwid.QGroupBox("PCA")
        formpca = qwid.QGridLayout()
        pcabox.setLayout(formpca)
        npc = qwid.QSpinBox()
        npc.setMinimum(1)
        npc.setMaximum(150)
        MainWindow.button_is_checked=true
        self.pcabtn = qwid.QPushButton("Go PCA")
        self.pcabtn.setCheckable(false)
        self.pcabtn.setEnabled(true)
        py"""
        $self.pcabtn.clicked.connect(lambda: $self.callpca($npc.value()))
        """
        self.overcheckpca = qwid.QCheckBox()
        self.fitsourcenamepca = "FITSNAME"
        self.fitsourcepathpca = qwid.QLineEdit("FITSOURCE FOR PCA")
        self.fitsourcebtnpca = qwid.QPushButton("Browse")
        formpca.addWidget(qwid.QLabel("Enter the path to fits source : "),0,0)
        formpca.addWidget(self.fitsourcepathpca,0,1)
        formpca.addWidget(self.fitsourcebtnpca,0,3) 
        self.namefitsourcepca="FITSNAME.fits"
        py"""
        $self.fitsourcebtnpca.clicked.connect(lambda: $self.open_file($self.fitsourcepathpca,$self.fitsourcenamepca))
        """
        formpca.addWidget(qwid.QLabel("Number of PC:"),1,0)
        formpca.addWidget(npc,1,1)
        formpca.addWidget(qwid.QLabel("Overwrite PCA saved cube?"),2,0)
        formpca.addWidget(self.overcheckpca,2,1)
        formpca.addWidget(self.pcabtn,3,0,1,2)
        pcabox.setObjectName("ColoredGroupBox")  
        pcabox.setStyleSheet("$stylesheet") 


        # SWO GROUP
        swobox = qwid.QGroupBox("SWO")
        formswo = qwid.QGridLayout()
        swobox.setLayout(formswo)
        npc = qwid.QSpinBox()
        npc.setMinimum(1)
        npc.setMaximum(2e9)
        MainWindow.button_is_checked=true
        self.swobtn = qwid.QPushButton("Go SWO")
        self.swobtn.setCheckable(false)
        self.swobtn.setEnabled(true)
        py"""
        $self.swobtn.clicked.connect(lambda: $self.callswo($npc.value()))
        """
        self.fitsourcenameswo = "FITSNAME"
        self.fitsourcepathswo = qwid.QLineEdit("FITSOURCE FOR SWO")
        self.fitsourcebtnswo = qwid.QPushButton("Browse")
        formswo.addWidget(qwid.QLabel("Enter the path to fits source : "),0,0)
        formswo.addWidget(self.fitsourcepathswo,0,1)
        formswo.addWidget(self.fitsourcebtnswo,0,3) 
        self.namefitsourcepca="FITSNAME.fits"
        py"""
        $self.fitsourcebtnswo.clicked.connect(lambda: $self.open_file($self.fitsourcepathswo,$self.fitsourcenameswo))
        """
        self.overcheckswo = qwid.QCheckBox()
        formswo.addWidget(qwid.QLabel("Overwrite SWO saved cube?"),1,0)
        formswo.addWidget(self.overcheckswo,1,1)
        formswo.addWidget(self.swobtn,2,0,1,2)
        swobox.setObjectName("ColoredGroupBox")  
        swobox.setStyleSheet("$stylesheet") 


        # CONVPCA
        convpcabox = qwid.QGroupBox("Convergence of the PCA")
        formconvpca = qwid.QGridLayout()
        convpcabox.setLayout(formconvpca)
        self.npchigh = qwid.QSpinBox()
        self.npchigh.setMinimum(1)
        self.npchigh.setMaximum(2e9)
        MainWindow.button_is_checked=true
        self.convpcabtn = qwid.QPushButton("Go ConvPCA")
        self.convpcabtn.setCheckable(false)
        self.convpcabtn.setEnabled(true)
        py"""
        $self.convpcabtn.clicked.connect(lambda: $self.callconvpca())
        """
        self.fitsourcenamecpca = "FITSNAME"
        self.fitsourcepathcpca = qwid.QLineEdit("FITSOURCE FOR CONV PCA")
        self.fitsourcebtncpca = qwid.QPushButton("Browse")
        formconvpca.addWidget(qwid.QLabel("Enter the path to fits source : "),0,0)
        formconvpca.addWidget(self.fitsourcepathcpca,0,1)
        formconvpca.addWidget(self.fitsourcebtncpca,0,3) 
        self.namefitsourcecpca="FITSNAME.fits"
        py"""
        $self.fitsourcebtncpca.clicked.connect(lambda: $self.open_file($self.fitsourcepathcpca,$self.fitsourcenamecpca))
        """
        self.overcheckconvpca = qwid.QCheckBox()
        formconvpca.addWidget(qwid.QLabel("Highest PC to study:"),1,0)
        formconvpca.addWidget(self.npchigh,1,1)
        formconvpca.addWidget(qwid.QLabel("Overwrite files?"),2,0)
        formconvpca.addWidget(self.overcheckconvpca,2,1)
        formconvpca.addWidget(self.convpcabtn,3,0,1,2)
        convpcabox.setObjectName("ColoredGroupBox")  
        convpcabox.setStyleSheet("$stylesheet")   



        # CV
        cvbox = qwid.QGroupBox("CV")
        formcv = qwid.QGridLayout()
        cvbox.setLayout(formcv)
        self.threshcv = qwid.QDoubleSpinBox()
        py"""$self.threshcv.setMinimum(float("-inf"))
        $self.threshcv.setMaximum(float("inf"))    """
        self.vshiftcv = qwid.QDoubleSpinBox()
        py"""$self.vshiftcv.setMinimum(float("-inf"))
        $self.vshiftcv.setMaximum(float("inf"))    """
        self.cvbtn = qwid.QPushButton("Go CV")
        self.cvbtn.setCheckable(false)
        self.cvbtn.setEnabled(true)
        py"""
        $self.cvbtn.clicked.connect(lambda: $self.callcv($self.npchigh.value()))
        """
        self.fitsourcenamecv = "FITSNAME"
        self.fitsourcepathcv = qwid.QLineEdit("FITSOURCE FOR CV")
        self.fitsourcebtncv = qwid.QPushButton("Browse")
        formcv.addWidget(qwid.QLabel("Enter the path to fits source : "),0,0)
        formcv.addWidget(self.fitsourcepathcv,0,1)
        formcv.addWidget(self.fitsourcebtncv,0,3) 
        self.namefitsourcecv="FITSNAME.fits"
        py"""
        $self.fitsourcebtncv.clicked.connect(lambda: $self.open_file($self.fitsourcepathcv,$self.fitsourcenamecv))
        """
        formcv.addWidget(qwid.QLabel("Intensity threshold (K): "),1,0)
        formcv.addWidget(self.threshcv,1,1)
        formcv.addWidget(qwid.QLabel("Shift in velocity (in same units as given above): "),2,0)
        formcv.addWidget(self.vshiftcv,2,1)
        self.overcheckcv = qwid.QCheckBox()
        formcv.addWidget(qwid.QLabel("Overwrite CV data?"),3,0)
        formcv.addWidget(self.overcheckcv,3,1)
        formcv.addWidget(self.cvbtn,4,0,1,2)
        cvbox.setObjectName("ColoredGroupBox")  
        cvbox.setStyleSheet("$stylesheet") 

        # CVI
        cvibox = qwid.QGroupBox("CVI")
        formcvi = qwid.QGridLayout()
        cvibox.setLayout(formcvi)
        self.cvibtn = qwid.QPushButton("Go CVI")
        self.cvibtn.setCheckable(false)
        self.cvibtn.setEnabled(true)
        py"""
        $self.cvibtn.clicked.connect(lambda: $self.callcvi($self.npchigh.value()))
        """
        self.fitsourcenamecvi = "FITSNAME"
        self.fitsourcepathcvi = qwid.QLineEdit("FITSOURCE FOR CVI")
        self.fitsourcebtncvi = qwid.QPushButton("Browse")
        formcvi.addWidget(qwid.QLabel("Enter the path to fits source (cv): "),0,0)
        formcvi.addWidget(self.fitsourcepathcvi,0,1)
        formcvi.addWidget(self.fitsourcebtncvi,0,3) 
        self.namefitsourcecvi="FITSNAME.fits"
        py"""
        $self.fitsourcebtncvi.clicked.connect(lambda: $self.open_file($self.fitsourcepathcvi,$self.fitsourcenamecvi))
        """
        self.lags = qwid.QLineEdit("2,3,4,10,15")
        self.difftype = qwid.QLineEdit("abs/relative")
        formcvi.addWidget(qwid.QLabel("Lags values (separate by a , )"),1,0)
        formcvi.addWidget(self.lags,1,1)
        formcvi.addWidget(qwid.QLabel("Type of differences : abs for absolute, or relative)"),2,0)
        formcvi.addWidget(self.difftype,2,1)
        self.overcheckcvi = qwid.QCheckBox()
        formcvi.addWidget(qwid.QLabel("Overwrite CVI data?"),3,0)
        formcvi.addWidget(self.overcheckcvi,3,1)
        formcvi.addWidget(self.cvibtn,4,0,1,2)
        cvibox.setObjectName("ColoredGroupBox")  
        cvibox.setStyleSheet("$stylesheet") 


        # Sp(l)
        splbox = qwid.QGroupBox("Structure Function")
        formspl = qwid.QGridLayout()
        splbox.setLayout(formspl)
        self.splbtn = qwid.QPushButton("Go Spl")
        self.splbtn.setCheckable(false)
        self.splbtn.setEnabled(true)
        py"""
        $self.splbtn.clicked.connect(lambda: $self.callspl($self.npchigh.value()))
        """
        self.fitsourcenamespl = "FITSNAME"
        self.fitsourcepathspl = qwid.QLineEdit("FITSOURCE FOR SP(l)")
        self.fitsourcebtnspl = qwid.QPushButton("Browse")
        formspl.addWidget(qwid.QLabel("Enter the path to fits source (cvi allangles): "),0,0)
        formspl.addWidget(self.fitsourcepathspl,0,1)
        formspl.addWidget(self.fitsourcebtnspl,0,3) 
        self.namefitsourcespl="FITSNAME.fits"
        py"""
        $self.fitsourcebtnspl.clicked.connect(lambda: $self.open_file($self.fitsourcepathspl,$self.fitsourcenamespl))
        """
        self.spl = qwid.QLineEdit("1,2,3,4,5,6")
        formspl.addWidget(qwid.QLabel("Orders to compute (separate by a , )"),1,0)
        formspl.addWidget(self.spl,1,1)
        self.overcheckspl = qwid.QCheckBox()
        formspl.addWidget(qwid.QLabel("Overwrite SPL data?"),2,0)
        formspl.addWidget(self.overcheckspl,2,1)
        formspl.addWidget(self.splbtn,3,0,1,2)
        splbox.setObjectName("ColoredGroupBox")  
        splbox.setStyleSheet("$stylesheet")  

        layout.addWidget(kinggroup,1,0,1,2)
        layout.addWidget(swobox,2,1,1,1.5)
        layout.addWidget(pcabox,2,0,1,1.5)
        layout.addWidget(convpcabox,3,0,1,2)
        layout.addWidget(cvbox,4,0)
        layout.addWidget(cvibox,4,1)
        layout.addWidget(splbox,5,0,1,2)


        widget = qwid.QWidget()
        widget.setLayout(layout)
        self.setCentralWidget(widget)
        #self.show()



    end 

    function open_file(self,paramname,filenameonly)
        filename, _ = py"""$qwid.QFileDialog.getOpenFileName(None,str("Open Fits file for PCA "),'',str("*.fits"))"""
        paramname.setText(filename)
        filenameonly = qpro.QUrl.fromLocalFile(filename)
    end
    
    function open_dir(self,paramname)
        pathname = py"""$qwid.QFileDialog.getExistingDirectory(None,str("Select saving folder"))"""
        paramname.setText(pathname)
    end

    function callpca(self) # THIS IS A SLOT
        Unveil.pca("")
    end

    function callswo(self,nn) # THIS IS A SLOT
        println("PCA running...")
        println("$nn")
    end

    function callconvpca(self)#,fitspath,filename,pathtosave,savename,noisecan,unitvel,highpc,blank,overwrite) # THIS IS A SLOT
        Unveilapp.convpca(self.fitsourcepathcpca.text(),self.fitsourcenamecpca,self.savepath.text(),self.savename.text(),[self.noisebtn1.value(),self.noisebtn2.value()],self.units.text(),self.npchigh.value(),self.blankvalbtn.value(),self.overcheckconvpca.isChecked())
    end

    function callcv(self,nn) # THIS IS A SLOT
        println("PCA running...")
        println("$nn")
    end

    function callcvi(self,nn) # THIS IS A SLOT
        println("PCA running...")
        println("$nn")
    end

    function callspl(self,nn) # THIS IS A SLOT
        println("PCA running...")
        println("$nn")
    end

end


app = qwid.QApplication(sys.argv)
wind = MainWindow()
wind.show()
sys.exit(app.exec())





