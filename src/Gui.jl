module Gui
# using Pkg
# Pkg.add("PyCall")
# Pkg.add("DelimitedFiles")
# Pkg.add("https://gricad-gitlab.univ-grenoble-alpes.fr/delcamps/unveil#Graphical_Interface")
using PyCall
using DelimitedFiles
using Unveil
#pyimport_conda("QtWidgets","PyQt5") 
#@pyimport PyQt5
@pyimport PyQt5
@pyimport sys
qwid = pyimport("PyQt5.QtWidgets")
qpro = pyimport("PyQt5.QtCore")
qgui = pyimport("PyQt5.QtGui")



@pydef mutable struct MainWindow <: qwid.QMainWindow 
    #class MainWindow(QMainWindow):#class MainWindow(QMainWindow):
    function __init__(self)  #def __init__(self):
        # copy!(qwid, pyimport("PyQt6.QtWidgets", "PyQt6"))
        # copy!(qpro, pyimport("PyQt6.QtCore", "PyQt6"))
        # copy!(qgui, pyimport("PyQt6.QtGui", "PyQt6"))

        py"""super($MainWindow,$self).__init__()"""
        #qwid.QMainWindow.__init__(self)
        #super(MainWindow,self).__init__()
       # self.ui()


    end


    function ui(self)
        self.scroll = qwid.QScrollArea() 
        backcolid = "#3b1eb4" #"#5D5D5D"
        bordsize = 2
        self.setStyleSheet("background-color: #F3FBFF")

        stylesheet = "QGroupBox {border: $(bordsize)px solid $(backcolid); background-color: white; font: bold 20px ; margin-top: 4.5ex;subcontrol-origin:  margin } QGroupBox::title{subcontrol-origin: margin ; subcontrol-position : top middle; padding: 25px 0px -10px 0px;  color: $(backcolid)} QWidget{background-color: white} QPushButton::pressed {background-color:black;} QPushButton::!enabled {background-color:#BABAC5; color:black} QPushButton::hover {outline-color: transparent; border: 1.5px solid; border-color: rgb(151, 195, 243);} QCheckBox {color : #FF4741}"
        self.setWindowTitle("Unveil")

        # GLOBAL LAYOUT, CONTAINS EVERYTHING
        layout = qwid.QGridLayout()

        #self.setStyleSheet("background-color: #DADAFC;")
        #self.setLayout(layout)


        #GROUP ESSENTIAL PARAMETER
        kinggroup = qwid.QGroupBox("Essential parameters")
        form_layout = qwid.QGridLayout()
        #form_layout.setStyleSheet("background-color: #DADAFC")

        #kinggroup.setObjectName("ColoredGroupBox")  
        kinggroup.setStyleSheet("QGroupBox{border: $(bordsize)px solid $(backcolid); background: white; font: bold 20px; margin-top: 3ex; subcontrol-origin: margin} QGroupBox::title{subcontrol-origin: margin; subcontrol-position : top middle; padding: 10px 0px -10px 0px; color: $(backcolid)} QWidget{background-color: white} QPushButton::pressed {background-color:black} QPushButton::!enabled {background-color:#BABAC5; color:black} QPushButton::hover{
            outline-color: transparent;
            border: 1px solid;
            border-color: rgb(151, 195, 243);} QCheckBox::checked {
            color : #FF4741}")  
        self.savepath = qwid.QLineEdit("")#("PATH TO SAVE DATA")
        self.savepathbtn = qwid.QPushButton("Browse")
        form_layout.addWidget(qwid.QLabel("Enter the path to saving directory : "),1,0)
        form_layout.addWidget(self.savepath,1,1)
        form_layout.addWidget(self.savepathbtn,1,3) 
        py"""
        $self.savepathbtn.clicked.connect(lambda: $self.open_dir($self.savepath))
        #$self.savepath.mousePressEvent($self.savepath.setText(""))#$self.mousePressed($self.savepath))
        """

        self.savename = qwid.QLineEdit("")#"NAME FOR SAVED DATA")
        form_layout.addWidget(qwid.QLabel("Enter the generic name to give to saved data WITHOUT SPACE: "),2,0)
        form_layout.addWidget(self.savename,2,1)
        self.noisebtn1 = qwid.QSpinBox()
        self.noisebtn1.setMinimum(1)
        self.noisebtn1.setMaximum(2e9)
        form_layout.addWidget(qwid.QLabel("Indices of first channel with noise "),3,0)
        form_layout.addWidget(self.noisebtn1,3,1)
        self.noisebtn2 = qwid.QSpinBox()
        self.noisebtn2.setMinimum(2)
        self.noisebtn2.setMaximum(2e9)
        form_layout.addWidget(qwid.QLabel("Indices of first channel with noise "),4,0)
        form_layout.addWidget(self.noisebtn2,4,1)
        self.blankvalbtn = qwid.QSpinBox()
        self.blankvalbtn.setMinimum(-1e9)
        self.blankvalbtn.setMaximum(1e9)   
        form_layout.addWidget(qwid.QLabel("Blank value of your data "),5,0)
        form_layout.addWidget(self.blankvalbtn,5,1)


        self.units = qwid.QLineEdit("")#"Units")
        form_layout.addWidget(qwid.QLabel("Units of the velocity axis (see header, should be m/s or km/s) "),6,0)
        form_layout.addWidget(self.units,6,1)
        kinggroup.setLayout(form_layout)





        # PCA GROUP
        pcabox = qwid.QGroupBox("PCA")
        formpca = qwid.QGridLayout()
        pcabox.setLayout(formpca)
        self.npc = qwid.QSpinBox()
        self.npc.setMinimum(1)
        self.npc.setMaximum(150)
        MainWindow.button_is_checked=true
        self.pcabtn = qwid.QPushButton("Go PCA")
        self.pcabtn.setCheckable(false)
        self.pcabtn.setEnabled(true)
        py"""
        $self.pcabtn.clicked.connect(lambda: $self.callpca())
        """
        self.overcheckpca = qwid.QCheckBox()
        self.fitsourcenamepca = "FITSNAME"
        self.fitsourcepathpca = qwid.QLineEdit("")#"FITSOURCE FOR PCA")
        self.fitsourcebtnpca = qwid.QPushButton("Browse")
        formpca.addWidget(qwid.QLabel("Enter the path to fits source : "),0,0)
        formpca.addWidget(self.fitsourcepathpca,0,1)
        formpca.addWidget(self.fitsourcebtnpca,0,3) 
        self.namefitsourcepca="FITSNAME.fits"
        py"""
        $self.fitsourcebtnpca.clicked.connect(lambda: $self.open_file($self.fitsourcepathpca,$self.fitsourcenamepca))
        """
        formpca.addWidget(qwid.QLabel("Number of PC:"),1,0)
        formpca.addWidget(self.npc,1,1)
        formpca.addWidget(qwid.QLabel("Overwrite PCA saved cube?"),2,0)
        formpca.addWidget(self.overcheckpca,2,1.5)
        formpca.addWidget(self.pcabtn,3,0,1,2)
        #pcabox.setObjectName("ColoredGroupBox")  
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
        
        py"""
        $self.swobtn.clicked.connect(lambda: $self.callswo())
        """
        self.fitsourcenameswo = "FITSNAME"
        self.fitsourcepathswo = qwid.QLineEdit("")#FITSOURCE FOR SWO")
        self.fitsourcebtnswo = qwid.QPushButton("Browse")
        formswo.addWidget(qwid.QLabel("Enter the path to fits source : "),0,0)
        formswo.addWidget(self.fitsourcepathswo,0,1)
        formswo.addWidget(self.fitsourcebtnswo,0,3) 
        self.namefitsourcepca="FITSNAME.fits"
        py"""
        $self.fitsourcebtnswo.clicked.connect(lambda: $self.open_file($self.fitsourcepathswo,$self.fitsourcenameswo))
        """
        self.step = qwid.QSpinBox()
        formswo.addWidget(qwid.QLabel("Square size pixel unit to use around false window results"))
        formswo.addWidget(self.step,1,1)
        self.overcheckswo = qwid.QCheckBox()
        formswo.addWidget(qwid.QLabel("Overwrite SWO saved cube?"),2,0)
        formswo.addWidget(self.overcheckswo,2,1)
        formswo.addWidget(self.swobtn,3,0,1,2)
        #swobox.setObjectName("ColoredGroupBox")  
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
        self.fitsourcepathcpca = qwid.QLineEdit("")#FITSOURCE FOR CONV PCA")
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
        #convpcabox.setObjectName("ColoredGroupBox")  
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
        $self.cvbtn.clicked.connect(lambda: $self.callcv())
        """
        self.fitsourcenamecv = "FITSNAME"
        self.fitsourcepathcv = qwid.QLineEdit("")#FITS FOR CV")
        self.fitsourcebtncv = qwid.QPushButton("Browse")
        formcv.addWidget(qwid.QLabel("Enter the path to the fits treated/reconstructed : "),0,0)
        formcv.addWidget(self.fitsourcepathcv,0,1)
        formcv.addWidget(self.fitsourcebtncv,0,3) 
        self.namefitsourcecv="FITSNAME.fits"
        py"""
        $self.fitsourcebtncv.clicked.connect(lambda: $self.open_file($self.fitsourcepathcv,$self.fitsourcenamecv))
        """
        self.fitsourcenameoricv = "FITSNAME"
        self.fitsourcepathoricv = qwid.QLineEdit("")#FITSOURCE FOR CV")
        self.fitsourcebtnoricv = qwid.QPushButton("Browse")
        formcv.addWidget(qwid.QLabel("Enter the path to the fits not treated/reconstructed : "),1,0)
        formcv.addWidget(self.fitsourcepathoricv,1,1)
        formcv.addWidget(self.fitsourcebtnoricv,1,3) 
        self.namefitsourcecv="FITSNAME.fits"
        py"""
        $self.fitsourcebtnoricv.clicked.connect(lambda: $self.open_file($self.fitsourcepathoricv,$self.fitsourcenameoricv))
        """
        formcv.addWidget(qwid.QLabel("Intensity threshold (K): "),2,0)
        formcv.addWidget(self.threshcv,2,1)
        formcv.addWidget(qwid.QLabel("Shift in velocity (in same units as given above): "),3,0)
        formcv.addWidget(self.vshiftcv,3,1)
        self.overcheckcv = qwid.QCheckBox()
        formcv.addWidget(qwid.QLabel("Overwrite CV data?"),4,0)
        formcv.addWidget(self.overcheckcv,4,1)
        formcv.addWidget(self.cvbtn,5,0,1,2)
        #cvbox.setObjectName("ColoredGroupBox")  
        cvbox.setStyleSheet("$stylesheet") 


        # CVI
        cvibox = qwid.QGroupBox("CVI")
        formcvi = qwid.QGridLayout()
        cvibox.setLayout(formcvi)
        self.cvibtn = qwid.QPushButton("Go CVI")
        self.cvibtn.setCheckable(false)
        self.cvibtn.setEnabled(true)
        py"""
        $self.cvibtn.clicked.connect(lambda: $self.callcvi())
        """
        self.fitsourcenamecvi = "FITSNAME"
        self.fitsourcepathcvi = qwid.QLineEdit("")#FITSOURCE FOR CVI")
        self.fitsourcebtncvi = qwid.QPushButton("Browse")
        formcvi.addWidget(qwid.QLabel("Enter the path to fits source (cv): "),0,0)
        formcvi.addWidget(self.fitsourcepathcvi,0,1)
        formcvi.addWidget(self.fitsourcebtncvi,0,3) 
        self.namefitsourcecvi="FITSNAME.fits"
        py"""
        $self.fitsourcebtncvi.clicked.connect(lambda: $self.open_file($self.fitsourcepathcvi,$self.fitsourcenamecvi))
        """
        self.lags = qwid.QLineEdit("2,3,4,10,15")
        self.difftype = qwid.QLineEdit("")#"abs/relative")
        formcvi.addWidget(qwid.QLabel("Lags values (separate by a , )"),1,0)
        formcvi.addWidget(self.lags,1,1)
        formcvi.addWidget(qwid.QLabel("Type of differences : abs for absolute, or relative"),2,0)
        formcvi.addWidget(self.difftype,2,1)
        self.overcheckcvi = qwid.QCheckBox()
        formcvi.addWidget(qwid.QLabel("Overwrite CVI data?"),3,0)
        formcvi.addWidget(self.overcheckcvi,3,1)
        formcvi.addWidget(self.cvibtn,4,0,1,2)
        #cvibox.setObjectName("ColoredGroupBox")  
        cvibox.setStyleSheet("$stylesheet") 



        # Sp(l)
        splbox = qwid.QGroupBox("Structure Function")
        formspl = qwid.QGridLayout()
        splbox.setLayout(formspl)
        self.splbtn = qwid.QPushButton("Go Spl")
        self.splbtn.setCheckable(false)
        self.splbtn.setEnabled(true)
        py"""
        $self.splbtn.clicked.connect(lambda: $self.callspl())
        """
        self.fitsourcenamespl = "FITSNAME"
        self.fitsourcepathspl = qwid.QLineEdit("")#FITSOURCE FOR SP(l)")
        self.fitsourcebtnspl = qwid.QPushButton("Browse")
        formspl.addWidget(qwid.QLabel("Enter the path to fits source (cvi allangles): "),0,0)
        formspl.addWidget(self.fitsourcepathspl,0,1)
        formspl.addWidget(self.fitsourcebtnspl,0,3) 
        self.namefitsourcespl="FITSNAME.fits"
        py"""
        $self.fitsourcebtnspl.clicked.connect(lambda: $self.open_file($self.fitsourcepathspl,$self.fitsourcenamespl))
        """
        self.orders = qwid.QLineEdit("1,2,3,4,5,6")
        formspl.addWidget(qwid.QLabel("Orders to compute (separate by a , )"),1,0)
        formspl.addWidget(self.orders,1,1)
        self.method = qwid.QLineEdit("moninyaglom/hily")
        formspl.addWidget(qwid.QLabel("Which method to use? moninyaglom/hily"),2,0)
        formspl.addWidget(self.method,2,1)
        self.overcheckspl = qwid.QCheckBox()
        formspl.addWidget(qwid.QLabel("Overwrite SPL data?"),3,0)
        formspl.addWidget(self.overcheckspl,3,1)
        formspl.addWidget(self.splbtn,4,0,1,2)
        #splbox.setObjectName("ColoredGroupBox")  
        splbox.setStyleSheet("$stylesheet")  

        # SAVE AND IMPORT PARAMETER FILE
        self.saveparambtn = qwid.QPushButton("Export parameters into a file")
        self.importparambtn = qwid.QPushButton("Import parameters from a file")
        self.saveparambtn2 = qwid.QPushButton("Export parameters into a file")
        self.importparambtn2 = qwid.QPushButton("Import parameters from a file")
        py"""
        $self.saveparambtn.clicked.connect(lambda: $self.exportfile())
        $self.importparambtn.clicked.connect(lambda: $self.importfile())
        $self.saveparambtn2.clicked.connect(lambda: $self.exportfile())
        $self.importparambtn2.clicked.connect(lambda: $self.importfile())
        """

        self.exitbtn = qwid.QPushButton("Exit")
        self.exitbtn.setStyleSheet("QPushButton {color:#FF4741} QPushButton::hover {outline-color: transparent; border: 1.5px solid; border-color: #FF4741;}")
        py"""
        $self.exitbtn.clicked.connect(lambda: $self.close())
        """
        self.saveparambtn2.setStyleSheet("QPushButton::pressed {background-color:black;} QPushButton::!enabled {background-color:#BABAC5; color:black} QPushButton::hover {outline-color: transparent; border: 1.5px solid; border-color: rgb(151, 195, 243);}")
        self.importparambtn2.setStyleSheet("QPushButton::pressed {background-color:black;} QPushButton::!enabled {background-color:#BABAC5; color:black} QPushButton::hover {outline-color: transparent; border: 1.5px solid; border-color: rgb(151, 195, 243);}")
        self.saveparambtn.setStyleSheet("QPushButton::pressed {background-color:black;} QPushButton::!enabled {background-color:#BABAC5; color:black} QPushButton::hover {outline-color: transparent; border: 1.5px solid; border-color: rgb(151, 195, 243);}")
        self.importparambtn.setStyleSheet("QPushButton::pressed {background-color:black;} QPushButton::!enabled {background-color:#BABAC5; color:black} QPushButton::hover {outline-color: transparent; border: 1.5px solid; border-color: rgb(151, 195, 243);}")


        # ADD EVERY WIDGET ON THE GLOBAL LAYOUT

        layout.addWidget(self.saveparambtn2,0,0,1,1.5)
        layout.addWidget(self.importparambtn2,0,1.5,1,1.5)
        layout.addWidget(kinggroup,1,0,1,2)
        layout.addWidget(swobox,2,1,1,1.5)
        layout.addWidget(pcabox,2,0,1,1.5)
        layout.addWidget(convpcabox,3,0,1,2)
        layout.addWidget(cvbox,4,0)
        layout.addWidget(cvibox,4,1)
        layout.addWidget(splbox,5,0,1,2)
        layout.addWidget(self.saveparambtn,6,0,1,1.5)
        layout.addWidget(self.importparambtn,6,1.5,1,1.5)
        layout.addWidget(self.exitbtn,7,0,1,3)
        self.scroll.setWidgetResizable(true)
        self.setCentralWidget(self.scroll)
        self.setGeometry(600, 100, 1300, 800)

        widget = qwid.QWidget()

        widget.setLayout(layout)
        self.scroll.setWidget(widget)


        #self.show()



    end 

    function mousePressed(self,paramname)
        paramname.setText("")
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

    function importfile(self)
        pathimportfile,_ = py"""$qwid.QFileDialog.getOpenFileName(None,str("Select parameter file to import"),'',str("*.txt"))"""
        self.readparamfile(pathimportfile)
    end

    function readparamfile(self,pathfile)
        file = readdlm("$pathfile",comments=true)
        self.savepath.setText("$(file[1,2])") 
        self.savename.setText("$(file[2,2])") 
        self.noisebtn1.setValue(file[3,2]) 
        self.noisebtn2.setValue(file[4,2]) 
        self.blankvalbtn.setValue(file[5,2]) 
        self.units.setText("$(file[6,2])") 
        self.fitsourcepathpca.setText("$(file[7,2])") 
        self.npc.setValue(file[8,2]) 
        self.fitsourcepathswo.setText("$(file[9,2])") 
        self.fitsourcepathcpca.setText("$(file[10,2])") 
        self.npchigh.setValue(file[11,2]) 
        self.fitsourcepathcv.setText("$(file[12,2])") 
        self.fitsourcepathoricv.setText("$(file[13,2])") 
        self.threshcv.setValue(file[14,2]) 
        self.vshiftcv.setValue(file[15,2]) 
        self.fitsourcepathcvi.setText("$(file[16,2])") 
        self.lags.setText("$(file[17,2])") 
        self.difftype.setText("$(file[18,2])") 
        self.fitsourcepathspl.setText("$(file[19,2])") 
        self.orders.setText("$(file[20,2])") 
        self.method.setText("$(file[21,2])") 


    end

    function exportfile(self)
        filename,_ = py"""$qwid.QFileDialog.getSaveFileName(None,str("Save parameter file"),str($self.savename.text()+ ".txt"),str("*.txt"))"""
        self.saveparamfile(filename)
    end 

    function saveparamfile(self,path)
        open("$path","w") do io
            towrite = [
                "# ESSENTIAL PARAMETERS",
                "SAVINGDIR    \"$(self.savepath.text())\"",
                "SAVENAME     \"$(self.savename.text())\"",
                "FIRSTNOISE    $(self.noisebtn1.value())",
                "SECONDNOISE    $(self.noisebtn2.value())",
                "BLANKVALUE     $(self.blankvalbtn.value())",
                "UNIT           \"$(self.units.text())\"",
                "#",
                "# PCA",
                "FITSOURCE   \"$(self.fitsourcepathpca.text())\"",
                "PCNUMBER    $(self.npc.value())",
                "#",
                "# SWO",
                "FITSOURCE   \"$(self.fitsourcepathswo.text())\"",
                "#",
                "# CONVERGENCE OF THE PCA",
                "FITSOURCE   \"$(self.fitsourcepathcpca.text())\"",
                "HIGHESTPC   $(self.npchigh.value())",
                "#",
                "# CV",
                "FITSFORCV   \"$(self.fitsourcepathcv.text())\"",
                "FITSOURCE  \"$(self.fitsourcepathoricv.text())\"",
                "INTHRESHOLD    $(self.threshcv.value())",
                "VELOCITYSHIFT  $(self.vshiftcv.value())",
                "#",
                "# CVI",
                "FITSFORCVI   \"$(self.fitsourcepathcvi.text())\"",
                "LAGS        \"$(self.lags.text())\"",
                "DIFFTYPE    \"$(self.difftype.text())\"",
                "#",
                "# STRUCTUREFUNCTION",
                "FITSFORSPL   \"$(self.fitsourcepathspl.text())\"",
                "ORDERS        \"$(self.orders.text())\"",
                "METHOD        \"$(self.method.text())\"",
            ]



            writedlm(io,towrite,quotes=false)
        end
    end

    function callpca(self) # THIS IS A SLOT
        self.pcabtn.setEnabled(false)
        self.pcabtn.setText("PCA running....")
        qwid.QApplication.processEvents()
        Unveil.pca(
            self.fitsourcepathpca.text(),
            self.fitsourcenamepca,
            self.savepath.text(),
            self.savename.text(),
            self.units.text(),
            self.npc.value(),
            self.blankvalbtn.value(),
            self.overcheckpca.isChecked()
        )
        self.pcabtn.setEnabled(true)
        self.pcabtn.setText("Go PCA")
        qwid.QApplication.processEvents()
    end
    
    # function changebtnfa(self)
    #     self.pcabtn.setEnabled(false)
    # end
    # function changebtntr(self)
    #     self.pcabtn.setEnabled(true)
    # end
    function callswo(self) # THIS IS A SLOT
        self.swobtn.setEnabled(false)
        self.swobtn.setText("SWO running....")
        qwid.QApplication.processEvents()
        Unveil.swo(
            self.fitsourcepathswo.text(),
            self.fitsourcenameswo,
            self.savepath.text(),
            self.savename.text(),
            self.units.text(),
            self.step.value(),
            self.blankvalbtn.value(),
            [self.noisebtn1.value(),self.noisebtn2.value()],
            self.overcheckswo.isChecked()
        )
        self.swobtn.setEnabled(true)
        self.swobtn.setText("Go SWO")
        qwid.QApplication.processEvents()


    end

    function callconvpca(self)#,fitspath,filename,pathtosave,savename,noisecan,unitvel,highpc,blank,overwrite) # THIS IS A SLOT.
        self.convpcabtn.setEnabled(false)
        self.convpcabtn.setText("ConvPCA running....")
        qwid.QApplication.processEvents()
        Unveil.convpca(
            self.fitsourcepathcpca.text(),
            self.fitsourcenamecpca,
            self.savepath.text(),
            self.savename.text(),
            [self.noisebtn1.value(),self.noisebtn2.value()],
            self.units.text(),
            self.npchigh.value(),
            self.blankvalbtn.value(),
            self.overcheckconvpca.isChecked()
        )
        self.convpcabtn.setEnabled(true)
        self.convpcabtn.setText("Go ConvPCA")
        qwid.QApplication.processEvents()
    end


    function callcv(self) # THIS IS A SLOT
        self.cvbtn.setEnabled(false)
        self.cvbtn.setText("CV running....")
        qwid.QApplication.processEvents()
        Unveil.cv( 
            self.fitsourcepathcv.text(),
            self.fitsourcenamecv,
            self.savepath.text(),
            self.fitsourcepathoricv.text(),
            self.savename.text(),
            self.units.text(),
            self.threshcv.value(),
            [self.noisebtn1.value(),self.noisebtn2.value()],
            self.vshiftcv.value(),
            self.blankvalbtn.value(),
            self.overcheckcv.isChecked()
            )
        self.cvbtn.setEnabled(true)
        self.cvbtn.setText("Go CV")
        qwid.QApplication.processEvents()
    end



    function callcvi(self) # THIS IS A SLOT
        self.cvibtn.setEnabled(false)
        self.cvibtn.setText("CVI running....")
        qwid.QApplication.processEvents()
        Unveil.cvi(
            self.fitsourcepathcvi.text(),
            self.fitsourcenamecvi,
            self.savepath.text(),
            self.savename.text(),
            self.blankvalbtn.value(),
            self.lags.text(),
            self.difftype.text(),
            self.overcheckcvi.isChecked()
        )
        self.cvibtn.setEnabled(true)
        self.cvibtn.setText("Go CVI")
        qwid.QApplication.processEvents()
    end



    function callspl(self) # THIS IS A SLOT
        self.splbtn.setEnabled(false)
        self.splbtn.setText("Spl running....")
        Unveil.structure_functions(
            self.fitsourcepathspl.text(),
            self.fitsourcenamespl,
            self.savepath.text(),
            self.savename.text(),
            self.orders.text(),
            self.method.text(),
            self.blankvalbtn.value(),
            self.overcheckspl.isChecked()
        )
        self.splbtn.setEnabled(true)
        self.splbtn.setText("Go Spl")
        qwid.QApplication.processEvents()
    end




end


function opengui()
app = qwid.QApplication(sys.argv)
wind = MainWindow()
wind.show()
app.exec()
end




end