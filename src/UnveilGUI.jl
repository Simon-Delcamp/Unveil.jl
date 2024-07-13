module UnveilGUI


include("Dataprep.jl")              # Read and write fits
include("Analysis.jl")              # PCA metric, power_spectra, rms...
include("PCA.jl")                   # Functions for method PCA
include("SWO.jl")                   # Functions for method SWO
include("Graphic.jl")               # Plotting 
include("CVI.jl")                   # Functions for CVI computations
include("Structure_functions.jl")   # Functions for Structure functions computations


using .Dataprep
using .SWO
using .Graphic
using .Analysis
using .Structure_functions
using .CVI
using .PCA
using Plots
using ProgressBars
using StatsBase

export pca
export swo
export convpca
export cvcvi
export cvi 
export multipca
export combinecv
export prodvarfile
export structure_functions
export compmethod_stcfct

using PyCall

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
        UnveilGUI.pca(
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
        # Unveil.swo(
        #     self.fitsourcepathswo.text(),
        #     self.fitsourcenameswo,
        #     self.savepath.text(),
        #     self.savename.text(),
        #     self.units.text(),
        #     self.step.value(),
        #     self.blankvalbtn.value(),
        #     [self.noisebtn1.value(),self.noisebtn2.value()],
        #     self.overcheckswo.isChecked()
        # )
        self.swobtn.setEnabled(true)
        self.swobtn.setText("Go SWO")
        qwid.QApplication.processEvents()


    end

    function callconvpca(self)#,fitspath,filename,pathtosave,savename,noisecan,unitvel,highpc,blank,overwrite) # THIS IS A SLOT.
        self.convpcabtn.setEnabled(false)
        self.convpcabtn.setText("ConvPCA running....")
        qwid.QApplication.processEvents()
        # Unveil.convpca(
        #     self.fitsourcepathcpca.text(),
        #     self.fitsourcenamecpca,
        #     self.savepath.text(),
        #     self.savename.text(),
        #     [self.noisebtn1.value(),self.noisebtn2.value()],
        #     self.units.text(),
        #     self.npchigh.value(),
        #     self.blankvalbtn.value(),
        #     self.overcheckconvpca.isChecked()
        # )
        self.convpcabtn.setEnabled(true)
        self.convpcabtn.setText("Go ConvPCA")
        qwid.QApplication.processEvents()
    end


    function callcv(self) # THIS IS A SLOT
        self.cvbtn.setEnabled(false)
        self.cvbtn.setText("CV running....")
        qwid.QApplication.processEvents()
        # Unveil.cv( 
        #     self.fitsourcepathcv.text(),
        #     self.fitsourcenamecv,
        #     self.savepath.text(),
        #     self.fitsourcepathoricv.text(),
        #     self.savename.text(),
        #     self.units.text(),
        #     self.threshcv.value(),
        #     [self.noisebtn1.value(),self.noisebtn2.value()],
        #     self.vshiftcv.value(),
        #     self.blankvalbtn.value(),
        #     self.overcheckcv.isChecked()
        #     )
        self.cvbtn.setEnabled(true)
        self.cvbtn.setText("Go CV")
        qwid.QApplication.processEvents()
    end



    function callcvi(self) # THIS IS A SLOT
        self.cvibtn.setEnabled(false)
        self.cvibtn.setText("CVI running....")
        qwid.QApplication.processEvents()
        # Unveil.cvi(
        #     self.fitsourcepathcvi.text(),
        #     self.fitsourcenamecvi,
        #     self.savepath.text(),
        #     self.savename.text(),
        #     self.blankvalbtn.value(),
        #     self.lags.text(),
        #     self.difftype.text(),
        #     self.overcheckcvi.isChecked()
        # )
        self.cvibtn.setEnabled(true)
        self.cvibtn.setText("Go CVI")
        qwid.QApplication.processEvents()
    end



    function callspl(self) # THIS IS A SLOT
        self.splbtn.setEnabled(false)
        self.splbtn.setText("Spl running....")
        # Unveil.structure_functions(
        #     self.fitsourcepathspl.text(),
        #     self.fitsourcenamespl,
        #     self.savepath.text(),
        #     self.savename.text(),
        #     self.orders.text(),
        #     self.method.text(),
        #     self.blankvalbtn.value(),
        #     self.overcheckspl.isChecked()
        # )
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






#########################################################################
###             FOR GUI             ###
#########################################################################





# function convpca(FITSPATH,FILENAME,PATHTOSAVE,SAVENAME,NOISECAN,UNITVELOCITY,HIGHESTPC,BLANK,OVERWRITE; plot=true)

#     # Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
#     cube,VELOCITYVECTOR,DATADIMENSION,VELOCITYINCREMENT,HEAD = Dataprep.read_fits_ppv("$(FITSPATH)",UNITVELOCITY ; check=false)

#     # Prepare directories where plots and data will be saved.
#     Dataprep.directory_prep(PATHTOSAVE)



#     # Replace any NaN value into a missing value and deleted them (can't do PCA on missing values with that package).
#     cube = Dataprep.replace_nantomissing(cube)

#     ismis = 0
#     if any(ismissing,cube) 
#         ismis = 1
#         cube,missingplaces1D,missingplaces2D  = Dataprep.pca_prep(cube,DATADIMENSION)
#         cube                                  = convert(Array{Float64},cube)
#         DATADIMENSION_NOMISSING               = Dataprep.read_dim(cube)
#     else
#         cube                                 = reshape(cube,DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
#         cube                                 = convert(Array{Float64},cube)
#         DATADIMENSION_NOMISSING              = (DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
#     end

#     SIGMAT = Analysis.rms_cube(cube,NOISECAN)[2]

#     #First PCA
#     println("Perform PCA")
#     M, Yt, VARPERCENT,cubereconstructed = PCA.pca(cube,HIGHESTPC)

#     # Projection matrix
#     proj = PCA.proj(M)
#     mom1,mom2,mom3,mom4 = Analysis.fourmoments(proj,dim=2)


#     if ismis == 1
#         proj = Dataprep.addblank(proj,missingplaces2D[:,1:HIGHESTPC],BLANK,(DATADIMENSION[1],DATADIMENSION[2],HIGHESTPC))
#     end
#     proj = reshape(proj,(DATADIMENSION[1],DATADIMENSION[2],HIGHESTPC))
#     Dataprep.write_fits("$(FITSPATH)","$(SAVENAME)_projectionmatrix","$(PATHTOSAVE)/Data/",proj,(DATADIMENSION_NOMISSING[1],DATADIMENSION_NOMISSING[2],HIGHESTPC),BLANK,finished=true,overwrite=OVERWRITE)
#     xvector = range(1,HIGHESTPC)#[1:HIGHESTPC]



#     proj = 0
#     cubereconstructed = 0
#     cube = 0
#     missingplaces1D = 0.0 
#     missingplaces2D = 0
#     M = 0
#     Yt = 0
#     GC.gc()  
#     newname = "$(SAVENAME)_mom"
#     if (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/$(SAVENAME)_mom.pdf")==true)
#         println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
#         count = 0
#         for ix=1:size((findall.("$(SAVENAME)_mom",readdir("$(PATHTOSAVE)/Figures/"))))[1]
#             if size(findall("$(SAVENAME)_mom",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
#                 count += 1
#             end
#         end
#         newname = "$(newname)_$(count)"
#     end 


#     metric = Analysis.metricPCA(mom1,mom2,mom3,mom4,abs(VELOCITYINCREMENT))#,SIGMAT)#abs(VELOCITYINCREMENT))
#     #metric = Analysis.metricPCA(mom1,mom2,mom3,mom4)#,SIGMAT)#abs(VELOCITYINCREMENT))
#     #metric = Analysis.metricPCA(mom1,mom2,mom3,mom4,abs(VELOCITYINCREMENT))
#     #println(metric)
#     println("Metric calculated")
#     #BL#Graphic.distribcv_multipc(mom1[2:end-1],mom2[2:end-1],mom3[2:end-1],mom4[2:end-1],metric[2:end-1],xvector[2:end-1])
#     newname = "$(SAVENAME)_metric"
#     if (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/$(SAVENAME)_metric.pdf")==true)
#         println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
#         count = 0
#         for ix=1:size((findall.("$(SAVENAME)_metric",readdir("$(PATHTOSAVE)/Figures/"))))[1]
#             if size(findall("$(SAVENAME)_metric",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
#                 count += 1
#             end
#         end
#         newname = "$(newname)_$(count)"

#     end 
#     Dataprep.write_dat([metric mom1./0.05 mom2 mom3 mom4.-3 xvector],"$PATHTOSAVE/Data/","$(SAVENAME)_metricPCA",overwrite=OVERWRITE,more=["$FILENAME","Metric  Mom1   Mom2   Mom3   Mom4   PCs"])

# end #convpca





# function cv(FITSPATH,FITSNAME,PATHTOSAVE,FITSOURCE,SAVENAME,UNITVELOCITY,THRESHOLD,NOISECAN,VSHIFT,BLANK,OVERWRITE)
#     # Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
#     cube,VELOCITYVECTOR,DATADIMENSION,VELOCITYINCREMENT,HEAD = Dataprep.read_fits_ppv("$FITSPATH",UNITVELOCITY ; check=false)
#     if haskey(HEAD,"METHOD")==1 
#         if HEAD["METHOD"]=="PCA"
#             METH = HEAD["NBPC"]
#             METH = "$(METH)PC"
#         elseif HEAD["METHOD"]=="SWO"
#             METH = "SWO"
#         end
#     else
#         METH = "raw"
#     end


#     # Prepare directories where plots and data will be saved.
#     Dataprep.directory_prep(PATHTOSAVE)


#     cubesource = Dataprep.read_fits_ppv("$FITSOURCE",UNITVELOCITY ; check=false)[1]


#     # If no threshold, change it to the blank value and SIGMAT to 1, because function moment_one_field will blank every values lower than SIGMAT*THRESHOLD
#     if THRESHOLD==0
#         THRESHOLD = BLANK
#         SIGMAT = 1
#     else
#         cubesource = Dataprep.replace_nantomissing(cubesource)
#         cubesource = Dataprep.replace_blanktomissing(cubesource,BLANK)
#         cubesource = Dataprep.pca_prep(cubesource,DATADIMENSION)[1]
#         cubesource = convert(Array{Float64},cubesource)

#         SIGMAT     = Analysis.rms_cube(cubesource,NOISECAN)[2]
#     end

#     # Replace any NaN value into a missing value and deleted them (can't do PCA on missing values with that package).
#     cube = Dataprep.replace_nantomissing(cube)

#     ismis = 0
#     if any(ismissing,cube) 
#         ismis = 1
#         cube,missingplaces1D,missingplaces2D  = Dataprep.pca_prep(cube,DATADIMENSION)
#         cube                                  = convert(Array{Float64},cube)
#         DATADIMENSION_NOMISSING               = Dataprep.read_dim(cube)
#     else
#         cube                                 = reshape(cube,DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
#         cube                                 = convert(Array{Float64},cube)
#         DATADIMENSION_NOMISSING              = (DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
#     end

#     println("------ CV CALCULATION ------")
#     VELOCITYVECTOR = Dataprep.shiftspec(VELOCITYVECTOR,VSHIFT)
#     cvmap = CVI.moment_one_field(cube,SIGMAT,THRESHOLD,VELOCITYVECTOR,BLANK) # Calculate the first velocity moment order on data reconstructed
#     cube  = 0.0 
#     GC.gc()

#     if ismis == 1
#         cvmap = Dataprep.addblank(cvmap,missingplaces2D[:,1],BLANK,(DATADIMENSION[1],DATADIMENSION[2]))
#     end
#     cvmap = reshape(cvmap,(DATADIMENSION[1],DATADIMENSION[2]))

#     cvmap .= cvmap.+VSHIFT
#     VELOCITYVECTOR = Dataprep.shiftspec(VELOCITYVECTOR,-VSHIFT)

#     Dataprep.write_fits("$(FITSPATH)","CV_$(SAVENAME)_$(METH)","$(PATHTOSAVE)/Data/",cvmap,(DATADIMENSION_NOMISSING[1],DATADIMENSION_NOMISSING[2]),BLANK,finished=true,overwrite=OVERWRITE,more=["THRESH",THRESHOLD])
#     #cvmap = 0.0

#     println("CV map saved in $(PATHTOSAVE)/Data/CV_$(SAVENAME)_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")


# end #cv





# function cvi(FITSPATH,FITSNAME,PATHTOSAVE,SAVENAME,BLANK,LAG,DIFFTYPE,OVERWRITE)
#     if length(LAG)!=1
#         LAG = [parse(Int, ss) for ss in split(LAG,",")]
#         mult = true
#     else
#         mult = false
#     end 

#     cvmap,HEAD,DATADIMENSION = Dataprep.read_fits_pp("$FITSPATH")
#     if haskey(HEAD,"METHOD")==1 
#         if HEAD["METHOD"]=="PCA"
#             METH = HEAD["NBPC"]
#             METH = "$(METH)PC"
#         elseif HEAD["METHOD"]=="SWO"
#             METH = "SWO"
#         end
#     else
#         METH = "raw"
#     end

#     # Prepare directories where plots and data will be saved.
#     Dataprep.directory_prep(PATHTOSAVE)

#     cvmap = Dataprep.replace_nantomissing(cvmap)   
#     cvmap = Dataprep.replace_blanktomissing(cvmap,BLANK)

#     println("------ CVI CALCULATION ------")
#     if DIFFTYPE=="relative" 
#         cviallangle,cvimap_averaged,NANGLE = CVI.construct_cvimap(cvmap,LAG,(DATADIMENSION[1],DATADIMENSION[2]),diff="relative")
#         DIFFTYPE = "rel"

#     elseif DIFFTYPE=="abs" 
#         cviallangle,cvimap_averaged,NANGLE = CVI.construct_cvimap(cvmap,LAG,(DATADIMENSION[1],DATADIMENSION[2]),diff="absolute")
#         DIFFTYPE = "abs"

#     else
#         error("Not good argument in DIFFTYPE (should be abs or relative)")
#     end
#     cvmap = 0
#     GC.gc()
#     if  mult==true
#         cvimap_averaged = reshape(cvimap_averaged,DATADIMENSION[1],DATADIMENSION[2],size(LAG)[1])
#         for lag=1:size(LAG)[1]
#             cvimap_averaged[:,1:LAG[lag],lag] .= missing
#             cvimap_averaged[1:LAG[lag],:,lag] .= missing
#             cvimap_averaged[(DATADIMENSION[1]-LAG[lag])+1:DATADIMENSION[1],:,lag] .= missing
#             cvimap_averaged[:,(DATADIMENSION[2]-LAG[lag])+1:DATADIMENSION[2],lag] .= missing
#         end
#         cvimap_averaged = Dataprep.replace_missingtoblank(cvimap_averaged,BLANK)
#         cvimap_averaged = Dataprep.blank_equal(cvimap_averaged,0.0,BLANK)
#         cvimap_averaged = convert(Array{Float64},cvimap_averaged)

#         cviallangle = reshape(cviallangle,DATADIMENSION[1],DATADIMENSION[2],maximum(NANGLE),size(LAG)[1])
#         for lag=1:size(LAG)[1]
#             cviallangle[:,1:LAG[lag],:,lag] .= missing
#             cviallangle[1:LAG[lag],:,:,lag] .= missing
#             cviallangle[(DATADIMENSION[1]-LAG[lag])+1:DATADIMENSION[1],:,:,lag] .= missing
#             cviallangle[:,(DATADIMENSION[2]-LAG[lag])+1:DATADIMENSION[2],:,lag] .= missing
#         end
#         cviallangle = reshape(cviallangle,DATADIMENSION[1]*DATADIMENSION[2],maximum(NANGLE),size(LAG)[1])
#         cviallangle = Dataprep.replace_nantoblank(cviallangle,BLANK)
#         cviallangle = Dataprep.replace_missingtoblank(cviallangle,BLANK)
#         cviallangle = Dataprep.blank_equal(cviallangle,0.0,BLANK)
#         cviallangle = convert(Array{Float64},cviallangle)

#         # REDUCING THE SIZE OF THE ALLANGLE CVI FITS BY REMOVING SOME BLANKING VALUES AND CREATING A 2D ARRAY INSTEAD OF A 3D
#         sizee = Array{Float64}(undef,size(LAG)[1])
#         for lx=1:size(LAG)[1]
#             temp     = size(Dataprep.delete_allnotvalue(cviallangle[:,:,lx],BLANK))[1] |> Int64 #Data without missing value
#             sizee[lx]=temp
#         end 
#         maxi = maximum(sizee) |> Int64
#         cviallanglereduced = Array{Float64}(undef,maxi,size(LAG)[1])
#         cviallanglereduced .= BLANK
#         for lx=1:size(LAG)[1]
#             tr = sizee[lx] |> Int64
#             cviallanglereduced[1:tr,lx] .= Dataprep.delete_allnotvalue(cviallangle[:,:,lx],BLANK)
#         end
#         println("Start saving")
#         Dataprep.write_fits("$(FITSPATH)","CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)","$(PATHTOSAVE)/Data/",cvimap_averaged,(DATADIMENSION[1],DATADIMENSION[2],size(LAG)[1]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG])
#         println("CVI map saved in $(PATHTOSAVE)/Data/CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")
#        #Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)","$(PATHTOSAVE)/Data",cviallangle,(DATADIMENSION[1]*DATADIMENSION[2],maximum(NANGLE),size(LAG)[1]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG],cvi=true)
#         Dataprep.write_fits("$(FITSPATH)","CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)","$(PATHTOSAVE)/Data",cviallanglereduced,(maxi,size(LAG)[1]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG],cvi=true)
#         println("CVI map with all angles values saved in the $(PATHTOSAVE)/Data/CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")
    
#     else
#         cvimap_averaged = reshape(cvimap_averaged,DATADIMENSION[1],DATADIMENSION[2])
#         cvimap_averaged[:,1:LAG] .= missing
#         cvimap_averaged[1:LAG,:] .= missing
#         cvimap_averaged[(DATADIMENSION[1]-LAG)+1:DATADIMENSION[1],:] .= missing
#         cvimap_averaged[:,(DATADIMENSION[2]-LAG)+1:DATADIMENSION[2]] .= missing
#         cvimap_averaged = Dataprep.replace_missingtoblank(cvimap_averaged,BLANK)
#         cvimap_averaged = Dataprep.blank_equal(cvimap_averaged,0.0,BLANK)
#         cvimap_averaged = convert(Array{Float64},cvimap_averaged)

#         cviallangle = reshape(cviallangle,DATADIMENSION[1],DATADIMENSION[2],maximum(NANGLE))
#         cviallangle[:,1:LAG,:] .= missing
#         cviallangle[1:LAG,:,:] .= missing
#         cviallangle[(DATADIMENSION[1]-LAG)+1:DATADIMENSION[1],:,:] .= missing
#         cviallangle[:,(DATADIMENSION[2]-LAG)+1:DATADIMENSION[2],:] .= missing
        
#         cviallangle = reshape(cviallangle,DATADIMENSION[1]*DATADIMENSION[2],maximum(NANGLE))
#         cviallangle = Dataprep.replace_nantoblank(cviallangle,BLANK)
#         cviallangle = Dataprep.replace_missingtoblank(cviallangle,BLANK)
#         cviallangle = Dataprep.blank_equal(cviallangle,0.0,BLANK)
#         cviallangle = convert(Array{Float64},cviallangle)
#         #cviallangle = Dataprep.delete_allnotvalue(cviallangle,BLANK)

#         # REDUCING THE SIZE OF THE ALLANGLE CVI FITS BY REMOVING SOME BLANKING VALUES AND CREATING A 2D ARRAY INSTEAD OF A 3D
#         temp     = size(Dataprep.delete_allnotvalue(cviallangle[:,:],BLANK))[1] |> Int64 #Data without missing value
#         sizee=temp
 
#         maxi = maximum(sizee) |> Int64
#         cviallanglereduced = Array{Float16}(undef,maxi)
#         cviallanglereduced .= BLANK
#         tr = sizee |> Int64
#         cviallanglereduced[1:tr] .= Dataprep.delete_allnotvalue(cviallangle[:,:],BLANK)
#         println("Start saving")
#         Dataprep.write_fits("$(FITSPATH)","CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)","$(PATHTOSAVE)/Data/",cvimap_averaged,(DATADIMENSION[1],DATADIMENSION[2]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG])
#         println("CVI map saved in $(PATHTOSAVE)/Data/CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")

#         Dataprep.write_fits("$(FITSPATH)","CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)","$(PATHTOSAVE)/Data",cviallangle,(DATADIMENSION[1]*DATADIMENSION[2],size(cviallangle)[2]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG],cvi=true)
#         println("CVI map with all angles values saved in the $(PATHTOSAVE)/Data/CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")
#     end

# end #function cvi



# function pca(FITSPATH,FILENAME,PATHTOSAVE,SAVENAME,UNITVELOCITY,NBPC,BLANK,OVERWRITE)
#     (NBPC == 0) && (NBPC="raw")


#     # Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
#     cube,VELOCITYVECTOR,DATADIMENSION,VELOCITYINCREMENT,HEAD = Dataprep.read_fits_ppv("$(FITSPATH)",UNITVELOCITY ; check=false)

#     # Prepare directories where plots and data will be saved.
#     Dataprep.directory_prep(PATHTOSAVE)

#     # Replace any NaN value into a missing value and deleted them (can't do PCA on missing values with that package)
#     cube = Dataprep.replace_nantomissing(cube)
#     cube = Dataprep.replace_blanktomissing(cube,BLANK)
    
#     #cube = Dataprep.replace_nosignal(cube,DATADIMENSION,VELOCITYVECTOR,BLANK,SIGMAMAP)


#     ismis = 0
#     if any(ismissing,cube) 
#         ismis = 1
#         cube,missingplaces1D,missingplaces2D  = Dataprep.pca_prep(cube,DATADIMENSION)
#         cube                                  = convert(Array{Float64},cube)
#         DATADIMENSION_NOMISSING               = Dataprep.read_dim(cube)
#     else
#         cube                                 = reshape(cube,DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
#         cube                                 = convert(Array{Float64},cube)
#         DATADIMENSION_NOMISSING              = (DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
#     end

#     # Perform the first PCA analysis (same notation as in the MultivariateStats doc)
#     println("Perform PCA")
#     M, Yt, VARPERCENT,cubereconstructed = PCA.pca(cube,NBPC)
#     Dataprep.write_fits("$(FITSPATH)","Yt_$(NBPC)PC","$PATHTOSAVE/Data/",Yt,(NBPC,DATADIMENSION[3]),BLANK,overwrite=OVERWRITE,more=["NBPC",NBPC,"VARPERC",VARPERCENT[NBPC]*100,"METHOD","PCA"])
#     println("Matrix of PCs saved in in $(PATHTOSAVE)/Data/Yt_$(NBPC)PC_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")

#     if ismis == 1
#         mmean = Dataprep.addblank(M.mean,missingplaces2D[:,1],BLANK,(DATADIMENSION[1],DATADIMENSION[2]))
#         #projec            = Dataprep.addblank(PCA.proj(M),missingplaces2D[:,1:NBPC],BLANK,DATADIMENSION) 
#         mmean = reshape(mmean,(DATADIMENSION[1],DATADIMENSION[2]))
#     else
#         mmean = reshape(M.mean,(DATADIMENSION[1],DATADIMENSION[2]))
#     end
#     Dataprep.write_fits("$(FITSPATH)","mmean_$(NBPC)PC","$PATHTOSAVE/Data/",mmean,(DATADIMENSION[1],DATADIMENSION[2]),BLANK,overwrite=OVERWRITE,more=["NBPC",NBPC,"VARPERC",VARPERCENT[NBPC]*100,"METHOD","PCA"])
#     s = open("$(PATHTOSAVE)/Data/Yt_$(NBPC)PC.bin", "w+")
#     write(s,Yt)
#     close(s)
#     Ytpath = "$(PATHTOSAVE)/Data/Yt_$(NBPC)PC.bin"


#     proj = PCA.proj(M)
#     if ismis == 1
#         proj = Dataprep.addblank(proj,missingplaces2D[:,1:NBPC],BLANK,(DATADIMENSION[1],DATADIMENSION[2],NBPC))
#     end
#     proj = reshape(proj,(DATADIMENSION[1],DATADIMENSION[2],NBPC))
#     Dataprep.write_fits("$(FITSPATH)","$(SAVENAME)_projectionmatrix","$(PATHTOSAVE)/Data/",proj,(DATADIMENSION_NOMISSING[1],DATADIMENSION_NOMISSING[2],NBPC),BLANK,finished=true,overwrite=OVERWRITE)
#     xvector = range(1,NBPC)#[1:HIGHESTPC]

#     # Cleaning memory
#     cube = 0.0 
#     Yt   = 0.0
#     head = 0.0
#     GC.gc()
#     if ismis == 1
#         cubereconstructed = Dataprep.addblank(cubereconstructed,missingplaces2D,BLANK,DATADIMENSION)
#     end
#     cubereconstructed = reshape(cubereconstructed,DATADIMENSION)
#     #projec            = reshape(PCA.proj(M),DATADIMENSION)
#     #cubereconstructed = Dataprep.blank_equal(cubereconstructed,BLANK,0)
#     println("Saving Fits")
#     Dataprep.write_fits("$(FITSPATH)","RECONSTRUCTED_$(SAVENAME)_$(NBPC)PC","$PATHTOSAVE/Data/",cubereconstructed,DATADIMENSION,BLANK,overwrite=OVERWRITE,more=["NBPC",NBPC,"VARPERC",VARPERCENT[NBPC]*100,"METHOD","PCA"])
#     println("Data reconstructed from PCA saved in $(PATHTOSAVE)/Data/RECONSTRUCTED_$(SAVENAME)_$(NBPC)PC_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")


#     cubereconstructed = 0
#     GC.gc()
# end   #pca








# function structure_functions(FITSPATH,FILENAME,PATHTOSAVE,SAVENAME,ORDERSTXT,meth,BLANK,OVERWRITE;limi=0,limf=0)
#     ORDERS = [parse(Int, ss) for ss in split(ORDERSTXT,",")]


#     # Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
#     cvicube,DATADIMENSION,HEAD = Dataprep.read_fits_cvi("$(FITSPATH)" ; check=false)
#     haskey(HEAD,"THRESH") && (THRESHOLD = HEAD["THRESH"])
#     haskey(HEAD,"THRESH") || (THRESHOLD = 0)
#     if haskey(HEAD,"METHOD")==1 
#         if HEAD["METHOD"]=="PCA"
#             METH = HEAD["NBPC"]
#             METH = "$(METH)PC"
#             METHV = HEAD["NBPC"]
#         elseif HEAD["METHOD"]=="SWO"
#             METH = "SWO"
#             METHV = -1
#         end
#     else
#         METHV = 0
#         METH = "raw"
#     end

#     LAG = [parse(Int,ss) for ss in split(HEAD["LAG"][2:end-1],",")]

#     # Prepare directories where plots and data will be saved.
#     Dataprep.directory_prep(PATHTOSAVE)

#     cvicube = Dataprep.replace_nantomissing(cvicube)
#     cvicube = Dataprep.replace_blanktomissing(cvicube,BLANK)
#     cvicube = Dataprep.replace_blanktomissing(cvicube,0)

#     #println(cvicube)
#     if meth=="moninyaglom"
#         sct = Structure_functions.fct_sct(cvicube,LAG,ORDERS)  
#     elseif meth=="hily"
#         sct = Structure_functions.fct_sct_int(cvicube,LAG,ORDERS) 
#     else
#         error("The method given as an option when calling the structure_functions function is not correct. Please, use 'moninyaglom' or 'hily'")
#     end
#     nl = cat(0,LAG ; dims=1)
#     nsct = cat(reshape(nl,(1,size(nl)[1])),cat(ORDERS,sct ; dims=2);dims=1)

#     Dataprep.write_dat(nsct,"$(PATHTOSAVE)/Data/","$(SAVENAME)_Sp(l)_$(METH)", more=["METHOD $(METH) ; FILE : $(SAVENAME) ; Intensity threshold during CV computation : $(THRESHOLD). Each column is a lag, each row an order. First column give the orders, first row the lags. For information : p=$(ORDERS), and l=$(LAG) ;  " ], overwrite=OVERWRITE)


#     if limi==0 || limf==0
#         println(" ")
#         println("On which Lag to fit ? Give the indices in the array Lag of the var file (first indice=1). Size of the array : $(size(LAG)[1])")
#         println("First indice : ")
#         CANALINF   = parse(Int64,readline())
#         println("Second indice : ")
#         CANALSUP   = parse(Int64,readline())
#         CANALTOFIT = CANALINF:CANALSUP
#     else
#         CANALTOFIT = limi:limf
#     end    
#     zeta = Structure_functions.xhi_fct_p(ORDERS[:],sct[:,CANALTOFIT])


#     Dataprep.write_dat(cat([METHV 0 0 0],cat(zeta,ORDERS,dims=2),dims=1),"$(PATHTOSAVE)/Data/","$(SAVENAME)_stcfct_$(METH)", more=["METHOD $(METH) ; FILE : $(SAVENAME) ; Intensity threshold during CV computation : $(THRESHOLD). ROW are results for differents orders which are given at the last column. First column is the exponant, second column is the factor A : Sp(l)=A*S3(l)^B. Also, at first row and first column is the method used : <0 for SWO, 0 for raw cube, >0 for PCA. The value gives the number of PCs for PCA)" ], overwrite=OVERWRITE)


# end #function structure_function








# function swo(FITSPATH,FILENAME,PATHTOSAVE,SAVENAME,UNITVELOCITY,STEP,BLANK,NOISECAN,OVERWRITE)
#     println("Perform SWO")
#     # Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
#     cube,VELOCITYVECTOR,DATADIMENSION,VELOCITYINCREMENT,HEAD = Dataprep.read_fits_ppv("$(FITSPATH)",UNITVELOCITY ; check=false)
    
#     # Prepare directories where plots and data will be saved.
#     Dataprep.directory_prep(PATHTOSAVE)
    
#     # Replace any NaN value into a missing value and deleted them (can't do PCA on missing values with that package).
#     cube = Dataprep.replace_nantomissing(cube)
#     ismis = 0
#     if any(ismissing,cube) 
#         ismis = 1
#         cube,missingplaces1D,missingplaces2D  = Dataprep.pca_prep(cube,DATADIMENSION)
#         cube                                  = convert(Array{Float64},cube)
#         DATADIMENSION_NOMISSING               = Dataprep.read_dim(cube)
#     else
#         cube                                 = reshape(cube,DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
#         cube                                 = convert(Array{Float64},cube)
#         DATADIMENSION_NOMISSING              = (DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
#     end
    
#     # if meth=="swo"
#     maskinterv,mask,posimap = SWO.newswo(cube,DATADIMENSION_NOMISSING,VELOCITYVECTOR,NOISECAN)

#     if ismis == 1
#         maskinterv = Dataprep.addblank(maskinterv,missingplaces2D,BLANK,DATADIMENSION)
#         posimapinf = Dataprep.addblank(posimap[:,1],missingplaces2D[:,1],BLANK,(DATADIMENSION[1],DATADIMENSION[2]))
#         posimapsup = Dataprep.addblank(posimap[:,2],missingplaces2D[:,1],BLANK,(DATADIMENSION[1],DATADIMENSION[2]))
#         #maskintervpety = Dataprep.addblank(maskintervpety,missingplaces2D,BLANK,DATADIMENSION)
#     else 
#         posimapinf = posimap[:,1]
#         posimapsup = posimap[:,2]
#     end

#     maskinterv = reshape(maskinterv,DATADIMENSION)
#     maskinterv = Dataprep.blank_equal(maskinterv,BLANK,0)
#     posimap = Array{Float64}(undef, (DATADIMENSION[1],DATADIMENSION[2],2))
#     posimap .= BLANK
#     posimap[:,:,1] .= reshape(posimapinf,(DATADIMENSION[1],DATADIMENSION[2])) 
#     posimap[:,:,2] .= reshape(posimapsup,(DATADIMENSION[1],DATADIMENSION[2]))
#     posimap = Dataprep.replace_blanktomissing(posimap,BLANK)


#     # These next ~40 rows allow to treat spectra for which SWO didn't found an optimised window. Will average the positions that SWO found for the 5x5 pixels around the pixel without a window. Thus, can't work if multiples spectra doesn't have a window around them, but it is logical. 
#     # Would be better to include it in a dedicated function.
#     for px=1:size(maskinterv)[1]
#         for py=1:size(maskinterv)[2]
#             if maskinterv[px,py,3]!=0 && maskinterv[px,py,3]!=BLANK
#                 if px>STEP && px<DATADIMENSION[1]-STEP && py>STEP && py<DATADIMENSION[2]-STEP
#                     posi = floor(moment(collect(skipmissing(posimap[px-STEP:px+STEP,py-STEP:py+STEP,1])),1,0)) |> Int64
#                     posf = floor(moment(collect(skipmissing(posimap[px-STEP:px+STEP,py-STEP:py+STEP,2])),1,0))+1 |> Int64
#                 elseif px<STEP && py>STEP && py<DATADIMENSION[2]-STEP
#                     posi = floor(moment(collect(skipmissing(posimap[px:px+STEP,py-STEP:py+STEP,1])),1,0)) |> Int64
#                     posf = floor(moment(collect(skipmissing(posimap[px:px+STEP,py-STEP:py+STEP,2])),1,0))+1 |> Int64
#                 elseif py<STEP && px>STEP &&  px<DATADIMENSION[1]-STEP 
#                     posi = floor(moment(collect(skipmissing(posimap[px-STEP:px+STEP,py:py+STEP,1])),1,0)) |> Int64
#                     posf = floor(moment(collect(skipmissing(posimap[px-STEP:px+STEP,py:py+STEP,2])),1,0))+1 |> Int64
#                 elseif px>DATADIMENSION[1]-STEP && py>STEP && py<DATADIMENSION[2]-STEP
#                     posi = floor(moment(collect(skipmissing(posimap[px-STEP:px,py-STEP:py+STEP,1])),1,0)) |> Int64
#                     posf = floor(moment(collect(skipmissing(posimap[px-STEP:px,py-STEP:py+STEP,2])),1,0))+1 |> Int64
#                 elseif py>DATADIMENSION[2]-STEP && px>STEP &&  px<DATADIMENSION[1]-STEP
#                     posi = floor(moment(collect(skipmissing(posimap[px-STEP:px+STEP,py-STEP:py,1])),1,0)) |> Int64
#                     posf = floor(moment(collect(skipmissing(posimap[px-STEP:px+STEP,py-STEP:py,2])),1,0))+1 |> Int64
#                 elseif px==STEP && py>STEP && py<DATADIMENSION[2]-STEP
#                     posi = floor(moment(collect(skipmissing(posimap[px-STEP+1:px+STEP,py-STEP:py+STEP,1])),1,0)) |> Int64
#                     posf = floor(moment(collect(skipmissing(posimap[px-STEP+1:px+STEP,py-STEP:py+STEP,2])),1,0))+1 |> Int64
#                 elseif px==DATADIMENSION[1]-STEP && py>STEP && py<DATADIMENSION[2]-STEP
#                     posi = floor(moment(collect(skipmissing(posimap[px-STEP:px+STEP-1,py-STEP:py+STEP,1])),1,0)) |> Int64
#                     posf = floor(moment(collect(skipmissing(posimap[px-STEP:px+STEP-1,py-STEP:py+STEP,2])),1,0))+1 |> Int64
#                 elseif py==STEP && px>STEP && px<DATADIMENSION[1]-STEP
#                     posi = floor(moment(collect(skipmissing(posimap[px-STEP:px+STEP,py-STEP+1:py+STEP,1])),1,0)) |> Int64
#                     posf = floor(moment(collect(skipmissing(posimap[px-STEP:px+STEP,py-STEP+1:py+STEP,2])),1,0))+1 |> Int64
#                 elseif py==DATADIMENSION[2]-STEP && px>STEP && px<DATADIMENSION[1]-STEP
#                     posi = floor(moment(collect(skipmissing(posimap[px-STEP:px+STEP,py-STEP:py+STEP-1,1])),1,0)) |> Int64
#                     posf = floor(moment(collect(skipmissing(posimap[px-STEP:px+STEP,py-STEP:py+STEP-1,2])),1,0))+1 |> Int64
#                 elseif px==STEP && py==STEP
#                     posi = floor(moment(collect(skipmissing(posimap[px-STEP+1:px+STEP,py-STEP+1:py+STEP,1])),1,0)) |> Int64
#                     posf = floor(moment(collect(skipmissing(posimap[px-STEP+1:px+STEP,py-STEP+1:py+STEP,2])),1,0))+1 |> Int64
#                 elseif px==STEP && py==DATADIMENSION[2]-STEP
#                     posi = floor(moment(collect(skipmissing(posimap[px-STEP+1:px+STEP,py-STEP:py+STEP-1,1])),1,0)) |> Int64
#                     posf = floor(moment(collect(skipmissing(posimap[px-STEP+1:px+STEP,py-STEP:py+STEP-1,2])),1,0))+1 |> Int64
#                 elseif px==DATADIMENSION[2]-STEP && py==STEP
#                     posi = floor(moment(collect(skipmissing(posimap[px-STEP:px+STEP-1,py-STEP+1:py+STEP,1])),1,0)) |> Int64
#                     posf = floor(moment(collect(skipmissing(posimap[px-STEP:px+STEP-1,py-STEP+1:py+STEP,2])),1,0))+1 |> Int64
#                 elseif px==DATADIMENSION[2]-STEP && py==DATADIMENSION[2]-STEP
#                     posi = floor(moment(collect(skipmissing(posimap[px-STEP:px+STEP-1,py-STEP:py+STEP-1,1])),1,0)) |> Int64
#                     posf = floor(moment(collect(skipmissing(posimap[px-STEP:px+STEP-1,py-STEP:py+STEP-1,2])),1,0))+1 |> Int64
#                 elseif py>=DATADIMENSION[2]-STEP && px<=STEP 
#                     posi = floor(moment(collect(skipmissing(posimap[px:px+STEP,py-STEP+1:py,1])),1,0)) |> Int64
#                     posf = floor(moment(collect(skipmissing(posimap[px:px+STEP,py-STEP+1:py,2])),1,0))+1 |> Int64
#                 elseif py>=DATADIMENSION[2]-STEP && px>=DATADIMENSION[2]-STEP
#                     posi = floor(moment(collect(skipmissing(posimap[px-STEP:px,py-STEP+1:py,1])),1,0)) |> Int64
#                     posf = floor(moment(collect(skipmissing(posimap[px-STEP:px,py-STEP+1:py,2])),1,0))+1 |> Int64
#                 elseif py<=STEP && px<=STEP 
#                     posi = floor(moment(collect(skipmissing(posimap[px:px+STEP,py:py+STEP,1])),1,0)) |> Int64
#                     posf = floor(moment(collect(skipmissing(posimap[px:px+STEP,py:py+STEP,2])),1,0))+1 |> Int64
#                 elseif py<=STEP && px>=DATADIMENSION[2]-STEP
#                     posi = floor(moment(collect(skipmissing(posimap[px-STEP:px,py:py+STEP,1])),1,0)) |> Int64
#                     posf = floor(moment(collect(skipmissing(posimap[px-STEP:px,py:py+STEP,2])),1,0))+1 |> Int64
#                 else 
#                     posi = 2
#                     posf = DATADIMENSION[3]-1
#                 end
                
#                 if posf<=0
#                     posf = DATADIMENSION[3]-1
#                 end
                
#                 if posi<=0
#                     posi = 2
#                 end
#                 maskinterv[px,py,1:posi] .= 0
#                 maskinterv[px,py,posf:end] .=0
#             end
#         end
#     end

#     Dataprep.write_fits("$(FITSPATH)","RECONSTRUCTED_$(SAVENAME)_SWO","$PATHTOSAVE/Data/",maskinterv,DATADIMENSION,BLANK,overwrite=OVERWRITE,more=["METHOD","SWO"])
#     println("Data reconstructed from SWO method saved in $(PATHTOSAVE)/Data/RECONSTRUCTED_$(SAVENAME)_SWO_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")


# end















end # module UnveilGUI
