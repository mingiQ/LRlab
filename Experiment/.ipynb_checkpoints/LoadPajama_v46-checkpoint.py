# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 21:20:58 2018

@author: Tailung　

version: 4.6
   review 7/30/2024
   - Rewrite "Gate_density" function due to Pandas drops "append" function
   - Create "HallAnalysis" and "SdHAnalysis" class
   - Change "Hall" function output format to dictionary
   - Add Gauss function 
   - Add "LoadDats" function: load multiples data and output a single dictionary

version: 4.5
  review: 6/14/2024
  Rewrite "LoadDatsX" because Pandas drops "Dataframe.append" function

version: 4.4
  review: 12/2/2021
  Modify "Sym_Correct" function


version: 4.3
  review: 11/07/2021
  Add line style option in "plotxy" function

verson: 4.2
  review : 11/01/2021
  modify "Load_QCodesDB" function

version：　４．１
　—　Ｍｏｄｉｆｙ　＂twin_xy＂　ｔｏ　ｓｕｐｐｏｒｔ　ｍｕｌｔｉｐｌｅ　ｃｕｒｖｅｓ　ｏｎ　ｅａｃｈ　ａｘｉｓ


version: 4.0

Review at 8/6/2021
  -- Add Load PPMS data format and Output Pandas data format

verstion: 3.2

Review at 10/27/2019
  -- Add Data cursor class: Cursor and SnaptoCursor

Review at 7/18/2019
  -- Add "Load_data_QCodesDB" function

Review at 6/11/2019
  -- Add 1D-WAL function in dirty limit :WAL_1D_Dirty(B,l_a,l_so):
  -- Add symetric correction function : Sym_Correct
  
Review at 6/05/2019
  -- Add twin_xy function -- twin_xy
  -- Add extract Hall density n(Vg) from Rxy(B,Vg) 2D data-- Gate_density
"""

########################
##Package
########################
import pandas as pd
from scipy import stats
import scipy.constants as sc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib import cm
from matplotlib.widgets import Cursor
from matplotlib import rcParams

#rcParams['font.family'] = 'serif'
#rcParams['font.serif'] = ['Times New Roman']

rcParams['xtick.top'] = rcParams['ytick.right'] = True

from matplotlib.ticker import AutoMinorLocator,FormatStrFormatter,LinearLocator
#from mpldatacursor import datacursor
#########################

#############################################
## Print in Notebook with markdown format  ##
#############################################
from IPython.display import Markdown, display
def printmd(string:str, color=None):
    """
    Display markdown string
    color option: 'green', 'blue', 'red', 'yellow'
    """
    colorstr = "<span style='color:{}'>{}</span>".format(color, string)
    display(Markdown(colorstr))


##########################################
## Load data from PPMS Dataformat       ##
##########################################

def LoadPPMS(filename:str='', # PPMS Data file name
             Columns=[], # Column name list which we would like to extract
             Key_col = '', # Subset data with Key_col column
             Subset=(), # Subset base on "Key_col" value, extract data between (Value min,max) or () load all
             Cols_name:list=[] # rename column names or [] do not change
             ):
    """
    Load PPMS data format and output a pandas data format\n
    Skip first 30 rows\n
    ===========\n
    Parameters:\n
    filename: PPMS data file name\n
    Columns : Column name list which we would like to extract\n
    Key_col : Specific Column name which subset is extracted
    Subset duplet: Subset base on "Key_col" value, extract data between (Value min,max) or () load all\n
    Cols_name : rename column names or [] do not change\n
    \n
    Return Pandas dataframe \n
    """
    if Columns != []:
        df = pd.read_csv(filename,skiprows=30).filter(items= Columns, axis=1)
        df = df.dropna()
        if Key_col != '':
            df = df.loc[(df[Key_col] > Subset[0]) & (df[Key_col] < Subset[1])]
            
        if Cols_name != []:
            df.columns=Cols_name
            
        df = df.reset_index(drop=True)
        
    else:
        df = pd.read_csv(filename,skiprows=30)
        
        
    return df


##########################################
## Load data from QCodes Dataset        ##
##########################################

def Load_QCodesDB(Run_id:int, # exp. run ID in database
                  Cols:list=[]):  # rename column names
    """
    Load single experiment run from QCodes Dataset into a pandas data frame\n
    ===========\n
    How to load Qcodes Database
       from qcodes import initialise_or_create_database_at
       db_file_path = os.path.join(os.getcwd(), 'M10-23-20&M10-15-20.db')
       initialise_or_create_database_at(db_file_path)
    How to load Exp. from Database
       from qcodes.dataset.experiment_container import load_experiment
       exp = load_experiment(ExpID)
    How to check Runs
       from qcodes.dataset.experiment_container import experiments
       experiments()
    ===========\n
    Parameters:\n
    Run_id : experiment run ID in dataset. Run ID can be checked by "experiments()" command \n
    Cols : array of column names.\n
    \n
    Return Pandas dataframe \n
    """
    from qcodes.dataset.data_set import load_by_id
    #tmp=load_by_id(Run_id).get_data_as_pandas_dataframe()
    tmp=load_by_id(Run_id).to_pandas_dataframe_dict()
    i=0
    for s in tmp.values(): # replace column name
        i+=1
        if i==1:
            df=s
        else:
            df = pd.concat([df,s], axis=1)
        
    df=df.reset_index()
    if Cols != []:
        df.columns=Cols
    
    del i,s,tmp
    return df



##########################################
## Load data from QCodes program        ##
##########################################
def Load_QCodes(filename:str,
                Cols:list=[]
               ):
    r"""
    Load Qcode data file (not Database) into a pandas data frame\n
    Parameters:\n
    ===========\n
    Filename: Qcodes data file\n
    Cols : array of column names. Default: the 2nd row in data\n
    \n
    Return Pandas dataframe \n\n
    Example:\n
    ===========\n
    cols_name=['Vg (V)',r'Rxx (k$\Omega$)',r'Rxy ($\Omega$)','Ig (pA)']\n
    df002=Load_QCodes('data.dat',cols_name)
    """
    df=pd.read_csv(filename,sep='\t',skiprows=[0,2],header=0)
    if Cols != []:
        df.columns=Cols

    return df

#################################################
## Load Jimmy's pajama data  into X,Y,Z column ##
#################################################
def Load_XYZ(filename,Cols):
    """
    Load Jimmy's pajama data into one matrix in pandas data frame\n
    Parameters:\n
    ===========\n
    Filename: pajama data file\n
    Col_labels : array of column names, Example ['V','I','R']\n
    1st Column: major variable (x) ; 2nd Column : secondary variable(y); 3rd Column : Data\n
    \n
    Return Pandas dataframe
    """
    df=pd.read_csv(filename,sep='\n',names=['A'])
    foo = lambda x: pd.Series([i for i in x.split('\t')[:-1]])
    data = df['A'].apply(foo) # row line slice
        
    boo = lambda x: pd.Series([i for i in x.split('|')[:len(Cols)]])
    # Extract index list
    Out = pd.DataFrame(columns=Cols)
    
    for i in range(len(data.index)):
        rev2 = data.iloc[i,:].dropna(axis=0,how='any').apply(boo)
        rev2.columns=Cols
        Out=Out.append(rev2,ignore_index=True)
    
    for i in range(len(rev2.columns)):
        Out.iloc[:,i]=pd.to_numeric(Out.iloc[:,i])    
    
    return Out


##########################################
## Load Jimmy's pajama data into Matrix ##
##########################################
def Load_Matrix(filename,X_label='',Y_label='',dataName='',z=2,x=0,y=1):
    """
    Load Jimmy's pajama data into one matrix in pandas data frame\n
    Parameters:\n
    ===========\n
    Filename: pajama data file\n
    z : Column index for z (default : 2)\n
    x : Column index for x (default : 0)\n
    y : Column index for y (default : 1)\n\n
    X_label : name of X\n
    Y_label : name of y\n
    dataName : name of dataframe\n
    \n
    Return Pandas dataframe
    """
    df=pd.read_csv(filename,sep='\n',names=['A'])
    foo = lambda x: pd.Series([i for i in reversed(x.split('\t'))])
    data = df['A'].apply(foo)
    data.drop(columns=[0],inplace=True)
    
    boo = lambda x: pd.Series([i for i in (x.split('|'))])
    # Extract index list
    rev2 = data.iloc[:,0].apply(boo)
    idx=pd.to_numeric(rev2.iloc[:,y])  #Extract index list ; y
    
    cols=[]
    for i in range(len(data.columns)):
        rev2 = data.iloc[:,i].apply(boo)
        rev2.drop(columns=[3],inplace=True)
        cols=np.append(cols,float(rev2.iloc[0,x])) #Extract column list ; x
        data.iloc[:,i]=pd.to_numeric(rev2.iloc[:,z]) # Z-data
        
    data.columns=cols
    data.index=idx
    
    if X_label != '':
        data.columns.name = X_label
    
    if Y_label != '':
        data.index.name = Y_label
    
    if dataName != '':
        data.name = dataName
    
    return data


############################################################
#Load LRScan in Orgin format: *.DAT
#First 2 rows in file are Title and unit individually
###########################################################    
def LoadDat(filename:str):# file
    """
    Load LR data into one pandas data frame\n
    Parameters:\n
    ===========\n
    Filename of LR data file (*.dat)\n
    \n
    Return Pandas dataframe with title and unit in Column index\n
    Each column index has the following format : "title (unit)"
    """
# af=1: attach filename in front of title; af!=0 don't attach filename
# Assemble title of data in the format of "title (unit)" from first 2 rows
    f = open(filename, 'r')
    Title= f.readline()[:-2].split(",")
    Unit= f.readline()[:-2].split(",")
    f.close
    for i in range(len(Title)):
        Title[i]=" (".join((Title[i], Unit[i]))+")"  ##Construct Title_Unit as column name
#load data into pandas data frame
    return pd.read_csv(filename,sep=',',skiprows=2,names=Title)


def LoadDats(fileformat:str, # fileformat
             file_index_tuple:tuple=(), # file index tuple
             ):
    """
    Load Multiple LR data files
    Output dictionary format
    \n
    Parameters:\n
    ===========\n
    Fileformat: (Example) '%03d.DAT', %03d is the zero fills 3 digit integer variable\n
    file_index_tuple: file index tuple\n
    """
    Result_dict={}
    for file_idx in file_index_tuple:
        #Construct column title and units 
        f = open((fileformat %(file_idx)), 'r')
        Title= f.readline()[:-2].split(",")
        Unit= f.readline()[:-2].split(",")
        f.close
        for i in range(len(Title)):
            Title[i]=" (".join((Title[i], Unit[i]))+")"
        # Save data into output dictionary
        Result_dict["df%03d"%(file_idx)]= pd.read_csv((fileformat %(file_idx)),sep=',',skiprows=2,names=Title)

    return Result_dict



#### load multiple data file and merge into one single pandas dataframe

def LoadDats1(fileformat:str, # fileformat
              f_i:int,f_f:int, # Initial and final file index number
              B_i:float,B_step:float, # Initial B-field value and step 
              B_colName:str='B (T)'): # Colunm name for the 2nd independent parameter
    """
    Load Multiple LR data files and merge into one single pandas data frame\n
    Add the 2nd independent parameter (ex. B-field) as the 1st column\n
    \n
    Parameters:\n
    ===========\n
    Fileformat: (Example) '%03d.DAT', %03d is the zero fills 3 digit integer variable\n
    f_i and f_f: Initial and final file index number, files index should be consecutive integer number.\n
    B_i : Initial value for the 2nd independent parameter (ex. B-field, T)\n
    B_step:  the 2nd independent parameter step (delta B-field)\n
    B_colName : Column Name for  the 2nd independent parameter (B-field, default 'B (T))\n
    """
    #Extract column title and units from first file
    f = open((fileformat %(f_i)), 'r')
    Title= f.readline()[:-2].split(",")
    Unit= f.readline()[:-2].split(",")
    f.close
    #Construct Title
    # Unit as column name
    for i in range(len(Title)):
        Title[i]=" (".join((Title[i], Unit[i]))+")"

    #Construct the column for output dataframe
    df = pd.read_csv((fileformat %(f_i)),sep=',',skiprows=2,names=Title)
    B_df = pd.DataFrame({B_colName:np.full(df.shape[0], B_i)})
    
    for j in range((f_i+1),(f_f+1)):
        Tmp = pd.read_csv((fileformat %(j)),sep=',',skiprows=2,names=Title)
        df = pd.concat([df,Tmp])
        B_s = B_i+(j-f_i)*B_step
        B_df = pd.concat([B_df,
                          pd.DataFrame({B_colName:np.full(Tmp.shape[0], B_s)})
                         ])

    df.insert(0,B_colName,B_df)
    df.reset_index(drop=True, inplace=True)
    return df


def LoadDats2(fileformat:str, # fileformat
              f_i:int, # Initial file index number
              f_f:int, # Final file index number
              B_tuple:tuple=(), #  The corresponding value of the 2nd independent parameter 
              B_colName:str='B (T)'):# Colunm name for the 2nd independent parameter
    """
    Load Multiple LR data files and merge into one single pandas data frame\n
    Add the 2nd independent parameter (ex. B-field) as the 1st column\n
    \n
    Parameters:\n
    ===========\n
    Fileformat: (Example) '%03d.DAT', %03d is the zero fills 3 digit integer variable\n
    f_i and f_f: Initial and final file index number, files index should be consecutive integer number.\n
    B_tuple : The corresponding value of the 2nd independent parameter (ex. B-field, T)\n
    B_colName : Column Name for  the 2nd independent parameter (B-field, default 'B (T))\n
    """
    #Extract column title and units from first file
    f = open((fileformat %(f_i)), 'r')
    Title= f.readline()[:-2].split(",")
    Unit= f.readline()[:-2].split(",")
    f.close
    #Construct Title
    # Unit as column name
    for i in range(len(Title)):
        Title[i]=" (".join((Title[i], Unit[i]))+")"

    #Construct the column for output dataframe
    df = pd.read_csv((fileformat %(f_i)),sep=',',skiprows=2,names=Title)
    B_df = pd.DataFrame({B_colName:np.full(df.shape[0], B_tuple[0])})
    i=1
    for j in range((f_i+1),(f_f+1)):
        Tmp = pd.read_csv((fileformat %(j)),sep=',',skiprows=2,names=Title)
        df = pd.concat([df,Tmp])
        B_df = pd.concat([B_df,
                          pd.DataFrame({B_colName:np.full(Tmp.shape[0], B_tuple[i])})
                         ])
        i=i+1

    df.insert(0,B_colName,B_df)
    df.reset_index(drop=True, inplace=True)
    return df

def LoadDats3(fileformat:str, # fileformat
              file_index_tuple:tuple=(), # file index tuple
              B_tuple:tuple=(), #  The corresponding value of the 2nd independent parameter 
              B_colName:str='B (T)' # Colunm name for the 2nd independent parameter
             ):
    """
    Load Multiple LR data files and merge into one single pandas data frame\n
    Add the 2nd independent parameter (ex. B-field) as the 1st column\n
    \n
    Parameters:\n
    ===========\n
    Fileformat: (Example) '%03d.DAT', %03d is the zero fills 3 digit integer variable\n
    file_index_tuple: file index tuple\n
    B_tuple : The 2nd independent parameter array (ex. B-field, T)\n
    B_colName : Column Name for  the 2nd independent parameter (B-field, default 'B (T))\n
    """
    #Extract column title and units from first file
    f = open((fileformat %(file_index_tuple[0])), 'r')
    Title= f.readline()[:-2].split(",")
    Unit= f.readline()[:-2].split(",")
    f.close
    #Construct Title
    # Unit as column name
    for i in range(len(Title)):
        Title[i]=" (".join((Title[i], Unit[i]))+")"

    #Construct the dataframe for 1st file 
    df = pd.read_csv((fileformat %(file_index_tuple[0])),sep=',',skiprows=2,names=Title)
    # create the corresponding constant B-field array for 1st file
    B_df = pd.DataFrame({B_colName:np.full(df.shape[0], B_tuple[0])})

    for (file_idx, B_field) in zip(file_index_tuple[1:], B_tuple[1:]): # loop from the 2nd file
        Tmp = pd.read_csv((fileformat %(file_idx)),sep=',',skiprows=2,names=Title)
        df = pd.concat([df,Tmp]) 
        # corresponding B-field array 
        B_df = pd.concat([B_df,
                          pd.DataFrame({B_colName:np.full(Tmp.shape[0], B_field)})
                         ])

    df.insert(0,B_colName,B_df)
    del B_df
    df.reset_index(drop=True, inplace=True)
    return df




def MergeDataFrame(Data_dic:dict, # Dictionary contain identical pandas dataframes 
                   Select_Col:list, # select columns which will merge
                   IndParam:tuple=(), #  The corresponding value of the 2nd independent parameter 
                   IndParam_ColName:str='' # Colunm name for the 2nd independent parameter
                  ):
    """
    Merge pandas dataframe from dictionary with selected columns\n
    Then, output one single dataframe\n
    Add the 2nd independent parameter (ex. B or Vg) as the 1st column of output dataframe\n
    
    Parameters:\n
    ===========\n
    Data_dic: Input dictionary which contain multiple and identical pandas dataframes\n
    Select_Col: string list of columns which will merge\n
    IndParam: tuple, the corresponding value of the 2nd independent parameter\n 
    IndParam_ColName: string, column name of the 2nd independent parameter\n 

    Output:\n
    ===========\n
    Dataframe
    """
    i = 0
    for _,Ind_param in zip(Data_dic.values(),IndParam):
        if i ==0:
            Out_df = _[Select_Col]
            IndParam_df = pd.DataFrame({IndParam_ColName:np.full(_.shape[0], Ind_param)})
        else:
            Out_df = pd.concat([Out_df,_[Select_Col]])
            IndParam_df = pd.concat([IndParam_df,
                                     pd.DataFrame({IndParam_ColName:np.full(_.shape[0], Ind_param)})
                                    ])
        i+=1

    Out_df.insert(0,IndParam_ColName,IndParam_df)
    Out_df.reset_index(drop=True, inplace=True)
    return Out_df
            

###############################################################################
##  Ploting function
###############################################################################


#####################################
## Plot X,Y,Z color map
#####################################
def XYZColorMap(DataSet:pd.DataFrame,
                z:int=2,x:int=0,y:int=1,
                Zmin=None,Zmax=None,Zstp_nmb=100,
                title:str='',
                zlabel:str='',xlabel:str='',ylabel:str='',
                N=200, #contour N automatically-chosen levels.
                color=cm.viridis):
    """
    X,Y,Z Color coutour plot from Pandas dataframe\n
    Parameters:\n
    ===========\n
    DataSet: Pandas dataframe\n
    z : Column index for z (default : 2)\n
    x : Column index for x (default : 0)\n
    y : Column index for y (default : 1)\n
    Zmin (default: None), Zmax (default: None), Zstp_nmb (default: 100)\n
    title: title of figure\n
    zlabel,xlabel,ylabel (default is read from column name)\n
    N:contour N automatically-chosen levels.(default: 200)\n
    color:(default) cm.viridis\n
      'Miscellaneous', [\n
        'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',\n
        'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg', 'hsv',\n
        'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar']\n
      'Diverging', [\n
        'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',\n
        'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic'])\n\n\n
    Return:\n
    ===========\n
    Return fig,ax
    
    Example:\n
    ===========\n
    from matplotlib import cm\n
    fig,ax=tw.XYZColorMap(df006,z=2,x=1,y=0,\n
                          title="Run ID#006; Exp.:4K_001 (M06-19-19.2_E-A1HB3)",\n
                          color=cm.RdBu)
    """
    if Zmax == None:
        Zmax = DataSet.iloc[:,z].values.max()
    
    if Zmin == None:
        Zmin = DataSet.iloc[:,z].values.min()
        
    levels = np.linspace(Zmin,Zmax, num=Zstp_nmb)
    
    fig=plt.figure()
    ax=fig.add_subplot(111)
    cax=ax.tricontourf(DataSet.iloc[:,x],DataSet.iloc[:,y],DataSet.iloc[:,z]
                    ,N #contour N automatically-chosen levels.
                    ,cmap=color
                    ,vmax=Zmax, vmin=Zmin
                    ,extend = 'both'
                    ,levels=levels)
    
    # Axis label    
    if xlabel == '':
        xlabel = DataSet.columns[x]
    if ylabel == '':
        ylabel = DataSet.columns[y]    
    if zlabel == '':
        zlabel = DataSet.columns[z]
    ax.set_xlabel(xlabel, fontsize = 16)
    ax.set_ylabel(ylabel, fontsize = 16)
    
    #Color bar
    cbar=fig.colorbar(cax)
    cbar.set_label(zlabel,rotation=270 ,labelpad=14 ,fontsize = 14)
    cbar.ax.tick_params(axis='y', direction='in')
    
    #Axis tick
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    
    #plt title
    if title != '' :
        ax.set_title(title, style='italic')
    
    #Data Cursor
    #dc1 = datacursor(ax, xytext=(15, -15), bbox=None)
   
    # set useblit = True on gtkagg for enhanced performance
    #Cursor(ax,useblit=True, color='black', linewidth=0.5,linestyle='--')
    
    #fig.savefig('plot.png', dpi=600,transparent=True)

    plt.show()
    return fig,ax


#################################
## Plot Matrix into color map  ##
#################################
#DataSet : Panda Data frame
def MatrixColorMap (DataSet:pd.DataFrame,
                    title:str='',
                    zlable:str='',xlable:str='',ylable:str='',
                    Zmax=None,Zmin=None,
                    color=cm.coolwarm):
    """
    Color coutour plot of DataSet\n
    Parameters:\n
    ===========\n
    DataSet: Pandas dataframe\n
    Xlabel:\n
    Ylabel:\n
    Zlabel:\n
    \n
    Zmax : Max color scale in plot. 'None': Autoscale\n
    Zmin : Min color scale in plot. 'None': Autoscale\n
    \n
    Title:\n
    Color: (default) cm.coolwarm\n
      'Miscellaneous', [\n
        'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',\n
        'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg', 'hsv',\n
        'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar']\n
      'Diverging', [\n
        'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',\n
        'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic'])
    """
    fig =plt.figure()
    ax=fig.add_subplot(111)
    cax=ax.imshow(DataSet, interpolation='bilinear'
            ,cmap=color,aspect='auto'
            ,origin='high'
            ,extent=[DataSet.columns[0]   # Xmin
                     ,DataSet.columns[DataSet.shape[1]-1]  #Xmax
                     ,DataSet.index[0]    
                     ,DataSet.index[DataSet.shape[0]-1]]#[Xmin,Xmax,Ymin,Ymax]
            ,vmax= Zmax ,vmin=Zmin    # Z range
              )
    
    # Axis label    
    if xlable == '':
        xlable = DataSet.columns.name
    if ylable == '':
        ylable = DataSet.index.name
    
    if xlable != None:
        ax.set_xlabel(xlable, fontsize = 16)
    if ylable != None:
        ax.set_ylabel(ylable, fontsize = 16)
        
    if title != '' :
        ax.set_title(title, style='italic')
    #Cursor
    datacursor(ax, bbox=dict(fc='white'),
                 arrowprops=dict(arrowstyle='simple', fc='white', alpha=0.5))
    
    #Color bar
    if zlable == '' :
        zlable = DataSet.name
        
    cbar=fig.colorbar(cax)
    cbar.set_label(zlable
                 ,rotation=270 ,labelpad=12 ,fontsize = 14)
    cbar.ax.tick_params(axis='y', direction='in')
    
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    
    plt.show()



############################
## Plot X,Y,Z 3D plot     ##
############################
def XYZ3DPlot(DataSet:pd.DataFrame,
              z:int=2,x:int=0,y:int=1,
              title:str='',
              zlabel:str='',xlabel:str='',ylabel:str='',
              color=cm.coolwarm):
    """
    X,Y,Z Color 3D plot\n
    Parameters:\n
    ===========\n
    DataSet: Pandas dataframe\n
    z : Column index for z (default : 2)\n
    x : Column index for x (default : 0)\n
    y : Column index for y (default : 1)\n
    title:\n
    xlabel,ylabel,zlabel:\n
    color
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #ax = fig.gca(projection='3d')
    surf = ax.plot_trisurf(DataSet.iloc[:,0],DataSet.iloc[:,1],DataSet.iloc[:,2]
                    ,cmap=color)
    ax.zaxis.set_major_locator(LinearLocator(5))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.01f'))
    
    # Axis label    
    if xlabel == '':
        xlabel = DataSet.columns[x]
    if ylabel == '':
        ylabel = DataSet.columns[y]    
    if zlabel == '':
        zlabel = DataSet.columns[z]
    ax.set_xlabel(xlabel, fontsize = 14)
    ax.set_ylabel(ylabel, fontsize = 14)
    
    #Color bar
    cb=fig.colorbar(surf, spacing='proportional', aspect=15)
    cb.set_label(zlabel, rotation=270,labelpad=9)
    cb.ax.tick_params(axis='y', direction='in')
    
    #plt title
    if title != '' :
        ax.set_title(title, style='italic')

    fig.tight_layout()
    #plt.show()
    
    return fig,ax

#######################################
##  Twin lines plots in single panel ##
#######################################
    
def twin_xy (DataSet_L:tuple=(),DataSet_R:tuple=(), # DataSet_L : Pandas dataframe for left axis ; 
             X_L_index:tuple=(),Y_L_index:tuple=(),
             X_R_index:tuple=(),Y_R_index:tuple=(),
             Xlim:tuple=(),Y_left_lim:tuple=(),Y_right_lim:tuple=(),
             lables_L:tuple=(),lables_R:tuple=(),
             Xlable:str='',Y_L_lable:str='',Y_R_lable:str='',
             lo:int=0,
             title:str='',
             legend=False,
            figsize=(8,5)):
    """
    Lines plots of pair columns from different DataSet\n
    Parameters:\n
    ===========\n
    DataSet_L : tuple of Pandas dataframe for left axis. ex: (df01, df02)\n
    DataSet_R : tuple of Pandas dataframe for right axis. ex: (df03, df04)\n
    X_index: tuple of columns index of DataSet (x_left,x2_right)\n
    Y_index: tuple of columns index of DataSet (y_left,y_rgith)\n
    Xlim: tuple (x_min, x_max), default: autoscale\n
    Y_left_lim: tuple (y_min, y_max), default: autoscale\n
    Y_right_lim: tuple (y_min, y_max), default: autoscale\n
    lables: tupe of strings ('X_label','Y_left_label','Y_right_label')
    lo: location of legend (0-10, 0:'best')\n
    title: title of figure. (default) Empty string\n
    Example:\n
    ===========\n
    fig,ax_l,ax_r = twin_xy((df,df),\n
                            X_index= (0,0), # x_left,x_right\n
                            Y_index= (1,2)) # y_left,y_right\n
                            title= 'Run ID#004; Exp.:4K_001 (M06-19-19.2_E-A1HB3)')
    """    
    fig, ax = plt.subplots(1,1, figsize=figsize)
    #ls_list=['-','--','-.',':','.','loosely dashdotted']   #
    ls_list = ['-', '--', '-.', ':',  'solid', 'dashed', 'dashdot', 'dotted']
    ###
    color_L = 'tab:blue'
    i=0
    for j in range(len(DataSet_L)):
        if lables_L==():
            b=DataSet_L[j].columns[Y_L_index[j]]
        else:
            b=lables_L[j]
        
        ax.plot(DataSet_L[j].iloc[:,X_L_index[j]],
                DataSet_L[j].iloc[:,Y_L_index[j]],
                label=b,marker='.', ls=ls_list[i],color=color_L)
        i += 1

            
    # Axis label    
    if Xlable == '':
        Xlable = DataSet_L[0].columns[X_L_index[0]]
        
    if Y_L_lable == '':
        Y_L_lable = DataSet_L[0].columns[Y_L_index[0]]
    
    ax.set_xlabel(Xlable, fontsize = 14,fontweight='bold')
    ax.set_ylabel(Y_L_lable, fontsize = 14,fontweight='bold',color= color_L)
    
    #plot x- or y-lim
    if Xlim != ():
        ax.set_xlim([Xlim[0],Xlim[1]])
    if Y_left_lim != ():
        ax.set_ylim([Y_left_lim[0],Y_left_lim[1]])    
    
    #plt title
    if title != '' :
        ax.set_title(title, style='italic',fontweight='bold', wrap=True)
        
    #ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    
    ax.tick_params(which='both', direction='in')
    
    
    ax.spines['left'].set_color(color_L)
    ax.tick_params(axis='y', colors=color_L)
    ax.spines['right'].set_visible(False)

    del b,i
    
    #Plot right axis
    ax2= ax.twinx() # instantiate a second axes that shares the same x-axis
    color_R = 'tab:red'
    i=0
    for j in range(len(DataSet_R)):
        if lables_R==():
            b=DataSet_R[j].columns[Y_R_index[j]]
        else:
            b=lables_R[j]
        
        ax2.plot(DataSet_R[j].iloc[:,X_R_index[j]],
                 DataSet_R[j].iloc[:,Y_R_index[j]],
                 label=b,marker='.',ls=ls_list[i],color=color_R)
        i += 1
    
    
    if Y_R_lable == '':
        Y_R_lable = DataSet_R[0].columns[Y_R_index[0]]
    
    ax2.set_ylabel(Y_R_lable, fontsize = 14,rotation=270,
                   labelpad=18,fontweight='bold',color= color_R)
    
    if Y_right_lim != ():
        ax2.set_ylim([Y_right_lim[0],Y_right_lim[1]])
    
    ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
    
    ax2.tick_params(which='both', direction='in')
    
    ax2.spines['right'].set_color(color_R)
    ax2.tick_params(axis='y', colors=color_R)
    ax2.spines['left'].set_visible(False)

    if legend:
        ax.legend(loc=lo,frameon=False) 
        ax2.legend(loc=lo, frameon=False)
    
    fig.tight_layout()
    #plt.show()
    
    
    return fig,ax,ax2
    
    
    

##################################
##  Lines plot of pair columns  ##
##################################

def plotxy (DataSet:tuple=(),
            X_index:tuple=(),Y_index:tuple=(),
            LnLegend:tuple=(),
            colors:tuple=(),
            Xlim:tuple=(),Ylim:tuple=(),
            Xlable:str='',Ylable:str='',
            ls:tuple=(), # line style
            lo:int=0,
            title:str='',
            figsize=(8,5)
           ):
    r"""
    Lines plots of pair columns from different DataSet\n
    Parameters:\n
    ===========\n
    DataSet: tuple of dataframe\n
    X_index: tuple of columns index of DataSet (x1,x2,x3,...)\n
    Y_index: tuple of columns index of DataSet (y1,y2,y3,...)\n
    LnLegend: Lengend for each curve ('plot1','plot2',...)\n
    Xlim: tuple (x_min, x_max), default: autoscale\n
    Ylim: tuple (y_min, y_max), default: autoscale\n
    Xlabel: (default) Column name from the 1st curve\n
    Ylabel: (default) Column name from the 1st curve\n
    ls : tuple of line style for each trace\n
        Line style : ['-','--','-.',':','.','loosely dashdotted']\n
    Legend: (default) None\n
    lo: location of legend (0-10, 0:'best')\n
    title: title of figure. (default) Empty string\n
    colors: tuple of color for each trace\n
      ‘b’	blue\n
      ‘g’	green\n
      ‘r’	red\n
      ‘c’	cyan\n
      ‘m’	magenta\n
      ‘y’	yellow\n
      ‘k’	black\ 
      ‘w’	white\n\n
     default: '#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf'
    Example:\n
    ===========\n
    """    
    fig, ax = plt.subplots(1,1, figsize=figsize)
    for j in range(len(X_index)):
        if LnLegend==():
            b=DataSet[j].columns[Y_index[j]]
        else:
            b=LnLegend[j]
            
        if ls == (): # line style
            set_ls = 'solid' # solid line as default setting
        else:
            set_ls = ls[j]
        
        if colors==():
            ax.plot(DataSet[j].iloc[:,X_index[j]],
                    DataSet[j].iloc[:,Y_index[j]],
                    label=b,marker='.',ls=set_ls)
        else:
            ax.plot(DataSet[j].iloc[:,X_index[j]],
                    DataSet[j].iloc[:,Y_index[j]],
                    label=b,color=colors[j],marker='.',ls=set_ls)
            
    # Axis label    
    if Xlable == '':
        Xlable = DataSet[0].columns[X_index[0]]
        
    if Ylable == '':
        Ylable = DataSet[0].columns[Y_index[0]]
    
    ax.set_xlabel(Xlable, fontsize = 15,fontweight='bold')
    ax.set_ylabel(Ylable, fontsize = 15,fontweight='bold')
    
    #plot x- or y-lim
    if Xlim != ():
        ax.set_xlim([Xlim[0],Xlim[1]])
    if Ylim != ():
        ax.set_ylim([Ylim[0],Ylim[1]])    
    
    #plt title
    if title != '' :
        ax.set_title(title, style='italic') 
        
    #ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    
    ax.tick_params(which='both', direction='in')
    ax.legend(loc=lo,frameon=False) 
    
    #ax.axhline(0, color='black', lw=0.5)

    del b
    #plt.show()
    return fig,ax



##################################
##  Lines plot of pair columns  ##
##################################

def ScatterXY (DataSet:tuple=(),
               X_index:tuple=(),Y_index:tuple=(),
               LnLegend:tuple=(),
               Xlim:tuple=(),Ylim:tuple=(),
               Xlable:str='',Ylable:str='',
               lo:int=0,
               title:str=''
              ):
    """
    Lines plots of pair columns from different DataSet\n
    Parameters:\n
    ===========\n
    DataSet: tuple of dataframe\n
    X_index: tuple of columns index of DataSet (x1,x2,x3,...)\n
    Y_index: tuple of columns index of DataSet (y1,y2,y3,...)\n
    LnLegend: Lengend for each curve ('plot1','plot2',...)\n
    Xlim: tuple (x_min, x_max), default: autoscale\n
    Ylim: tuple (y_min, y_max), default: autoscale\n
    Xlabel: (default) Column name from the 1st curve\n
    Ylabel: (default) Column name from the 1st curve\n
    Legend: (default) None\n
    lo: location of legend (0-10, 0:'best')\n
    title: title of figure. (default) Empty string\n
    """    
    fig, ax = plt.subplots(1,1)
    for j in range(len(X_index)):
        if LnLegend==():
            b=DataSet[j].columns[Y_index[j]]
        else:
            b=LnLegend[j]
        ax.scatter(DataSet[j].iloc[:,X_index[j]]
               ,DataSet[j].iloc[:,Y_index[j]]
               ,s=25, marker="o"
               ,label=b)
            
    # Axis label    
    if Xlable == '' or Ylable == '':
        Xlable = DataSet[0].columns[X_index[0]]
        Ylable = DataSet[0].columns[Y_index[0]]
    
    ax.set_xlabel(Xlable, fontsize = 16)
    ax.set_ylabel(Ylable, fontsize = 16)
    
    #plot x- or y-lim
    if Xlim != ():
        ax.set_xlim([Xlim[0],Xlim[1]])
    if Ylim != ():
        ax.set_ylim([Ylim[0],Ylim[1]])    
        
    #ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    
    ax.tick_params(which='both', direction='in')
    ax.legend(loc=lo,frameon=False) 
    
    #plt title
    if title != '' :
        ax.set_title(title, style='italic') 
    
    #ax.axhline(0, color='black', lw=0.5)
    #ax.axvline(0, color='black', lw=0.5)

    del b
    #plt.show()
    return fig,ax


#%%
def MAR (fig,ax,delta_0=200.0,N=6):
    """
    Multiple Andreev reflection (MAR) analysis\n
    Input R or G vs V_sd plot\n
    Search for corresponding SC gap (delta)\n
    \n
    #######################\n
    Input parameters:\n
    fig: existed figure reference\n
    ax : existed axis reference \n
    delta_0=200 #initial delta value, unit: ueV\n
    N=6         #number of voltage line for fitting\n
    \n
    """
    
    plt.subplots_adjust(bottom= 0.25)
   
    #Initial voltage lines based on Detal_0
    ly=[]
    for n in range(1,(N+1)):
        Voltage = delta_0*2/n
        ly.append(ax.axvline(Voltage,lw=1.5, color='r',ls='dashed'))
    
    #Text box to show MAR formula
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    textstr = r'V$_{n}$ = 2$\Delta$/n\n$\Delta$ ='
    delta_unit = r' $\mu$eV'
    Tex=textstr+str(delta_0)+delta_unit
    t1=ax.text(0.1, 0.95, Tex, transform=ax.transAxes, fontsize=12,
             verticalalignment='top',bbox=props,)
    
    #fig = plt.gcf()
    ## slider widget
    axcolor = 'lightgoldenrodyellow'
    axEng = plt.axes([0.22, 0.1, 0.65, 0.04], facecolor=axcolor)
    sEng = Slider(axEng, r'$\Delta$ ($\mu$eV)'
                  ,(delta_0-100),(delta_0+100), valinit=delta_0, valstep=0.5)
    sEng.label.set_fontsize('large')
     
    
    def update(val):
        delta= sEng.val
        Tex=textstr+str(delta)+delta_unit
        t1.set_text(Tex)
        for n in range(1,(N+1)):
            ly[n-1].set_xdata(delta*2/n)
        #fig.canvas.draw_idle()
        fig.canvas.draw()
    
    sEng.on_changed(update)
    #lt.show()
    
#%% Extract Hall density n(Vg) from Rxy(Vg,B) 2D data
def Gate_density(DataSet:pd.DataFrame, 
                 Vg_col:int=1, # Vg column index, unit: V
                 B_col:int=0, # B column index, unit : T
                 Rxy_col:int=2 # Hall resistance column index, unit: Ohm
                 ):
                
    """
    Extract Hall density as function of Vg from Rxy(Vg,B) 2D data\n
    Input parameters:\n
    ======\n
    DataSet : Rxy(Vg,B) pandas data frame\n
    Vg_col  : Vg column index(default:1); unit: V\n
    B_col   : B field column index(default:0) ; unit : T\n
    Rxy_col   : Hall resistance column index (default:2) ; unit : Ohm\n
    Output:\n
    =======\n
    Output : a dataframe with following columns\n
     ['Vg (V)', 'HallCoef (Ohm/T)','n (x10^11 cm^-2)','un_HallCoef','un_n']\n
     un_HallCoef: uncertainty of Hall Coefficient\n
     un_n : uncertainty of 2D Hall density\n
    Example: \n
    =========\n
    df6_Hall=Gate_density(df6,Vg_col=0, B_col =1, Rxy_col=3)
    """
    import scipy.constants as sc
    Cols_name=['Vg (V)', r'Hall_coef ($\Omega$/T)',r'$n_{Hall} (10^{11}$ $cm^{-2}$)','un_HallCoef','un_n']
    df=pd.DataFrame(columns=Cols_name)
    
    # pick up one serious of Vg step from specific B-field value
    Bs = DataSet.iloc[0,B_col]  # extract one exsited B-field value
    
    for i in DataSet[DataSet.iloc[:,B_col]==Bs].iloc[:,Vg_col]: # loop one set of Vg serious
        tmp=DataSet[DataSet.iloc[:,Vg_col]==i]  # Rxy(B) at specific Vg
        f = stats.linregress(tmp.iloc[:,B_col], tmp.iloc[:,Rxy_col])  # Rxy(B) linear fit
        n=1/sc.e/f[0]/1e4/1e11 # n: density (x10^11 cm-2) from Rxy1
        un_n= np.abs(f[4]/sc.e/f[0]/f[0])/1e4/1e11 # f[4] : uncertainty of slope; un_n : uncertainty of density
        d = {Cols_name[0]: i,          # Vg
                Cols_name[1]: [f[0],], #Hall Coefficient
                Cols_name[2]: [n,],    # Hall density
                Cols_name[3]: [f[4],], #uncertainty of Hall Coefficient
                Cols_name[4]: [un_n,]} #uncertainty of Hall density
        if df.empty:
            df = pd.DataFrame(data=d) 
        else:
            df = pd.concat([df,pd.DataFrame(data=d)])
    
    return df




## Example for calculating mobility and its uncersity
## Extract Rxx(B=0) and save in result dataframe: df006_Hall
#tmp=df006[df006.iloc[:,1]==0]
#tmp.reset_index(inplace = True, drop = True)
#df006_Hall['Rxx ($\Omega$)']=tmp['Rxx (k$\Omega$)']*1000
##  Calculate mobility
#g=22 # geometry factor
## mobility
#df006_Hall[r'$\mu_{Hall}$ ($\times$ $10^{3}$ $cm^{2}/Vs$)']=df006_Hall['Hall_coef ($\Omega$/T)']/(df006_Hall['Rxx ($\Omega$)']/g)*1e4/1e3
##un_mu : uncertainty of mobility
#df006_Hall['un_mu'] = df006_Hall['un_HallCoef']/(df006_Hall['Rxx ($\Omega$)']/g)*1e4/1e3


#############################################################################
# WAL fitting function
# Assumption (1) Dirty limit : le << D 
## elastic scattering lenght is smaller than wire width
# Assumption (2) 1D 
##  lm >> D : magnet lenght is larger than wire width
##  Assuming wire width ~ 100 nm, Then B < 66 mT

def WAL_1D_Dirty(B,l_phi,l_so):
    """
    1D WAL model in dirty region\n
    =======\n
    Input:\n
    B: B-field (unit : T)\n
    l_phi : coherence length (unit: nm)\n
    l_so  : spin-orbit scatteringlength (unit: nm)\n
    =======\n
    Output:
        delta_G(B) , unit: (e^2/h) 
    """
    import scipy.constants as sc
    import numpy as np
    ## Wire dimension 
    L = 2e-6 # device length in unit: m
    D = 100e-9 # wire width in unit: m
    return -2/L*(1.5*np.power((1/np.square(l_phi*1e-9) + 4/(3*np.square(l_so*1e-9))+np.square(D*sc.e*B/sc.hbar)/3),-0.5)
                 -0.5*np.power((1/np.square(l_phi*1e-9)+np.square(D*sc.e*B/sc.hbar)/3),-0.5))

###############################################################################
#  Systemetric correction for Data refer to center point
#    SS : symetric part  [Data(+)+Data(-)]/2
#    AS : antisymmetric part  [Data(+)-Data(-)]/2
###############################################################################
def Sym_Correct(
                DataSet:pd.DataFrame,
                Center_RowIdx:int=79,
                X_Col_Idx:int=0,
                Y_Col_Idx:int=1
               ):
    """
    Symetric correction\n
    SS : symetric part [Data(+)+Data(-)]/2 \n
    AS : antisymmetric part [Data(+)-Data(-)]/2  \n
    Assume data points are symetric to cetner point: Y(X+) <-> Y(X-)\n
    If no, please do data interplolation first\n\n
    Input parameters:\n
    =======\n
    DataSet : Pandas data frame\n
    Symetric center point : Row Index\n
    X_Col_Idx : X-array column index(default:0 )\n
    Y_Col_Idx : Y-array column index(default:1 )\n
    Output:\n
    =======\n
     Output dataframe has 3 following columns\n
     [X, Y_SS, Y_AS]\n
     X: X array starts from center point and goes in X+ direction\n
     Y_SS : symetric part \n
     Y_AS : antisymmetric part
    """ 
    Cols_name=[DataSet.columns[X_Col_Idx],
           DataSet.columns[Y_Col_Idx]+'_SS',DataSet.columns[Y_Col_Idx]+'_AS']
    tmp=pd.DataFrame(columns=Cols_name)
    
    Total_stp = min(len(DataSet.iloc[Center_RowIdx:,Y_Col_Idx]),
                    len(DataSet.iloc[:Center_RowIdx,Y_Col_Idx]))
    
    for i in range(Total_stp):
        x= round(DataSet.iloc[Center_RowIdx+i,X_Col_Idx],4)
        y_p= DataSet.iloc[Center_RowIdx+i,Y_Col_Idx]
        y_n= DataSet.iloc[Center_RowIdx-i,Y_Col_Idx]
        y_ss = (y_p+y_n)/2
        y_as = (y_p-y_n)/2
        tmp=tmp.append({Cols_name[0]: x,  # x
                        Cols_name[1]: y_ss, # Symetric part
                        Cols_name[2]: y_as}, ignore_index=True) # Anti-symetric part
    
    return tmp


###########################################################################
#  Calculate 2D density and mobility from Rxy(B) data
################################################################
def Hall(
         DataSet:pd.DataFrame,
         B_Col_Idx:int=0,
         Rxy_Col_Idx:int=1,
         R0:float=1000,  # Rxx(B=0) , unit : Ohm
         g:int=20
        ):
    """
    Calculate 2D Hall density and mobility from Rxy(B) data\n\n
    Input parameter:\n
    ===================\n
    DataSet : Pandas dataframe\n
    B_Col_Idx :   X-array column index in dataframe(default:0 ), unit: T\n
    Rxy_Col_Idx : Y-array column index in dataframe(default:1 ), unit: Ohm\n
    R0 : Rxx(B=0 T), unit: Ohm\n
    g : geometry factor  (l/w)\n
    \n
    Output:\n
    ===================\n
    Result_dict: 
           density: 2D Hall density, unit: x10^11 cm^-2\n
           density_err
           mobility mobility, unit: x10^3 cm^2/Vs\n
           mobility_err
           density_unit
           mobility_unit
           LinregressFit: linregress fitting result
    """ 
    from scipy import stats
    import scipy.constants as sc
    f = stats.linregress(DataSet.iloc[:,B_Col_Idx], #B-array
                         DataSet.iloc[:,Rxy_Col_Idx])  # Rxy(B) array

    Result_dict= {"density": 1/sc.e/f[0]/1e4/1e11, #f[0]: slope; n: density (x10^11 cm-2) from Rxy(B) slope
                  "density_err": np.abs(f[4]/sc.e/f[0]/f[0])/1e4/1e11, # f[4] : uncertainty of slope; un_n : uncertainty of density
                  "mobility": f[0]/(R0/g)*1e4/1e3, #  $\mu$ ($10^{3}$ $cm^{2}$/Vs)'
                  "mobility_err": f[4]/(R0/g)*1e4/1e3, #un_mu : uncertainty of mobility
                  "density_unit": "$10^{11}$ $cm^{-2}$",
                  "mobility_unit": "$10^{3}$ $cm^{2}$/Vs",
                  "LinregressFit": f
                 }

    return Result_dict


class HallAnalysis(object):
    """
    Calculate 2D Hall density and mobility from Rxy(B) data\n\n
    Input parameter:\n
    ===================\n
    B_array :    unit: T\n
    Rxy_array : unit: Ohm\n
    R0 : Rxx(B=0 T), unit: Ohm, default: 500\n
    g : geometry factor  (l/w), default: 1\n
    \n
    Attributes\n
    ===================\n
    fit: linregress result
    
    Method:\n
    ===================\n
    density: Hall density, unit: x10^11 cm^-2\n
    mobility: Hall mobility, unit: x10^3 cm^2/Vs\n
    density_err: uncertainty of density\n
    mobility_err: uncertainty of mobility\n\n
    """ 
    from scipy import stats
    import scipy.constants as sc
    import numpy as np
    def __init__(self, B_array, Rxy_array, R0=500,g=1):
        self.B_array = B_array # B field array
        self.Rxy_array = Rxy_array # Hall resistance array (unit: Ohm)
        self.R0 = R0  # Rxx(B=0) (unit: Ohm)
        self.g = g # Hallbar geometry factor (l/w) 
        self.fit = stats.linregress(B_array, Rxy_array)

    def density(self):
        return 1/sc.e/self.fit.slope/1e4/1e11 #f[0]: slope; n: density (unit: x $10^{11}$ cm-2)

    def density_unit(self):
        return "$10^{11}$ $cm^{-2}$"

    def mobility(self):
        return self.fit.slope/(self.R0/self.g)*1e4/1e3 # $\mu$ (unit: x $10^{3}$ $cm^{2}$/Vs)'

    def mobility_unit(self):
        return "$10^{3}$ $cm^{2}$/Vs"

    def density_err(self):
        return np.abs(self.fit.stderr/sc.e/self.fit.slope/self.fit.slope)/1e4/1e11 # f[4] : uncertainty of slope; un_n : uncertainty of density

    def mobility_err(self):
        return self.fit.stderr/(self.R0/self.g)*1e4/1e3 #un_mu : uncertainty of mobility
        

class SdHAnalysis(object):
    """
    Calculate 2D density from Rxx(B) oscillation\n\n
    Input parameter:\n
    ===================\n
    B0_array :  Rxx(B) dip B-field location (unit: T)\n
    v_array :  assigned filling factor array\n
    \n
    Attributes\n
    ===================\n
    fit: linregress result of filling factor vs 1/B0
    
    Method:\n
    ===================\n
    density: Hall density, unit: x10^11 cm^-2\n
    density_err: uncertainty of density\n
    """ 
    from scipy import stats
    import scipy.constants as sc
    import numpy as np
    def __init__(self, B0_array, v_array):
        self.B0_array = B0_array # Rxx(B) dip B-field location (unit: T)
        self.v_array = v_array # assigned filling factor array 
        self.fit = stats.linregress(np.reciprocal(self.B0_array), self.v_array)

    def density(self):
        return self.fit.slope*sc.e/sc.h/1e4/1e11 # n: density (x10^11 cm-2) from slope

    def density_unit(self):
        return "$10^{11}$ $cm^{-2}$"

    def density_err(self):
        return self.fit.stderr*sc.e/sc.h/1e4/1e11 #   uncertainty of density


###############################################################################
#  1. search the B=0 point from delta_G(B) data
#  2. Pick up even data point from both side of B-field
#  3. Sort data with B-field
###############################################################################
def Center_group(
                 DataSet:pd.DataFrame,
                 X_idx:int = 7, # Column index of B-field
                 Y_idx:int = 6, # Column index of any data
                 B_range:float = 0.13, # desire subgroup B-field range, unit: T
                 B_step:float = 0.002
                ):  # B-resolution in data , unit : T
    """
    1. search the B=0 point (row index) from DataSet \n
    2. Pick up even data point from both side of B-field\n
    3. Sort data with B-field \n
    ===========\n
    DataSet : Pandas data frame\n
    X_idx : Column index of B-field \n
    Y_idx : Column index of desired data \n
    B_range : Desire subgroup B-field range, unit: T
    B_step : B-resolution in data , unit : T\n
    ============\n
    Output : dataframe with following columns \n
     [X-array, Y-array]\n
    """ 
    # search for row index of B=0 T
    center_idx= DataSet.index[round(DataSet.iloc[:,X_idx],3) == 0].tolist()[0]
    nm_stp= int(B_range/B_step) # calculate number of B-field steps
    df = DataSet.iloc[:,[X_idx,Y_idx]][(center_idx-nm_stp):(center_idx+1+nm_stp)]
    df= df.sort_values(by=DataSet.columns[X_idx])
    df= df.reset_index(drop=True)
    return df



###############################################################################
##  Data Cursor class
###############################################################################
class Cursor(object):
    def __init__(self, ax):
        self.ax = ax
        self.lx = ax.axhline(color='k')  # the horiz line
        self.ly = ax.axvline(color='k')  # the vert line

        # text location in axes coords
        self.txt = ax.text(0.7, 0.9, '', transform=ax.transAxes)

    def mouse_move(self, event):
        if not event.inaxes:
            return

        x, y = event.xdata, event.ydata
        # update the line positions
        self.lx.set_ydata(y)
        self.ly.set_xdata(x)

        self.txt.set_text('x=%1.2f, y=%1.2f' % (x, y))
        self.ax.figure.canvas.draw()


class SnaptoCursor(object):
    """
    Like Cursor but the crosshair snaps to the nearest x, y point.
    For simplicity, this assumes that *x* is sorted.
    """

    def __init__(self, ax, x, y):
        self.ax = ax
        self.lx = ax.axhline(color='k')  # the horiz line
        self.ly = ax.axvline(color='k')  # the vert line
        self.x = x
        self.y = y
        # text location in axes coords
        self.txt = ax.text(0.7, 0.9, '', transform=ax.transAxes)

    def mouse_move(self, event):
        if not event.inaxes:
            return

        x, y = event.xdata, event.ydata
        indx = min(np.searchsorted(self.x, x), len(self.x) - 1)
        x = self.x[indx]
        y = self.y[indx]
        # update the line positions
        self.lx.set_ydata(y)
        self.ly.set_xdata(x)

        self.txt.set_text('x=%1.2f, y=%1.2f' % (x, y))
        print('x=%1.2f, y=%1.2f' % (x, y))
        self.ax.figure.canvas.draw()

## Example
#t = np.arange(0.0, 1.0, 0.01)
#s = np.sin(2 * 2 * np.pi * t)
## Cursor
#fig, ax = plt.subplots()
#ax.plot(t, s, 'o')
#cursor = Cursor(ax)
#fig.canvas.mpl_connect('motion_notify_event', cursor.mouse_move)
## SnaptoCursor
#fig, ax = plt.subplots()
#ax.plot(t, s, 'o')
#snap_cursor = SnaptoCursor(ax, t, s)
#fig.canvas.mpl_connect('motion_notify_event', snap_cursor.mouse_move)
#
#plt.show()

#########################################
# find nearest value in the input array #
#########################################
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def gauss(x, H, A, x0, sigma): 
    return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))
