#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 10:33:11 2019

@author: rrabena
"""

import sys
sys.path.append('./ParaView-5.8.0-RC1-MPI-Linux-Python3.7-64bit/lib/python3.7/site-packages')
#from paraview.simple import *
#from math import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from cycler import cycler
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import csv
import os
from astropy.io import ascii
from matplotlib.ticker import AutoMinorLocator
import matplotlib.text as mtext
from durhampy.physics.convert import *
#print(os.getcwd())

#importation des packages et modules durhampy personnalisées crées
import durhampy
from durhampy.testsFunctions.tests import *
from durhampy.postProcessing.plotting import *
from durhampy.postProcessing.calculus import *
from durhampy.data.globals import *
from durhampy.data.units import *
from durhampy.physics.convert import *

#variables globales :

#global raie_dict, H2O_raie_dict, OI_raie_dict, pH2O_exceptions


from matplotlib import pyplot
import numpy as np
import yt

from pygrackle import \
    chemistry_data, \
    setup_fluid_container

from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    sec_per_Myr, \
    cm_per_mpc

"====================PREAMBULE======================"
"""Paramètres Matplotlib"""
mpl.rcdefaults()

Latex=False
font = {'family' : 'serif',
        'size'   : 9,
        'weight' : 'extra bold'}
plt.rc('text', usetex=Latex)
plt.rc('font', **font)
plt.rcParams['text.latex.preamble'] = [
    r'\usepackage[squaren,Gray]{SIunits}',
    r'\usepackage[version=4]{mhchem}']
plt.rcParams['axes.linewidth']=2


plt.rcParams['xtick.major.size'] = 10
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.minor.size'] = 5
plt.rcParams['xtick.minor.width'] = 2

plt.rcParams['ytick.major.size'] = 10
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['ytick.minor.size'] = 5
plt.rcParams['ytick.minor.width'] = 2


"""EXEMPLE D'INPUT"""
#on s'assure d'être dans le bon dossier
racine_dir = "/home/mrabenanahary/MPI-AMRVAC/amrvac/src/grackle/src/python/examples"
os.chdir(racine_dir)

#dossiers racines / modèles cx-ny-bz à ouvrir
input_dir = ['cj-n5-b1']

#liste des vitesses à ouvrir
vit_dir = [[5,10,19]]

#liste des âges à ouvrir
age_dir = [[[50,100,150,200,250,500],[50,100,150,200,250,500],[50,100,150,200,250,500]]]

#type d'ordonnée à afficher. Choix possibles : "Tn","Ti","Te"
choix_ordonnee = "Tn"

#type d'abscisse à afficher. Choix possibles : "tn","ti","z"
choix_abscisse = "ti"




"""1)EXTRACTION DE TEMPERATURE"""
def extraire_temperature(axis,
                         racineDir=racine_dir,
                         direc='cj-n4-b1',
                         vit=10,
                         ag=100,
                         choix_abscisse="ti",
                         choix_ordonnee="Tn",
                         xscale='log',
                         yscale='log',
                         simp_lbl="Y"):

    """

        Extrait 2 listes correspondant aux valeurs contenues dans les 2 colonnes
        d'un fichier-modèle de format nécessairement ASCII, et dont l'une des listes
        est celle du profil de température.

        Cette fonction prend un ensemble de paramètres entrés par l'utilisateur
        pour la guider dans le repérage du modèle de choc à traiter, et
        retourne deux listes:
            - une liste de valeurs d'une variable 'abscisse' associée à choix_abscisse,
            tout en respectant bien l'ordre dans lequel apparaît chaque ligne dans le
            fichier ASCII.
            - une liste de valeurs de températures 'ordonnée" associée à choix_ordonnee,
            et cela tout respectant à coup sûr l'ordre dans lequel apparaît les
            valeurs de la variable 'abscisse', i.e. aussi l'ordre d'apparition
            de ces valeurs dans le fichier ASCII

            (vérification faite de cet ordre le 12/11/2019)

    Utilisation:

        [line,modele] = extraire_temperature(**kwargs):

    Paramètres d'entrée:
    ===================
        Paramètres obligatoires:
        ------------------------

        AUCUN

        Paramètres facultatifs, **kwargs:
        ---------------------------------

        racineDir: chemin racine dans lequel la fonction va chercher
            les modèles. racineDir doit etre un string (indiquant un chemin d'accès)
            et doit être adaptée dans le script main.py pour chaque utilisateur/ordinateur.

        direc: Modèle-dossier  à traiter par la fonction, et nommé obligatoirement sous le
            format de nom cx-ny-bz. direc doit être un string.

        vit: vitesse de choc traitée par la fonction à partir
            d'un dossiers obligatoirement sous le format de nom v"a"-t"b". vit
            doit être un nombre réel

        age_dir: liste des âges de choc traitées par la fonction à partir
            des dossiers obligatoirement sous le format de nom v"a"-t"b". vit_dir
            ne peut prendre qu'une liste de réels à 2 dimensions.

        choix_abscisse: choix de la variable  à représenter en abscisse du
            graphique. Possibilités actuelles : "ti"
            Doit être un string.

        choix_ordonnee: choix de la variable à représenter en ordonnée du
            graphique. Possibilités actuelles : "Tn"
            Doit être un string.

        xscale: type de l'échelle en abscisse ("lin","log",...)

        yscale: type de l'échelle en ordonnée ("lin","log",...)

        simp_lbl: paramètre décidant si les labels affichés sur
        le graphique Matplotlib doivent retenir le préfixe de format cx-ny-bz
        du modèle (d'une façon à ce qu'il ne puisse pas y avoir mélange de labels
        entre les différentes courbes), ou si l'on peut occulter ce préfixe
        pour davantage de visibilité (permets aussi de s'assurer que le
        label entré pour chaque courbe correspond bien au modèle traité à chaque
        étape de chaque boucle)

    Paramètres de sortie:
    ===================

        line: variable contenant le retour de la fonction axis.plot(axis,**kwargs) de
            Matplotlib via une définition obligatoire de axis via les lignes
                plt.figure()
                ax = fig.add_subplot(1,1,1)
            Cette variable est nécessaire pour établir correctement
            les labels de la légende dans la fonction plot_temperature(**kwargs)
            de manière hiérarchisé (rangé dans la légende par modèle cx-ny-bz)


        modele: retourne le label du modèle du modèle de choc traité par cette
            fonction extraire_temperature(axis,**kwargs). Comme line, ce paramètre
            est utile pour la hiérarchisation correcte des labels des modèles
            de chocs dans plot_temperature(**kwargs).


    e.g. d'utilisation :
    ===================

            - plot_temperature(input_dir=['cj-n6-b1'],vit_dir = [[5,10,19]],age_dir = [[[50,100,150,200,250,500],[50,100,150,200,250,500],[50,100,150,200,250,500]]])
            - plot_temperature(input_dir=['cj-n6-b1'],vit_dir = [[5]],age_dir = [[[50,100,150,200,250,500]]])
            - plot_temperature(input_dir=['cj-n4-b1'],vit_dir = [[5]],age_dir = [[[50,100,200,500]]])
    """

    modele = 'v' + str(vit) + '-t' + str(ag)
    modele_chemin = racineDir + '/' + direc + \
    '/' + modele
    #on travaille dans le dossier c"x"-n"y"-b"z"/v"a"-t"b"
    """ 12/11/19 : j'ai checké, et les indices de colonne
    dans mhd_phys.out pour extraire:
        - la température Tn est i=14, resp. i=15,16 pour Ti,Te
    """
    os.chdir(modele_chemin)
    input_file=ascii.read('mhd_phys.out')
    if(choix_ordonnee =="Tn"):
        temperature = [row[14] for row in  input_file]
    if(choix_abscisse =="ti"):
        x_abscisse = [row[2] for row in  input_file]


    #ax.set_xlim(7,1000)
    [line] = axis.plot(np.array(x_abscisse),np.array(temperature),'-',label=modele)
    axis.set_xscale(xscale)
    axis.set_yscale(yscale)

    if(simp_lbl != "Y"):
        modele = direc +'-'+modele
        print(">Le modèle actuellement traité par extraire_temperature() est :\n",modele)

    return [line,modele]


"""2)PLOT DE TEMPERATURE"""
def plot_temperature(racine_dir=racine_dir, # racine dans lequel la fonction va chercher les modèles
                     input_dir = ['cj-n4-b1'],            # liste des modèles-dossiers cx-ny-bz à traiter
                     vit_dir = [[10]],                    # liste des vitesses de choc traitées dans v"a"-t"b"
                     age_dir = [[[50,100,200,500]]],      # liste des âges de choc traitées dans v"a"-t"b"
                     choix_abscisse="ti",                 # choix du type de variable-abscisse traitée (temps, distance, etc...)
                     choix_ordonnee="Tn",                 # choix du type de variable-ordonnée température traitée
                     xScale='log',                        # choix du type d'échelle en abscisse
                     yScale='log',                        # choix du type d'échelle en ordonnée
                     xmin=10,                             # échelle min en x
                     xmax=1500,                           # échelle max en x
                     do_you_save = "Y",                   # si on souhaite sauver le plot dans un fichier
                     save_format="pdf",                   # format du fichier output
                     ncolonnes=1,                         # nombre de colonnes de la légende du graphique
                     simplified_labels="Y"):
    """

        Affiche (et éventuellement sauvegarde) un graphique y=f(x)
        de chaque modèle entré en input.

        Cette fonction prend un ensemble de paramètres entrés par l'utilisateur
        pour la guider dans le repérage des modèles de chocs à traiter, et
        affiche (voire sauvegarde dans un fichier) la figure obtenue
        sous Matplotlib.

    Utilisation:

        plot_temperature(**kwargs):

    Paramètres d'entrée:
    ===================
        Paramètres obligatoires:
        ------------------------

        AUCUN

        Paramètres facultatifs, **kwargs:
        ---------------------------------

        racine_dir: chemin racine dans lequel la fonction va chercher
            les modèles. racine_dir doit etre un string (indiquant un chemin d'accès)
            et doit être adaptée dans le script main.py pour chaque utilisateur/ordinateur.

        input_dir: Liste des modèles-dossiers à traiter par la fonction, et nommés
            chacun obligatoirement sous le format de nom cx-ny-bz.
            input_dir ne peut prendre qu'une liste de string à 1 dimension.

        vit_dir: liste des vitesses de choc traitées par la fonction à partir
            des dossiers obligatoirement sous le format de nom v"a"-t"b". vit_dir
            ne peut prendre qu'une liste de réels à 2 dimensions.

        age_dir: liste des âges de choc traitées par la fonction à partir
            des dossiers obligatoirement sous le format de nom v"a"-t"b". vit_dir
            ne peut prendre qu'une liste de réels à 3 dimensions.

        choix_abscisse: choix de la variable  à représenter en abscisse du
            graphique. Possibilités actuelles : "ti"
            Doit être un string.

        choix_ordonnee: choix de la variable à représenter en ordonnée du
            graphique. Possibilités actuelles : "Tn"
            Doit être un string.

        xScale: type de l'échelle en abscisse ("lin","log",...)

        yScale: type de l'échelle en ordonnée ("lin","log",...)

        xmin: échellle minimale en abscisse. Prend une valeur de nombre réel.

        xmax: échellle maximale en abscisse. Prend une valeur de nombre réel.

        do_you_save: bascule pour sauvegarder ou non sous forme de fichier le
            graphique. Valeurs possibles : "Y" ou "N" (ou autre)

        save_format: nom du format de fichier sous lequel sauvegarder le graphique.
            Par défaut, on utilise le format pdf (valeur "pdf" entrée), mais on
            pourrait entrer par e.g. "png"

        ncolonnes: nombre de colonnes à afficher dans la légende du graphique.
            Prend un nombre entier.

        simplified_labels: paramètre décidant si les labels affichés sur
        le graphique Matplotlib doivent retenir le préfixe de format cx-ny-bz
        du modèle (d'une façon à ce qu'il ne puisse pas y avoir mélange de labels
        entre les différentes courbes), ou si l'on peut occulter ce préfixe
        pour davantage de visibilité (permets aussi de s'assurer que le
        label entré pour chaque courbe correspond bien au modèle traité à chaque
        étape de chaque boucle)

    Paramètres de sortie:
    ===================

    AUCUN

    e.g. d'utilisation :
    ===================

            - plot_temperature(input_dir=['cj-n6-b1'],vit_dir = [[5,10,19]],age_dir = [[[50,100,150,200,250,500],[50,100,150,200,250,500],[50,100,150,200,250,500]]])
            - plot_temperature(input_dir=['cj-n6-b1'],vit_dir = [[5]],age_dir = [[[50,100,150,200,250,500]]])
            - plot_temperature(input_dir=['cj-n4-b1'],vit_dir = [[5]],age_dir = [[[50,100,200,500]]])
    """

    if(do_you_save == "Y") : fileTitle=""
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.grid(True,axis='both',alpha=0.75,zorder=0)
    if choix_abscisse == "ti" :ax.set_xlabel(r'$t_i$ (ans)')
    if choix_ordonnee == "Tn" : ax.set_ylabel(r'$T_n$ (K) ')

    handl=[]
    labl=[]

    for i,directory in enumerate(input_dir):

        labl.append(r"\underline{\textbf{%s :}}" % directory)
        handl.append(plt.plot([],marker="", ls="")[0])
        if(do_you_save == "Y") :
            if(i==0):
                fileTitle += directory
                fileTitle += "-v"+str(np.min(vit_dir[i]))+"TO"+str(np.max(vit_dir[i]))
            else:
                fileTitle += "_"+directory
                fileTitle += "_v"+str(np.min(vit_dir[i]))+"TO"+str(np.max(vit_dir[i]))

        age_min,age_max=9.99e+100,0

        for j,vitesse in enumerate(vit_dir[i]):

            age_min_temp = np.min(age_dir[i][j])
            age_max_temp = np.max(age_dir[i][j])

            for age in age_dir[i][j]:

                age_min = np.min([age_min,age_min_temp])
                age_max = np.max([age_max,age_max_temp])

                [line,lbl] = extraire_temperature(axis = ax,
                                     racineDir = racine_dir,
                                     direc = directory,
                                     vit = vitesse,
                                     ag = age,
                                     choix_abscisse = "ti",
                                     choix_ordonnee = "Tn",
                                     xscale=xScale,
                                     yscale=yScale,
                                     simp_lbl=simplified_labels)
                if(simplified_labels := "Y"):
                    try:
                        label_plot_temperature = (directory+'-v' + str(vitesse) + '-t' + str(age))
                        labels_identiq = (label_plot_temperature ==lbl)
                        if(labels_identiq==False):
                            #print((directory+'-v' + str(vitesse) + '-t' + str(age)), lbl)
                            raise ValueError("Souci de labelisation du modèle 0:-v1:-t2:".format(directory, str(vitesse),str(age)))
                    except ValueError as errormsg:
                        print(">>>ERREUR : Le label actuel retourné par extraire_temperature() est :",lbl)
                        print("...tandis que le label actuel traité par plot_temperature() est :",label_plot_temperature)
                        print(errormsg)
                        print("\n")
                    else:
                        print("-->Aucun problème !! Le label {0:} de cette courbe affichée par extraire_temperature()\n correspond bien à celui {1:} du modèle de choc \n traité par plot_temperature()".format(lbl,(directory+'-v' + str(vitesse) + '-t' + str(age))))
                        print("\n")

                handl.append(line)
                labl.append(lbl)

        if(do_you_save == "Y") :
            if(i==0):
                fileTitle += "-t"+str(int(age_min))+"TO"+str(int(age_max))
            else:
                fileTitle += "_t"+str(int(age_min))+"TO"+str(int(age_max))

    plt.legend(handl,labl,fontsize=13,ncol=ncolonnes, bbox_to_anchor=(1, 0.5),loc="lower left")
    ax.set_xlim(xmin,xmax)
    if(do_you_save == "Y") : print("Saving "+fileTitle+"...")
    os.chdir(racine_dir)
    if(do_you_save == "Y") : plt.savefig("profils_temperatures/"+fileTitle+"."+save_format, bbox_inches='tight')
    plt.show()
    plt.close()
       #     extraire_temperature(racine_dir,directory,vitesse,age)






def extraire_raie(racine_dir=racine_dir,
                  molecule_choisie="CO",            #choix de la molécule dont est extraite le TDV de la raie
                  Jupper=16,                        ##si molecule_choisie=="CO", J_up de la transition
                  input_dir='cj-n4-b1',             #choix du modèle cx-ny-bz traité
                  vitesse=10,                       #choix de la vitesse de choc traité dans la nomenclature v"a"-t"b"
                  age=100,                          #choix de l'âge du choc traité dans la nomenclature v"a"-t"b"
                  CO_colonne_raie = "colonne_TDV",  #si molecule_choisie=="CO", alors ici est le header de la colonne de raie_dict d'où on on extrait la TDV
                  nuH2OGHz = 2547,                  #si molecule_choisie=="o/pH2O", fréquence en GHz de la transition
                  H2O_colonne_raie = "colonne_TDV", #si molecule_choisie=="o/pH2O", choix de la colonne a saisir dans o/pH2O_int_int.out
                  lambda_microns="63m",             #si molecule_choisie=="OI", choix de la colonne a saisir dans intensity.out correspondant à l'atome d'oxygène atomique
                  vup = 0,
                  tester_verifie_modele=True,
                  print_tests_results = False,
                  type_of_output = "default",
                  arguments_to_pass = {}):
    """1) Cette fonction a été testée et validée par Mialy comme marchant sans soucis \n (au moins pour l'extraction de la raie de la molécule CO) \n le 13/11/2019 à 11h35
       2) Cette fonction a été testée et validée par Mialy comme marchant sans soucis \n (au moins pour l'extraction de la raie des molécules CO et H2O) \n le 13/11/2019 à 16h49
       3) Cette fonction a été testée et validée par Mialy comme marchant sans soucis \n (au moins pour l'extraction de la raie des molécules CO, H2O et OI) \n le 01/12/2019 à 22h48
       4) Cette fonction a été testée et validée par Mialy comme marchant sans soucis \n (au moins pour l'extraction de la raie des molécules CO, H2O, OI et H2) \n le 02/12/2019 à 21h28
    """
    os.chdir(racine_dir)
    os.chdir(input_dir)
    modele = 'v' + str(vitesse) + '-t' + str(age)
    modele_complet = input_dir+'-'+modele
    os.chdir(modele)
    nomDeFonction = 'extraire_raie'

    if(molecule_choisie == "H2"): #0) partie traitant la molécule H2
        info_mhd_H2 = verifieModele(actual_dir=modele,inputDir=input_dir,vit=vitesse,ag=age, nom_de_fonction=nomDeFonction, return_info_mhd=True, print_tests=print_tests_results
               , valeurs_param_tester = {'input files':{},
                                          'shock parameters': {},
                                          'environmental parameters': {'Zeta(s-1)':5.00E-17},
                                          'grains parameters': {},
                                          'excitation & cooling': {},
                                          'integration parameters': {},
                                          'DVODE parameters': {},
                                          'outputs parameters': {},
                                          'H2 molecule': {},
                                          'elemental abundances (gas + mantles + PAH)': {}
                                          }          )
        #print(info_mhd_H2)
        H2_data_type = info_mhd_H2['outputs parameters']['H2_out']#type of data in H2_line.out -  'local' (erg/s/cm3) or 'integrated' (erg/s/cm2/sr)
        #print(H2_data_type)
        input_file = ascii.read('excit.out',header_start=0,data_start=1,fast_reader=False)
        #print(input_file)
        sortir = False
        for row in input_file:
            #print(row)
            if(row[0]==vup):
                if(row[1]==Jupper):
                    rTDV = row[-1]   #e.g. ln(N/g)
                    rXraie = row[-2] #Eup en K
                    sortir = True
                    if(sortir):
                        break
            #print(row)
    if(molecule_choisie == "CO"): #1) partie traitant la molécule CO
        if(tester_verifie_modele==True) : verifieModele(actual_dir=modele,inputDir=input_dir,vit=vitesse,ag=age, nom_de_fonction=nomDeFonction, print_tests=print_tests_results
               , valeurs_param_tester = {'input files':{},
                                          'shock parameters': {},
                                          'environmental parameters': {'Zeta(s-1)':5.00E-17},
                                          'grains parameters': {},
                                          'excitation & cooling': {},
                                          'integration parameters': {},
                                          'DVODE parameters': {},
                                          'outputs parameters': {},
                                          'H2 molecule': {},
                                          'elemental abundances (gas + mantles + PAH)': {}
                                          }          )
        input_file=ascii.read('CO_int_int.out',header_start=None,data_start=2)
        try:
            assert CO_colonne_raie in ["colonne_TDV","colonne_TDV2"]
        except AssertionError:
            print(">>Attention, dans le modèle {0:} dans la fonction extraire_raie(), \n CO_colonne_raie (la variable choisissant la colonne TDV dans \n CO_int_int.out) n'est pas dans [\"colonne_TDV\",\"colonne_TDV2\"]".format(modele_complet))
        i_Jupper = Jupper - 1 #pour saisir dans input_file la ligne correspondant à Jupper
        j_Jupper = raie_dict[molecule_choisie]["colonne_Jupper"]                    # pour saisir à la ligne i_Jupper la colonne correspondant à Jupper
        j_TDV = raie_dict[molecule_choisie][CO_colonne_raie]                        # pour saisir à la ligne i_Jupper la colonne correspondant au TDV
        rJupper = int(input_file[i_Jupper][j_Jupper])
        try:
            assert rJupper==Jupper
        except AssertionError:
            print(">>Attention, dans le modèle {0:} dans la fonction extraire_raie(), Jupper entré par l'utilisateur vaut {1:} et est != rJupper retrouvé par lecture du fichier CO_int_int.out, rJupper valant ici {2:} (ligne du fichier:\n {3:}\n)\n ".format(modele_complet,Jupper,rJupper,input_file[i_Jupper]))
        else:
            if(print_tests_results) : print(">Modèle {0:}\n dans extraire_raie():\n Aucun problème ! de correspondance de l'indice de ligne parcouru dans le fichier CO_int_int.out : \n Jupper==rJupper".format(modele_complet))
        print("\n")
        rTDV = input_file[i_Jupper][j_TDV] #TdV en K.km/s
        rXraie = rJupper                   #Jup


    elif ((molecule_choisie == "oH2O") or (molecule_choisie == "pH2O")): #2) partie traitant la molécule H2O
        if(tester_verifie_modele==True) : verifieModele(actual_dir=modele,inputDir=input_dir,vit=vitesse,ag=age, nom_de_fonction=nomDeFonction, print_tests=print_tests_results
                  , valeurs_param_tester = {'input files':{},
                                          'shock parameters': {},
                                          'environmental parameters': {'Zeta(s-1)':5.00E-17},
                                          'grains parameters': {},
                                          'excitation & cooling': {},
                                          'integration parameters': {},
                                          'DVODE parameters': {},
                                          'outputs parameters': {},
                                          'H2 molecule': {},
                                          'elemental abundances (gas + mantles + PAH)': {}
                                          })
        H2O_file = molecule_choisie + '_int_int.out'                                #nom du fichier de format 'o/pH2O_int_int.out' passe-partout corrigé des caractères non-ascii
        H2O_file_b = molecule_choisie + '_int_int_b.out'                            #nom du fichier de format 'o/pH2O_int_int_b.out' passe-partout corrigé des caractères non-ascii
        input_file=ascii.read(H2O_file_b)

        """Extraction de la fréquence et du TDV associé"""
        #on détermine dans o/pH2O_int_int_b.out la ligne correspondant à la fréquence entrée
        #<=> fréquence la plus proche de celle entrée par l'utilisateur"""
        i_freq = H2O_raie_dict[molecule_choisie][nuH2OGHz]                          #indice désignant la ligne (sans le Header) dans o/pH2O_int_int_b.out

        #on commence par déterminer dans 'o/pH2O_int_int_b.out' les indices de colonnes associées à la fréquence et au TDV
        #print(raie_dict)
        j_freq = raie_dict[molecule_choisie]["colonne_nuGHz"]                       # pour saisir à la ligne i_freq la colonne correspondant à la fréquence
        j_TDV = raie_dict[molecule_choisie][H2O_colonne_raie]

        #On extrait la fréquence et le TDV contenue dans le fichier 'o/pH2O_int_int_b.out'
        rTDV = input_file[i_freq][j_TDV]  #TdV en K.km/s                                           #TDV retournée par la fonction, si tout a bien été fait
        rXraie = input_file[i_freq][j_freq] #nu en GHz                                         #fréquence retournée par la fonction, si tout a bien été fait


        """==============================DEBUT (FACULTATIF)=============================="""
        """Test de bonne conversion des fichiers o/pH2O_int_int.out ---> o/pH2O_int_int_b.out"""
        with open(H2O_file,"rb") as input_file_2: #obligé de mettre le mode "rb" depuis le passage en Python3+ (bug Python)
            input_file_2.readline()                                                 #on zappe le commentaire du fichier o/pH2O_int_int.out introduit par la modification de B.G.
            input_file_2.readline()                                                 #on zappe le header du fichier o/pH2O_int_int.out

            #on zappe les lignes précédant i_freq premières
            for i_ligne in range(0,i_freq):
                input_file_2.readline()
            if(i_ligne>=11):input_file_2.readline()
            #on sauvegarde le contenu de la ligne i_freq
            ligne = input_file_2.readline()
            ligne_splitted = ligne.split()

            """tant que le bug lié à H2O (sorti par le script du code LVG d'Antoine)
            n'aura pas corrigé : retire ici la colonne 0 (celle qui présente des anomalies
            de caractères illisibles) :"""

            ligne_splitted_corriged =  np.double(ligne_splitted[1:])

            """On vérifie dans la suite que la ligne i_freq dans o/pH2O_int_int_b.out
            correspond est identique en fréquence et TdV à la ligne i_freq de o/pH2O_int_int_b.out :"""

            #1)Tests sur la fréquence
            #on teste la correspondance de la fréquence entrée par l'utilisateur et celle trouvée à l'issue de la lecture de 'o/pH2O_int_int_b.ou'
            test_freq_1 = True #ceci est-vrai si le test suivant est vérifié
            try:
                if (nuH2OGHz in pH2O_exceptions):
                    test_droite = int(nuH2OGHz[:-1])#nuH2OGHz
                else:
                    test_droite = nuH2OGHz
                assert test_droite==ligne_splitted_corriged[0]

            except AssertionError:
                print("a) Attention, dans le modèle {0:} dans la fonction extraire_raie(), la fréquence entrée par l'utilisateur vaut {6:} GHz et est != la fréquence retrouvée par lecture du fichier '{7:}', qui, elle, vaut ici {2:} GHz \n (ligne du fichier '{7:}' :\n {3:} ou\n {4:}\n ou {5:} \n)".format(modele_complet,rXraie,ligne_splitted_corriged[0],ligne,ligne_splitted,ligne_splitted_corriged,nuH2OGHz,H2O_file))
                test_freq_1 = False
            else:
                if(print_tests_results) : print(">Modèle {0:}\n dans extraire_raie():\n Aucun problème ! de correspondance entre les fréquences entrée par l'utilisateur et celle trouvée par la fonction dans '{1:}'".format(modele_complet,H2O_file))

            #on teste la correspondance des lignes entre les fichiers 'o/pH2O_int_int_b.out ' et 'o/pH2O_int_int_b.ou'
            test_freq_2 = True                                                      #ceci est-vrai si le test suivant est vérifié
            try:
                assert ligne_splitted_corriged[0]==rXraie
            except AssertionError:
                print("b) Attention, dans le modèle {0:} dans la fonction extraire_raie(), la fréquence trouvée dans '{8:}' vaut {1:} GHz et est != la fréquence retrouvée par lecture du fichier '{7:}', qui, elle, vaut ici {2:} GHz \n (ligne du fichier '{7:}' :\n {3:} ou\n {4:}\n ou {5:} \n ligne du fichier '{8:}' :\n {6:}\n )\n) ".format(modele_complet,rXraie,ligne_splitted_corriged[0],ligne,ligne_splitted,ligne_splitted_corriged,input_file[i_freq],H2O_file,H2O_file_b))
                test_freq_2 = False
            else:
                if(print_tests_results) : print(">>Modèle {0:}\n dans extraire_raie():\n Aucun problème ! de correspondance entre les fréquences trouvées dans  '{2:}' et '{1:}'".format(modele_complet,H2O_file,H2O_file_b))
            test3 = test_freq_1 and test_freq_2
            # si au moins l'un des 2 tests précédents est faux, on signale à l'utilisateur
            # de faire attention à l'erreur durant ce modèle
            if (test3 is False) : print("/!\|^|----ATTENTION, ERREUR CI-DESSUS!!!!")
            print("\n")

            #2)Tests sur la TDV

            #on teste la correspondance des lignes entre les fichiers 'o/pH2O_int_int_b.out ' et 'o/pH2O_int_int_b.ou'
            test_freq_3 = True #ceci est-vrai si le test suivant est vérifié
            try:
                assert ligne_splitted_corriged[1]==rTDV
            except AssertionError:
                print("b) Attention, dans le modèle {0:} dans la fonction extraire_raie(), la TDV trouvée dans '{8:}' vaut {1:} et est != la TDV retrouvée par lecture du fichier '{7:}', qui, elle, vaut ici {2:} \n (ligne du fichier '{7:}' :\n {3:} ou\n {4:}\n ou {5:} \n ligne du fichier '{8:}' :\n {6:}\n )\n ".format(modele_complet,rTDV,ligne_splitted_corriged[1],ligne,ligne_splitted,ligne_splitted_corriged,input_file[i_freq],H2O_file,H2O_file_b))
                test_freq_3 = False
            else:
                if(print_tests_results) : print(">>Modèle {0:}\n dans extraire_raie():\n Aucun problème ! de correspondance entre les TDV trouvées dans  '{2:}' et '{1:}'".format(modele_complet,H2O_file,H2O_file_b))

            """ -13/11/2019 à 15h57, les 3 tests précédents réalisés par Mialy ont montré que les strings affichés dans le try...except...else fonctionnent bien
            """

            # si au moins l'un des 2 tests précédents est faux, on signale à l'utilisateur
            # de faire attention à l'erreur durant ce modèle
            if (test_freq_3 is False) : print("/!\|^|----ATTENTION, ERREUR CI-DESSUS!!!!")
            print("\n")

            #print(rXraie,rTDV)
    elif (molecule_choisie == "OI"): #3) partie traitant l'atome OI
        if(tester_verifie_modele==True) : verifieModele(actual_dir=modele,inputDir=input_dir,vit=vitesse,ag=age, nom_de_fonction=nomDeFonction, print_tests=print_tests_results
                  , valeurs_param_tester = {'input files':{},
                                          'shock parameters': {},
                                          'environmental parameters': {'Zeta(s-1)':5.00E-17},
                                          'grains parameters': {},
                                          'excitation & cooling': {},
                                          'integration parameters': {},
                                          'DVODE parameters': {},
                                          'outputs parameters': {},
                                          'H2 molecule': {},
                                          'elemental abundances (gas + mantles + PAH)': {}
                                          })
        input_file=ascii.read('intensity.out',header_start=0,data_start=1)
        """Extraction de la longueur d'onde """
        #on retourne dans OI_raie_dict la longueur d'onde d'après Yang et al. 2017 de la raie de OI entrée
        i_last = -1                                                                 #indice désignant la ligne (sans le Header) dans intensity.out correspondant au dernier temps t_i du profil enregistré par le code de choc

        #on commence par déterminer dans 'intensity.out' les indices de colonnes associées au TDV de OI
        j_last = i_freq = OI_raie_dict[molecule_choisie][lambda_microns][0]         # pour saisir à la ligne i_last la colonne correspondant au TDV

        #On extrait la longueur d'onde en microns et le TDV contenu dans le fichier 'intensity.out'
        rTDV = input_file[i_last][j_last]                           #brillance de surface (en erg.cm-2.s-1.sr-1) retournée par la fonction, si tout a bien été fait
        rXraie =  OI_raie_dict[molecule_choisie][lambda_microns][1] #longueur d'onde retournée par la fonction, si tout a bien été fait
        """Attention, la quantité écrite dans intensity.out aux colonnes 'O(63m)' et 'O(145m)'
        sont en fait des brillances de surface en erg.cm-2.s-1.sr-1 (ou integrated intensity I_\nu \Delta \nu en anglais dans Flower 2015)
        et non pas des intensités intégrées TDV, contrairement aux colonnes du fichier CO_int_int.out ou o/pH2O_int_int.out
        qui sont bien des intensités intégrées en K.km.s-1 (ou integrated line temperature en anglais dans Flower 2015)"""
    if(molecule_choisie=='H2'):
        return rXraie,rTDV,H2_data_type
    else:
        return rXraie,rTDV



    """==============================FIN (FACULTATIF)=============================="""

def extraire_rapport(racineDir=racine_dir,
                  molecule_choisie_num="CO",            #choix de la molécule dont est extraite le TDV de la raie (numérateur)
                  molecule_choisie_den="CO",            #choix de la molécule dont est extraite le TDV de la raie (dénominateur)
                  Jupper_num=16,                        ##si molecule_choisie=="CO", J_up de la transition (numérateur)
                  Jupper_den=16,                        ##si molecule_choisie=="CO", J_up de la transition (dénominateur)
                  input_dir_num='cj-n4-b1',             #choix du modèle cx-ny-bz traité (numérateur)
                  input_dir_den='cj-n4-b1',             #choix du modèle cx-ny-bz traité (dénominateur)
                  vitesse_num=10,                       #choix de la vitesse de choc traité dans la nomenclature v"a"-t"b" (numérateur)
                  vitesse_den=10,                       #choix de la vitesse de choc traité dans la nomenclature v"a"-t"b" (dénominateur)
                  age_num=100,                          #choix de l'âge du choc traité dans la nomenclature v"a"-t"b" (numérateur)
                  age_den=100,                          #choix de l'âge du choc traité dans la nomenclature v"a"-t"b" (dénominateur)
                  CO_colonne_raie_num = "colonne_TDV",  #si molecule_choisie=="CO", alors ici est le header de la colonne de raie_dict d'où on on extrait la TDV (numérateur)
                  CO_colonne_raie_den = "colonne_TDV",  #si molecule_choisie=="CO", alors ici est le header de la colonne de raie_dict d'où on on extrait la TDV (dénominateur)
                  nu_num = 2547,                  #si molecule_choisie=="o/pH2O", fréquence en GHz de la transition (numérateur)
                  nu_den = 2547,                  #si molecule_choisie=="o/pH2O", fréquence en GHz de la transition (dénominateur)
                  H2O_colonne_raie_num = "colonne_TDV", #si molecule_choisie=="o/pH2O", choix de la colonne a saisir dans o/pH2O_int_int.out (numérateur)
                  H2O_colonne_raie_den = "colonne_TDV", #si molecule_choisie=="o/pH2O", choix de la colonne a saisir dans o/pH2O_int_int.out (dénominateur)
                  lambda_microns_num="63m",             #si molecule_choisie=="OI", choix de la colonne a saisir dans intensity.out correspondant à l'atome d'oxygène atomique (numérateur)
                  lambda_microns_den="63m",             #si molecule_choisie=="OI", choix de la colonne a saisir dans intensity.out correspondant à l'atome d'oxygène atomique (dénominateur)
                  vup_num = 0,
                  vup_den = 0,
                  testerVerifieModele=True,
                  printTestsResults = False,
                  typeOfOutput = "default",
                  argumentsToPass = {}):

    function_name = "extraire_rapport"

    modele_num = 'v' + str(vitesse_num) + '-t' + str(age_num)
    modele_num_complet = input_dir_num+'-'+modele_num

    raie_num = extraire_raie(racine_dir=racineDir,
                  molecule_choisie=molecule_choisie_num,            #choix de la molécule dont est extraite le TDV de la raie
                  Jupper=Jupper_num,                        ##si molecule_choisie=="CO", J_up de la transition
                  input_dir=input_dir_num,             #choix du modèle cx-ny-bz traité
                  vitesse=vitesse_num,                       #choix de la vitesse de choc traité dans la nomenclature v"a"-t"b"
                  age=age_num,                          #choix de l'âge du choc traité dans la nomenclature v"a"-t"b"
                  CO_colonne_raie = CO_colonne_raie_num,  #si molecule_choisie=="CO", alors ici est le header de la colonne de raie_dict d'où on on extrait la TDV
                  nuH2OGHz = nu_num,                  #si molecule_choisie=="o/pH2O", fréquence en GHz de la transition
                  H2O_colonne_raie = H2O_colonne_raie_num, #si molecule_choisie=="o/pH2O", choix de la colonne a saisir dans o/pH2O_int_int.out
                  lambda_microns=lambda_microns_num,             #si molecule_choisie=="OI", choix de la colonne a saisir dans intensity.out correspondant à l'atome d'oxygène atomique
                  vup = vup_num,
                  tester_verifie_modele=testerVerifieModele,
                  print_tests_results = printTestsResults,
                  type_of_output = typeOfOutput,
                  arguments_to_pass = argumentsToPass)

    modele_den = 'v' + str(vitesse_den) + '-t' + str(age_den)
    modele_den_complet = input_dir_den+'-'+modele_den

    raie_den = extraire_raie(racine_dir=racineDir,
                  molecule_choisie=molecule_choisie_den,            #choix de la molécule dont est extraite le TDV de la raie
                  Jupper=Jupper_den,                        ##si molecule_choisie=="CO", J_up de la transition
                  input_dir=input_dir_den,             #choix du modèle cx-ny-bz traité
                  vitesse=vitesse_den,                       #choix de la vitesse de choc traité dans la nomenclature v"a"-t"b"
                  age=age_den,                          #choix de l'âge du choc traité dans la nomenclature v"a"-t"b"
                  CO_colonne_raie = CO_colonne_raie_den,  #si molecule_choisie=="CO", alors ici est le header de la colonne de raie_dict d'où on on extrait la TDV
                  nuH2OGHz = nu_den,                  #si molecule_choisie=="o/pH2O", fréquence en GHz de la transition
                  H2O_colonne_raie = H2O_colonne_raie_den, #si molecule_choisie=="o/pH2O", choix de la colonne a saisir dans o/pH2O_int_int.out
                  lambda_microns=lambda_microns_den,             #si molecule_choisie=="OI", choix de la colonne a saisir dans intensity.out correspondant à l'atome d'oxygène atomique
                  vup = vup_den,
                  tester_verifie_modele=testerVerifieModele,
                  print_tests_results = printTestsResults,
                  type_of_output = typeOfOutput,
                  arguments_to_pass = argumentsToPass)

    #print(raie_num)
    #print(raie_den)

    rapport_tab = {}

    rapport_tab["raie numerateur"] = {}
    rapport_tab["raie numerateur"]["X"] = raie_num[0]
    rapport_tab["raie numerateur"]["Y"] = raie_num[1]
    rapport_tab["raie denominateur"] = {}
    rapport_tab["raie denominateur"]["X"] = raie_den[0]
    rapport_tab["raie denominateur"]["Y"] = raie_den[1]
    rapport_tab["ratio num/den"] = raie_num[1]/raie_den[1]
    rapport_tab["ratio den/num"] = raie_den[1]/raie_num[1]

    try:
        assert (np.abs(rapport_tab["ratio den/num"]-1./rapport_tab["ratio num/den"])<1e-6)
    except AssertionError:
        raise AssertionError(" Attention, dans la fonction {0:}(), le ratio num/den = {1:}({2:})/{3:}({4:}) vaut {5:}, le ratio den/num = {3:}({4:})/{1:}({2:}) vaut {6:} alors qu\' on  a 1/(num/den) valant {7:}".format(function_name,modele_num_complet,molecule_choisie_num,modele_den_complet,molecule_choisie_den,rapport_tab["ratio num/den"],rapport_tab["ratio den/num"],1./rapport_tab["ratio num/den"]))

    return rapport_tab

def plotRapport(figure, axes,
                  racineDir=racine_dir,
                  molecule_choisie_num="CO",            #choix de la molécule dont est extraite le TDV de la raie (numérateur)
                  molecule_choisie_den="CO",            #choix de la molécule dont est extraite le TDV de la raie (dénominateur)
                  Jupper_num=16,                        ##si molecule_choisie=="CO", J_up de la transition (numérateur)
                  Jupper_den=16,                        ##si molecule_choisie=="CO", J_up de la transition (dénominateur)
                  input_dir='cj-n4-b1',             #choix du modèle cx-ny-bz traité
                  vitesse=10,                       #choix de la vitesse de choc traité dans la nomenclature v"a"-t"b"
                  age=100,                          #choix de l'âge du choc traité dans la nomenclature v"a"-t"b"
                  CO_colonne_raie_num = "colonne_TDV",  #si molecule_choisie=="CO", alors ici est le header de la colonne de raie_dict d'où on on extrait la TDV (numérateur)
                  CO_colonne_raie_den = "colonne_TDV",  #si molecule_choisie=="CO", alors ici est le header de la colonne de raie_dict d'où on on extrait la TDV (dénominateur)
                  nu_num = 2547,                  #si molecule_choisie=="o/pH2O", fréquence en GHz de la transition (numérateur)
                  nu_den = 2547,                  #si molecule_choisie=="o/pH2O", fréquence en GHz de la transition (dénominateur)
                  H2O_colonne_raie_num = "colonne_TDV", #si molecule_choisie=="o/pH2O", choix de la colonne a saisir dans o/pH2O_int_int.out (numérateur)
                  H2O_colonne_raie_den = "colonne_TDV", #si molecule_choisie=="o/pH2O", choix de la colonne a saisir dans o/pH2O_int_int.out (dénominateur)
                  lambda_microns_num="63m",             #si molecule_choisie=="OI", choix de la colonne a saisir dans intensity.out correspondant à l'atome d'oxygène atomique (numérateur)
                  lambda_microns_den="63m",             #si molecule_choisie=="OI", choix de la colonne a saisir dans intensity.out correspondant à l'atome d'oxygène atomique (dénominateur)
                  vup_num = 0,
                  vup_den = 0,
                  testerVerifieModele=True,
                  printTestsResults = False,
                  typeOfOutput = "default",
                  argumentsToPass = {}):

    ratio = extraire_rapport(racineDir=racine_dir,
                  molecule_choisie_num=molecule_choisie_num,            #choix de la molécule dont est extraite le TDV de la raie (numérateur)
                  molecule_choisie_den=molecule_choisie_den,            #choix de la molécule dont est extraite le TDV de la raie (dénominateur)
                  Jupper_num=Jupper_num,                        ##si molecule_choisie=="CO", J_up de la transition (numérateur)
                  Jupper_den=Jupper_den,                        ##si molecule_choisie=="CO", J_up de la transition (dénominateur)
                  input_dir_num=input_dir,             #choix du modèle cx-ny-bz traité (numérateur)
                  input_dir_den=input_dir,             #choix du modèle cx-ny-bz traité (dénominateur)
                  vitesse_num=vitesse,                       #choix de la vitesse de choc traité dans la nomenclature v"a"-t"b" (numérateur)
                  vitesse_den=vitesse,                       #choix de la vitesse de choc traité dans la nomenclature v"a"-t"b" (dénominateur)
                  age_num=age,                          #choix de l'âge du choc traité dans la nomenclature v"a"-t"b" (numérateur)
                  age_den=age,                          #choix de l'âge du choc traité dans la nomenclature v"a"-t"b" (dénominateur)
                  CO_colonne_raie_num = CO_colonne_raie_num,  #si molecule_choisie=="CO", alors ici est le header de la colonne de raie_dict d'où on on extrait la TDV (numérateur)
                  CO_colonne_raie_den = CO_colonne_raie_den,  #si molecule_choisie=="CO", alors ici est le header de la colonne de raie_dict d'où on on extrait la TDV (dénominateur)
                  nu_num = nu_num,                  #si molecule_choisie=="o/pH2O", fréquence en GHz de la transition (numérateur)
                  nu_den = nu_den,                  #si molecule_choisie=="o/pH2O", fréquence en GHz de la transition (dénominateur)
                  H2O_colonne_raie_num = H2O_colonne_raie_num, #si molecule_choisie=="o/pH2O", choix de la colonne a saisir dans o/pH2O_int_int.out (numérateur)
                  H2O_colonne_raie_den = H2O_colonne_raie_den, #si molecule_choisie=="o/pH2O", choix de la colonne a saisir dans o/pH2O_int_int.out (dénominateur)
                  lambda_microns_num=lambda_microns_num,             #si molecule_choisie=="OI", choix de la colonne a saisir dans intensity.out correspondant à l'atome d'oxygène atomique (numérateur)
                  lambda_microns_den=lambda_microns_den,             #si molecule_choisie=="OI", choix de la colonne a saisir dans intensity.out correspondant à l'atome d'oxygène atomique (dénominateur)
                  vup_num = vup_num,
                  vup_den = vup_den,
                  testerVerifieModele=testerVerifieModele,
                  printTestsResults = printTestsResults,
                  typeOfOutput = typeOfOutput,
                  argumentsToPass = argumentsToPass)

    #axes.scatter(age,ratio['ratio num/den'])

    return ratio



"""=== 21/03/2020 : plot du graphique de gamma_N = f(T en Kelvin) ==="""



os.chdir('/home/mrabenanahary/MPI-AMRVAC/amrvac/src/grackle/src/python/examples')

file1 = ascii.read('Schure2009_cooling_curve_SPEX.dat')
file2 = ascii.read('DelganoMcCray1972_cooling_curve_DM.dat')

XX = [file1[i][0] for i in range(0,len(file1))]
YY = [file1[i][1] for i in range(0,len(file1))]
XXX = [file2[i][0] for i in range(0,len(file2))]
YYY = [file2[i][1] for i in range(0,len(file2))]

#xtab = [file2[i][0] for i in range(0,len(file2))]
xtab = [1e-4,1e-3,1e-2,1e-1]
ytab = [[file2[i][j] for j in range(1,len(file2[i]))] for i in range(0,len(file2))]

a = [list_tangent_slope(xtab,ytab[i]) for i in range(0,len(ytab))]
lambda_hd_over_f_i = np.array(ytab)#-np.log10(np.array(xtab))
lambda_hd_over_f_i = [[lambda_hd_over_f_i[i][j] for i in range(0,len(lambda_hd_over_f_i))] for j in range(0,len(lambda_hd_over_f_i[0]))]
xtab = [file2[i][0] for i in range(0,len(file2))]

pyplot.figure(figsize=(6,6))
"""pyplot.semilogy(XX, 10**np.array(YY),
              color="black",label=r'SPEX $\Lambda$ for solar metallicity above $T>10^4$ K (Schure et al. 2009)',)
"""
a1=[4,3,2,1]
a2=['red','green','blue','orange']
a4=['k','k','red','red','green','green',
    'blue','blue','orange','orange','pink','pink']
a5=['Grackle et al. 2017',
'AMRVAC:Grackle17+',
'Eq. chem. 1 (PyGrackle)',
'AMRVAC:Eq. chem. 1',
'Eq. chem. 2 (PyGrackle)',
'AMRVAC:Eq. chem. 2',
'Non-eq. atom. chem. 1 (PyGrackle)',
'AMRVAC:Non-eq. atom. chem. 1 (PyGrackle)',
'Non-eq. H2 chem. 1 (PyGrackle)',
'AMRVAC:Non-eq. H2 chem. 1 (PyGrackle)',
'Non-eq. H2 chem. 2 (PyGrackle)',
'AMRVAC:Non-eq. H2 chem. 2 (PyGrackle)',]
a6=[False,True,False]
a7=[1,0,2]
a3=[False,False,False,True]

a8=["CloudyData_noUVB.h5","CloudyData_noUVB.h5","CloudyData_noUVB.h5"]
a9=[0.1,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0]
a10=['-','None','--','None','--','None',':','None','dashdot','None','dashdot','None']
a13=['None','x','None','*','None','*','None','P','None','o','None','o']
a14=[0,10,0,10,0,10,0,10,0,10,0,10]
"""for i in range(0,4):
    pyplot.semilogy(xtab,10**np.array(lambda_hd_over_f_i[i]),
                  color=a2[i],label=r'$\Lambda$ for $T<10^4$ K \& $f_i = 10^{-%d}$ (Dalgarno \& McCray 1972)' % (a1[i]),)
"""

filenames = ["Curves - Grackle17+TH.csv",
            "Curves - Grackle17+AMRVAC.csv",
            "Curves - AtomiqueEquil1TH.csv",
            "Curves - AtomiqueEquil1AMRVAC.csv",
            "Curves - AtomiqueEquil2TH.csv",
            "Curves - AtomiqueEquil2AMRVAC.csv",
            "Curves - AtomiqueHorsEquil2TH.csv",
            "Curves - AtomiqueHorsEquil2AMRVAC.csv",
            "Curves - H2HorsEquil1TH.csv",
            "Curves - H2HorsEquil1AMRVAC.csv",
            "Curves - H2HorsEquil2TH.csv",
            "Curves - H2HorsEquil2AMRVAC.csv"]

abreak = [True,True,True,True,True,True,
          False,False,True,True,True,True]

for i,el in enumerate(filenames):  
    if(abreak[i]):continue
    data = ascii.read(el)
    Xdata = list(data[0][:])
    Ydata = list(data[1][:])
    pyplot.semilogy(np.log10(Xdata), Ydata,
                  color=a4[i],lw=a9[i], marker=a13[i],markersize=a14[i],ls=a10[i],label=r'{:s}'.format(a5[i]))
    pyplot.xlabel(r'$\log_{10}(T)$ [K]')
    pyplot.ylabel(r'$|\Lambda|$ [erg s$^{-1}$ cm$^{3}$]')


pyplot.legend(loc='best',
          fancybox=True, shadow=True)    
pyplot.tight_layout()
pyplot.savefig("cooling_curve_check_atom.pdf")
