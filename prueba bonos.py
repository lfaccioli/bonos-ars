#import ffn
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
#import pandas_datareader.data as web
import datetime as dt
#import fix_yahoo_finance as yf

#from pandas_datareader import data, wb
import requests
import urllib.parse
import scipy.optimize as optimize
from scipy.optimize import fsolve
#import eikon as ek
import configparser as cp
from pandas.tseries.offsets import BDay
import streamlit as st


#BAJO LOS DATOS DEL CER
cer_index3 = pd.read_csv("cerindex.csv")

#CODIGOS DE BONOS PARA BAJAR 


#Construyo DATAFRAMES

data_tf = pd.read_csv("datatf.csv")



data_cer = pd.read_csv("datacer.csv")



#Bonos TASA FIJA
to21_b = data_tf.drop(["ARTO233=BA","ARTO263=BA"],axis = 1)
to23_b = data_tf.drop(["ARTO213=BA","ARTO263=BA"],axis = 1)
to26_b = data_tf.drop(["ARTO233=BA","ARTO213=BA"],axis = 1)

###  BONOS CER ####
tx21_b = data_cer.drop(["ARTC213=BA","ART2X13=BA","ART2X23=BA","ARTX223=BA","ARTX233=BA","ARTX243=BA","ARTX263=BA"],axis = 1)
tc21_b = data_cer.drop(["ARTX213=BA","ART2X13=BA","ART2X23=BA","ARTX223=BA","ARTX233=BA","ARTX243=BA","ARTX263=BA"],axis = 1)
t2x1_b = data_cer.drop(["ARTX213=BA","ARTC213=BA","ART2X23=BA","ARTX223=BA","ARTX233=BA","ARTX243=BA","ARTX263=BA"],axis = 1)
t2x2_b = data_cer.drop(["ARTX213=BA","ARTC213=BA","ART2X13=BA","ARTX223=BA","ARTX233=BA","ARTX243=BA","ARTX263=BA"],axis = 1)
tx22_b = data_cer[["ARTX223=BA"]]
tx23_b = data_cer[["ARTX233=BA"]]
tx24_b = data_cer[["ARTX243=BA"]]
tx26_b = data_cer[["ARTX263=BA"]]


#CALCULADORAS

def to21(precio,dia, cupon = 0.1820):
    dia_hoy = dia + dt.timedelta(2)
    #Calculo si t + 2 es sabado o domingo
    if (dt.datetime.weekday(dia_hoy) == 5):
        dia_hoy = dia_hoy + dt.timedelta(2)
    elif (dt.datetime.weekday(dia_hoy) == 6):
        dia_hoy = dia_hoy + dt.timedelta(2)
    else:
        dia_hoy
    
    #Describo las fechas de los cupones del titulo
    c_1 = dt.datetime(2017,4,3)
    c_2 = dt.datetime(2017,10,3)
    c_3 = dt.datetime(2018,4,3)
    c_4 = dt.datetime(2018,10,3)
    c_5 = dt.datetime(2019,4,3)
    c_6 = dt.datetime(2019,10,3)
    c_7 = dt.datetime(2020,4,3)
    c_8 = dt.datetime(2020,10,3)
    c_9 = dt.datetime(2021,4,3)
    c_010 = dt.datetime(2021,10,3)
    
    cupones = [c_1,c_2,c_3,c_4,c_5,c_6,c_7,c_8,c_9,c_010]
    dif_cup = [0,(cupones[1] - cupones[0]).days / 365,(cupones[2] - cupones[1]).days / 365,(cupones[3] - cupones[2]).days / 365,(cupones[4] - cupones[3]).days / 365,(cupones[5] - cupones[4]).days / 365,(cupones[6] - cupones[5]).days / 365,(cupones[7] - cupones[6]).days / 365,(cupones[8] - cupones[7]).days / 365,(cupones[9] - cupones[8]).days / 365]
    dif_fechaact = [(cupones[0] - dia_hoy).days / 365,(cupones[1] - dia_hoy).days / 365,(cupones[2] - dia_hoy).days / 365,(cupones[3] - dia_hoy).days / 365,(cupones[4] - dia_hoy).days / 365,(cupones[5] - dia_hoy).days / 365,(cupones[6] - dia_hoy).days / 365,(cupones[7] - dia_hoy).days / 365,(cupones[8] - dia_hoy).days / 365,(cupones[9] - dia_hoy).days / 365]
    dif_fechas = [(cupones[0] - dia_hoy).days,(cupones[1] - dia_hoy).days,(cupones[2] - dia_hoy).days,(cupones[3] - dia_hoy).days,(cupones[4] - dia_hoy).days,(cupones[5] - dia_hoy).days,(cupones[6] - dia_hoy).days,(cupones[7] - dia_hoy).days,(cupones[8] - dia_hoy).days,(cupones[9] - dia_hoy).days]
    
    VR = 100
    cupon_efectivo = VR * cupon / 2
    cupon_año = VR * cupon
    amortizaciones = [0,0,0,0,0,0,0,0,0,100]
    cashflows = [0,0,0,0,0,0,0,0,0,0]
    
    chequeo_fechas = []
    
    
    datos = pd.DataFrame()
    
    datos["Fechas"] = cupones
    
    for i in cupones:
        if(dia_hoy > i):
            chequeo_fechas.append(1)
        else:
            chequeo_fechas.append(0)
    datos["chequeo"] = chequeo_fechas
    datos["cupon actual"] = np.where(datos["chequeo"].shift(1) == 1,1,0)
    
    datos["fechas desc"] = np.where(datos["chequeo"] == 1,0,dif_fechas)
    datos["dif fecha actual"] = np.where(datos["chequeo"] == 1,0,dif_fechaact)

    datos["VR"] = VR
    #datos["valor tecnico"] = np.round(valor_tecnico,2)
    datos["interes efectivo"] = np.round(cupon_efectivo,2)
    datos["amortizaciones"] = np.round(amortizaciones,2)
    datos["cf cupon actual"] = np.where(datos["cupon actual"] == 1, cupon_año * datos["fechas desc"] / 365,0)
    datos["cashflows"] = cashflows
    datos["cashflows"] = np.where(datos["chequeo"] == 1, 0,datos["interes efectivo"] + datos["amortizaciones"])
    datos["cf final"] = np.where(datos["cf cupon actual"] > 0,datos["cashflows"],datos["cashflows"])
    
    
    
    
    def vpn(tasa, flujos, dias):  
        return np.sum(flujos / (1 + tasa) ** (dias/360))
    
    def tir(flujos, dias, x0):
        return np.asscalar(fsolve(vpn, x0=x0, args=(flujos, dias)))
    
    flujos = datos["cf final"].values
    flujos=np.insert(flujos,0,-precio)
        
    dias = datos["fechas desc"].values
    dias = np.insert(dias,0,0)
    
    ytm = tir(flujos,dias, x0 = 0.01) * 100
    
    datos["TxC"] = datos["fechas desc"] * datos["cashflows"]
    datos["fact"] = ((1 + ytm / 100) ** datos["dif fecha actual"])
    datos["TxC/fact"] = datos["TxC"] / datos["fact"]
    suma_fact = np.sum(datos["TxC/fact"])
    
    dur_dias = suma_fact / precio
    mac_dur = dur_dias / 365
    mod_dur = mac_dur / (1 + ytm / 100)
    valor_tecnico = 100 + np.sum(datos["cf cupon actual"])
    
    pack = {"mac duration": np.round(mac_dur,2),
            "mod duration": np.round(mod_dur,2),
            "ytm": np.round(ytm,2),
            "valor tecnico": np.round(valor_tecnico,2)
        
    } 
    #return datos
    return pack["ytm"]
    #return np.round(ytm,2)
    
    
def to23(precio,dia, cupon = 0.16):
    dia_hoy = dia + dt.timedelta(2)
    #Calculo si t + 2 es sabado o domingo
    if (dt.datetime.weekday(dia_hoy) == 5):
        dia_hoy = dia_hoy + dt.timedelta(2)
    elif (dt.datetime.weekday(dia_hoy) == 6):
        dia_hoy = dia_hoy + dt.timedelta(2)
    else:
        dia_hoy
    
    #Describo las fechas de los cupones del titulo
    c_1 = dt.datetime(2017,4,17)
    c_2 = dt.datetime(2017,10,17)
    c_3 = dt.datetime(2018,4,17)
    c_4 = dt.datetime(2018,10,17)
    c_5 = dt.datetime(2019,4,17)
    c_6 = dt.datetime(2019,10,17)
    c_7 = dt.datetime(2020,4,17)
    c_8 = dt.datetime(2020,10,17)
    c_9 = dt.datetime(2021,4,17)
    c_010 = dt.datetime(2021,10,17)
    c_011 = dt.datetime(2022,4,17)
    c_012 = dt.datetime(2022,10,17)
    c_013 = dt.datetime(2023,4,17)
    c_014 = dt.datetime(2023,10,17)
    
    cupones = [c_1,c_2,c_3,c_4,c_5,c_6,c_7,c_8,c_9,c_010,c_011,c_012,c_013,c_014]
    dif_cup = [0,(cupones[1] - cupones[0]).days / 365,(cupones[2] - cupones[1]).days / 365,(cupones[3] - cupones[2]).days / 365,(cupones[4] - cupones[3]).days / 365,(cupones[5] - cupones[4]).days / 365,(cupones[6] - cupones[5]).days / 365,(cupones[7] - cupones[6]).days / 365,(cupones[8] - cupones[7]).days / 365,(cupones[9] - cupones[8]).days / 365,(cupones[10] - cupones[9]).days / 365,(cupones[11] - cupones[10]).days / 365,(cupones[12] - cupones[11]).days / 365,(cupones[13] - cupones[12]).days / 365]
    dif_fechaact = [(cupones[0] - dia_hoy).days / 365,(cupones[1] - dia_hoy).days / 365,(cupones[2] - dia_hoy).days / 365,(cupones[3] - dia_hoy).days / 365,(cupones[4] - dia_hoy).days / 365,(cupones[5] - dia_hoy).days / 365,(cupones[6] - dia_hoy).days / 365,(cupones[7] - dia_hoy).days / 365,(cupones[8] - dia_hoy).days / 365,(cupones[9] - dia_hoy).days / 365,(cupones[10] - dia_hoy).days / 365,(cupones[11] - dia_hoy).days / 365,(cupones[12] - dia_hoy).days / 365,(cupones[13] - dia_hoy).days / 365]
    dif_fechas = [(cupones[0] - dia_hoy).days,(cupones[1] - dia_hoy).days,(cupones[2] - dia_hoy).days,(cupones[3] - dia_hoy).days,(cupones[4] - dia_hoy).days,(cupones[5] - dia_hoy).days,(cupones[6] - dia_hoy).days,(cupones[7] - dia_hoy).days,(cupones[8] - dia_hoy).days,(cupones[9] - dia_hoy).days,(cupones[10] - dia_hoy).days,(cupones[11] - dia_hoy).days,(cupones[12] - dia_hoy).days,(cupones[13] - dia_hoy).days]
    
    VR = 100
    cupon_efectivo = VR * cupon / 2
    cupon_año = VR * cupon
    amortizaciones = [0,0,0,0,0,0,0,0,0,0,0,0,0,100]
    cashflows = [0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    
    chequeo_fechas = []
    
    
    datos = pd.DataFrame()
    
    datos["Fechas"] = cupones
    
    for i in cupones:
        if(dia_hoy > i):
            chequeo_fechas.append(1)
        else:
            chequeo_fechas.append(0)
    datos["chequeo"] = chequeo_fechas
    datos["cupon actual"] = np.where(datos["chequeo"].shift(1) == 1,1,0)
    
    datos["fechas desc"] = np.where(datos["chequeo"] == 1,0,dif_fechas)
    datos["dif fecha actual"] = np.where(datos["chequeo"] == 1,0,dif_fechaact)

    datos["VR"] = VR
    #datos["valor tecnico"] = np.round(valor_tecnico,2)
    datos["interes efectivo"] = np.round(cupon_efectivo,2)
    datos["amortizaciones"] = np.round(amortizaciones,2)
    datos["cf cupon actual"] = np.where(datos["cupon actual"] == 1, cupon_año * datos["fechas desc"] / 365,0)
    datos["cashflows"] = cashflows
    datos["cashflows"] = np.where(datos["chequeo"] == 1, 0,datos["interes efectivo"] + datos["amortizaciones"])
    datos["cf final"] = np.where(datos["cf cupon actual"] > 0,datos["cashflows"],datos["cashflows"])
    
    
    
    
    def vpn(tasa, flujos, dias):  
        return np.sum(flujos / (1 + tasa) ** (dias/360))
    
    def tir(flujos, dias, x0):
        return np.ndarray.item(fsolve(vpn, x0=x0, args=(flujos, dias)))
    
    flujos = datos["cf final"].values
    flujos=np.insert(flujos,0,-precio)
    
    dias = datos["fechas desc"].values
    dias = np.insert(dias,0,0)
    
    ytm = tir(flujos,dias, x0 = 0.01) * 100
    
    datos["TxC"] = datos["fechas desc"] * datos["cashflows"]
    datos["fact"] = ((1 + ytm / 100) ** datos["dif fecha actual"])
    datos["TxC/fact"] = datos["TxC"] / datos["fact"]
    suma_fact = np.sum(datos["TxC/fact"])
    
    dur_dias = suma_fact / precio
    mac_dur = dur_dias / 365
    mod_dur = mac_dur / (1 + ytm / 100)
    valor_tecnico = 100 + np.sum(datos["cf cupon actual"])
    
    pack = {"mac duration": np.round(mac_dur,2),
            "mod duration": np.round(mod_dur,2),
            "ytm": np.round(ytm,2),
            "valor tecnico": np.round(valor_tecnico,2)
        
    } 
    #return datos
    return pack["ytm"]
    #return np.round(ytm,2)
    
def to26(precio,dia, cupon = 0.155):
    dia_hoy = dia + dt.timedelta(2)
    #Calculo si t + 2 es sabado o domingo
    if (dt.datetime.weekday(dia_hoy) == 5):
        dia_hoy = dia_hoy + dt.timedelta(2)
    elif (dt.datetime.weekday(dia_hoy) == 6):
        dia_hoy = dia_hoy + dt.timedelta(2)
    else:
        dia_hoy
    
    #Describo las fechas de los cupones del titulo
    c_1 = dt.datetime(2017,4,17)
    c_2 = dt.datetime(2017,10,17)
    c_3 = dt.datetime(2018,4,17)
    c_4 = dt.datetime(2018,10,17)
    c_5 = dt.datetime(2019,4,17)
    c_6 = dt.datetime(2019,10,17)
    c_7 = dt.datetime(2020,4,17)
    c_8 = dt.datetime(2020,10,17)
    c_9 = dt.datetime(2021,4,17)
    c_010 = dt.datetime(2021,10,17)
    c_011 = dt.datetime(2022,4,17)
    c_012 = dt.datetime(2022,10,17)
    c_013 = dt.datetime(2023,4,17)
    c_014 = dt.datetime(2023,10,17)
    c_015 = dt.datetime(2024,4,17)
    c_016 = dt.datetime(2024,10,17)
    c_017 = dt.datetime(2025,4,17)
    c_018 = dt.datetime(2025,10,17)
    c_019 = dt.datetime(2026,4,17)
    c_020 = dt.datetime(2026,10,17)
    
    cupones = [c_1,c_2,c_3,c_4,c_5,c_6,c_7,c_8,c_9,c_010,c_011,c_012,c_013,c_014,c_015,c_016,
              c_017,c_018,c_019,c_020]
    dif_cup = [0,(cupones[1] - cupones[0]).days / 365,(cupones[2] - cupones[1]).days / 365,(cupones[3] - cupones[2]).days / 365,
               (cupones[4] - cupones[3]).days / 365,(cupones[5] - cupones[4]).days / 365,(cupones[6] - cupones[5]).days / 365,
               (cupones[7] - cupones[6]).days / 365,(cupones[8] - cupones[7]).days / 365,(cupones[9] - cupones[8]).days / 365,
               (cupones[10] - cupones[9]).days / 365,(cupones[11] - cupones[10]).days / 365,(cupones[12] - cupones[11]).days / 365,
               (cupones[13] - cupones[12]).days / 365,(cupones[14] - cupones[13]).days / 365,(cupones[15] - cupones[14]).days / 365,
              (cupones[16] - cupones[15]).days / 365,(cupones[17] - cupones[16]).days / 365,(cupones[18] - cupones[17]).days / 365,
              (cupones[19] - cupones[18]).days / 365]
    
    
    dif_fechaact = [(cupones[0] - dia_hoy).days / 365,(cupones[1] - dia_hoy).days / 365,(cupones[2] - dia_hoy).days / 365,
                    (cupones[3] - dia_hoy).days / 365,(cupones[4] - dia_hoy).days / 365,(cupones[5] - dia_hoy).days / 365,(cupones[6] - dia_hoy).days / 365,
                    (cupones[7] - dia_hoy).days / 365,(cupones[8] - dia_hoy).days / 365,(cupones[9] - dia_hoy).days / 365,
                    (cupones[10] - dia_hoy).days / 365,(cupones[11] - dia_hoy).days / 365,(cupones[12] - dia_hoy).days / 365,
                    (cupones[13] - dia_hoy).days / 365,(cupones[14] - dia_hoy).days / 365,
                   (cupones[15] - dia_hoy).days / 365,(cupones[16] - dia_hoy).days / 365,
                   (cupones[17] - dia_hoy).days / 365,(cupones[18] - dia_hoy).days / 365,
                   (cupones[19] - dia_hoy).days / 365]
    
    dif_fechas = [(cupones[0] - dia_hoy).days,(cupones[1] - dia_hoy).days,(cupones[2] - dia_hoy).days,
                  (cupones[3] - dia_hoy).days,(cupones[4] - dia_hoy).days,(cupones[5] - dia_hoy).days,
                  (cupones[6] - dia_hoy).days,(cupones[7] - dia_hoy).days,(cupones[8] - dia_hoy).days,
                  (cupones[9] - dia_hoy).days,(cupones[10] - dia_hoy).days,(cupones[11] - dia_hoy).days,
                  (cupones[12] - dia_hoy).days,(cupones[13] - dia_hoy).days,(cupones[14] - dia_hoy).days,
                 (cupones[15] - dia_hoy).days,(cupones[16] - dia_hoy).days,(cupones[17] - dia_hoy).days,
                  (cupones[18] - dia_hoy).days,(cupones[19] - dia_hoy).days
                 
                 ]
    
    VR = 100
    cupon_efectivo = VR * cupon / 2
    cupon_año = VR * cupon
    amortizaciones = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,100]
    cashflows = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    
    chequeo_fechas = []
    
    
    datos = pd.DataFrame()
    
    datos["Fechas"] = cupones
    
    for i in cupones:
        if(dia_hoy > i):
            chequeo_fechas.append(1)
        else:
            chequeo_fechas.append(0)
    datos["chequeo"] = chequeo_fechas
    datos["cupon actual"] = np.where(datos["chequeo"].shift(1) == 1,1,0)
    
    datos["fechas desc"] = np.where(datos["chequeo"] == 1,0,dif_fechas)
    datos["dif fecha actual"] = np.where(datos["chequeo"] == 1,0,dif_fechaact)

    datos["VR"] = VR
    #datos["valor tecnico"] = np.round(valor_tecnico,2)
    datos["interes efectivo"] = np.round(cupon_efectivo,2)
    datos["amortizaciones"] = np.round(amortizaciones,2)
    datos["cf cupon actual"] = np.where(datos["cupon actual"] == 1, cupon_año * datos["fechas desc"] / 365,0)
    datos["cashflows"] = cashflows
    datos["cashflows"] = np.where(datos["chequeo"] == 1, 0,datos["interes efectivo"] + datos["amortizaciones"])
    datos["cf final"] = np.where(datos["cf cupon actual"] > 0,datos["cashflows"],datos["cashflows"])
    
    
    
    
    def vpn(tasa, flujos, dias):  
        return np.sum(flujos / (1 + tasa) ** (dias/360))
    
    def tir(flujos, dias, x0):
        return np.ndarray.item(fsolve(vpn, x0=x0, args=(flujos, dias)))
    
    flujos = datos["cf final"].values
    flujos=np.insert(flujos,0,-precio)
    
    dias = datos["fechas desc"].values
    dias = np.insert(dias,0,0)
    
    ytm = tir(flujos,dias, x0 = 0.01) * 100
    
    datos["TxC"] = datos["fechas desc"] * datos["cashflows"]
    datos["fact"] = ((1 + ytm / 100) ** datos["dif fecha actual"])
    datos["TxC/fact"] = datos["TxC"] / datos["fact"]
    suma_fact = np.sum(datos["TxC/fact"])
    
    dur_dias = suma_fact / precio
    mac_dur = dur_dias / 365
    mod_dur = mac_dur / (1 + ytm / 100)
    valor_tecnico = 100 + np.sum(datos["cf cupon actual"])
    
    pack = {"mac duration": np.round(mac_dur,2),
            "mod duration": np.round(mod_dur,2),
            "ytm": np.round(ytm,2),
            "valor tecnico": np.round(valor_tecnico,2)
        
    } 
    #return datos
    return pack["ytm"]
    #return np.round(ytm,2)
    
def tx21(precio,dia,cer_hoy,cupon = 0.01, cer_emision = 19.2205):
    
    dia_hoy = dia + dt.timedelta(2)
    #Calculo si t + 2 es sabado o domingo
    if (dt.datetime.weekday(dia_hoy) == 5):
        dia_hoy = dia_hoy + dt.timedelta(2)
    elif (dt.datetime.weekday(dia_hoy) == 6):
        dia_hoy = dia_hoy + dt.timedelta(2)
    else:
        dia_hoy
    
    #Describo las fechas de los cupones del titulo
    c_1 = dt.datetime(2020,2,5)
    c_2 = dt.datetime(2020,8,5)
    c_3 = dt.datetime(2021,2,5)
    c_4 = dt.datetime(2021,8,5)
    
    cupones = [c_1,c_2,c_3,c_4]
    dif_cup = [0,(cupones[1] - cupones[0]).days / 365,(cupones[2] - cupones[1]).days / 365,(cupones[3] - cupones[2]).days / 365]
    dif_fechaact = [(cupones[0] - dia_hoy).days / 365,(cupones[1] - dia_hoy).days / 365,(cupones[2] - dia_hoy).days / 365,(cupones[3] - dia_hoy).days / 365]
    dif_fechas = [(cupones[0] - dia_hoy).days,(cupones[1] - dia_hoy).days,(cupones[2] - dia_hoy).days,(cupones[3] - dia_hoy).days]
    
    
    indexacion = cer_hoy / cer_emision
    
    valor_residual = 100
    valor_tecnico = valor_residual * indexacion
    cupon_efectivo = cupon * valor_tecnico / 2
    amort = valor_tecnico
    amortizaciones = [0,0,0,valor_tecnico]
    cashflows = [0,0,0,0]
    
    chequeo_fechas = []
    
    
    datos = pd.DataFrame()
    
    datos["Fechas"] = cupones
    
    for i in cupones:
        if(dia_hoy > i):
            chequeo_fechas.append(1)
        else:
            chequeo_fechas.append(0)
    datos["chequeo"] = chequeo_fechas
    #datos["cupon actual"] = np.where(datos["chequeo"].shift(1) == 1,1,0)
    
    #datos["dif fechas"] = dif_cup
    datos["fechas desc"] = np.where(datos["chequeo"] == 1,0,dif_fechas)
    datos["dif fecha actual"] = np.where(datos["chequeo"] == 1,0,dif_fechaact)
    

    
    
    
    
    
    datos["VR"] = valor_residual
    datos["valor tecnico"] = np.round(valor_tecnico,2)
    datos["interes efectivo"] = np.round(cupon_efectivo,2)
    datos["amortizaciones"] = np.round(amortizaciones,2)
    datos["cashflows"] = cashflows
    
    datos["cashflows"] = np.where(datos["chequeo"] == 1, 0,datos["interes efectivo"] + datos["amortizaciones"])
    
    
    def vpn(tasa, flujos, dias):  
        return np.sum(flujos / (1 + tasa) ** (dias/360))
    
    def tir(flujos, dias, x0):
        return np.ndarray.item(fsolve(vpn, x0=x0, args=(flujos, dias)))
    
    flujos = datos["cashflows"].values
    flujos=np.insert(flujos,0,-precio)
    
    dias = datos["fechas desc"].values
    dias = np.insert(dias,0,0)
    
    ytm = tir(flujos,dias, x0 = 0.01) * 100
    
    datos["TxC"] = datos["fechas desc"] * datos["cashflows"]
    datos["fact"] = ((1 + ytm / 100) ** datos["dif fecha actual"])
    datos["TxC/fact"] = datos["TxC"] / datos["fact"]
    suma_fact = np.sum(datos["TxC/fact"])
    
    dur_dias = suma_fact / precio
    mac_dur = dur_dias / 365
    mod_dur = mac_dur / (1 + ytm / 100)
    
    pack = {"mac duration": np.round(mac_dur,2),
            "mod duration": np.round(mod_dur,2),
            "ytm": np.round(ytm,2),
            "valor tecnico": np.round(valor_tecnico,2),
            "vencimiento": c_4
        
    } 
    #return datos
    return pack["ytm"]
    #return np.round(ytm,2)
    
    
def tx22(precio,dia,cer_hoy,cupon = 0.012, cer_emision = 20.0733):
    
    dia_hoy = dt.datetime.today() + dt.timedelta(2)
    #Calculo si t + 2 es sabado o domingo
    if (dt.datetime.weekday(dia_hoy) == 5):
        dia_hoy = dia_hoy + dt.timedelta(2)
    elif (dt.datetime.weekday(dia_hoy) == 6):
        dia_hoy = dia_hoy + dt.timedelta(2)
    else:
        dia_hoy
    
    #Describo las fechas de los cupones del titulo
    c_1 = dt.datetime(2020,3,18)
    c_2 = dt.datetime(2020,9,18)
    c_3 = dt.datetime(2021,3,18)
    c_4 = dt.datetime(2021,9,18)
    c_5 = dt.datetime(2022,3,18)
    
    cupones = [c_1,c_2,c_3,c_4,c_5]
    dif_cup = [0,(cupones[1] - cupones[0]).days / 365,(cupones[2] - cupones[1]).days / 365,(cupones[3] - cupones[2]).days / 365,(cupones[4] - cupones[3]).days / 365]
    dif_fechaact = [(cupones[0] - dia_hoy).days / 365,(cupones[1] - dia_hoy).days / 365,(cupones[2] - dia_hoy).days / 365,(cupones[3] - dia_hoy).days / 365,(cupones[4] - dia_hoy).days / 365]
    dif_fechas = [(cupones[0] - dia_hoy).days,(cupones[1] - dia_hoy).days,(cupones[2] - dia_hoy).days,(cupones[3] - dia_hoy).days,(cupones[4] - dia_hoy).days]
    
    
    indexacion = cer_hoy / cer_emision
    
    valor_residual = 100
    valor_tecnico = valor_residual * indexacion
    cupon_efectivo = cupon * valor_tecnico / 2
    amort = valor_tecnico
    amortizaciones = [0,0,0,0,valor_tecnico]
    cashflows = [0,0,0,0,0]
    
    chequeo_fechas = []
    
    
    datos = pd.DataFrame()
    
    datos["Fechas"] = cupones
    
    for i in cupones:
        if(dia_hoy > i):
            chequeo_fechas.append(1)
        else:
            chequeo_fechas.append(0)
    datos["chequeo"] = chequeo_fechas
    #datos["cupon actual"] = np.where(datos["chequeo"].shift(1) == 1,1,0)
    
    #datos["dif fechas"] = dif_cup
    datos["fechas desc"] = np.where(datos["chequeo"] == 1,0,dif_fechas)
    datos["dif fecha actual"] = np.where(datos["chequeo"] == 1,0,dif_fechaact)
    

    

    datos["VR"] = valor_residual
    datos["valor tecnico"] = np.round(valor_tecnico,2)
    datos["interes efectivo"] = np.round(cupon_efectivo,2)
    datos["amortizaciones"] = np.round(amortizaciones,2)
    datos["cashflows"] = cashflows
    
    datos["cashflows"] = np.where(datos["chequeo"] == 1, 0,datos["interes efectivo"] + datos["amortizaciones"])
    
    
    def vpn(tasa, flujos, dias):  
        return np.sum(flujos / (1 + tasa) ** (dias/360))
    
    def tir(flujos, dias, x0):
        return np.ndarray.item(fsolve(vpn, x0=x0, args=(flujos, dias)))
    
    flujos = datos["cashflows"].values
    flujos=np.insert(flujos,0,-precio)
    
    dias = datos["fechas desc"].values
    dias = np.insert(dias,0,0)
    
    ytm = tir(flujos,dias, x0 = 0.01) * 100
    
    datos["TxC"] = datos["fechas desc"] * datos["cashflows"]
    datos["fact"] = ((1 + ytm / 100) ** datos["dif fecha actual"])
    datos["TxC/fact"] = datos["TxC"] / datos["fact"]
    suma_fact = np.sum(datos["TxC/fact"])
    
    dur_dias = suma_fact / precio
    mac_dur = dur_dias / 365
    mod_dur = mac_dur / (1 + ytm / 100)
    
    pack = {"mac duration": np.round(mac_dur,2),
            "mod duration": np.round(mod_dur,2),
            "ytm": np.round(ytm,2),
            "valor tecnico": np.round(valor_tecnico,2)
        
    } 
    
    return pack["ytm"]


def t2x2(precio,dia,cer_hoy,cupon = 0.013, cer_emision = 21.1268552445605 ):
    dia_hoy = dia + dt.timedelta(2)
    #Calculo si t + 2 es sabado o domingo
    if (dt.datetime.weekday(dia_hoy) == 5):
        dia_hoy = dia_hoy + dt.timedelta(2)
    elif (dt.datetime.weekday(dia_hoy) == 6):
        dia_hoy = dia_hoy + dt.timedelta(2)
    else:
        dia_hoy
    
    #Describo las fechas de los cupones del titulo
    c_1 = dt.datetime(2020,9,20)
    c_2 = dt.datetime(2021,3,20)
    c_3 = dt.datetime(2021,9,20)
    c_4 = dt.datetime(2022,3,20)
    c_5 = dt.datetime(2022,9,20)
    
    cupones = [c_1,c_2,c_3,c_4,c_5]
    dif_cup = [0,(cupones[1] - cupones[0]).days / 365,(cupones[2] - cupones[1]).days / 365,(cupones[3] - cupones[2]).days / 365,(cupones[4] - cupones[3]).days / 365]
    dif_fechaact = [(cupones[0] - dia_hoy).days / 365,(cupones[1] - dia_hoy).days / 365,(cupones[2] - dia_hoy).days / 365,(cupones[3] - dia_hoy).days / 365,(cupones[4] - dia_hoy).days / 365]
    dif_fechas = [(cupones[0] - dia_hoy).days,(cupones[1] - dia_hoy).days,(cupones[2] - dia_hoy).days,(cupones[3] - dia_hoy).days,(cupones[4] - dia_hoy).days]
    
    
    indexacion = cer_hoy / cer_emision
    
    valor_residual = 100
    valor_tecnico = valor_residual * indexacion
    cupon_efectivo = cupon * valor_tecnico / 2
    amort = valor_tecnico
    amortizaciones = [0,0,0,0,valor_tecnico]
    cashflows = [0,0,0,0,0]
    
    chequeo_fechas = []
    
    
    datos = pd.DataFrame()
    
    datos["Fechas"] = cupones
    
    for i in cupones:
        if(dia_hoy > i):
            chequeo_fechas.append(1)
        else:
            chequeo_fechas.append(0)
    datos["chequeo"] = chequeo_fechas
    #datos["cupon actual"] = np.where(datos["chequeo"].shift(1) == 1,1,0)
    
    #datos["dif fechas"] = dif_cup
    datos["fechas desc"] = np.where(datos["chequeo"] == 1,0,dif_fechas)
    datos["dif fecha actual"] = np.where(datos["chequeo"] == 1,0,dif_fechaact)
    

    

    datos["VR"] = valor_residual
    datos["valor tecnico"] = np.round(valor_tecnico,2)
    datos["interes efectivo"] = np.round(cupon_efectivo,2)
    datos["amortizaciones"] = np.round(amortizaciones,2)
    datos["cashflows"] = cashflows
    
    datos["cashflows"] = np.where(datos["chequeo"] == 1, 0,datos["interes efectivo"] + datos["amortizaciones"])
    
    
    def vpn(tasa, flujos, dias):  
        return np.sum(flujos / (1 + tasa) ** (dias/360))
    
    def tir(flujos, dias, x0):
        return np.ndarray.item(fsolve(vpn, x0=x0, args=(flujos, dias)))
    
    flujos = datos["cashflows"].values
    flujos=np.insert(flujos,0,-precio)
    
    dias = datos["fechas desc"].values
    dias = np.insert(dias,0,0)
    
    ytm = tir(flujos,dias, x0 = 0.01) * 100
    
    datos["TxC"] = datos["fechas desc"] * datos["cashflows"]
    datos["fact"] = ((1 + ytm / 100) ** datos["dif fecha actual"])
    datos["TxC/fact"] = datos["TxC"] / datos["fact"]
    suma_fact = np.sum(datos["TxC/fact"])
    
    dur_dias = suma_fact / precio
    mac_dur = dur_dias / 365
    mod_dur = mac_dur / (1 + ytm / 100)
    
    pack = {"mac duration": np.round(mac_dur,2),
            "mod duration": np.round(mod_dur,2),
            "ytm": np.round(ytm,2),
            "valor tecnico": np.round(valor_tecnico,2)
        
    } 
    
    return pack["ytm"]


def tx23(precio,dia,cer_hoy,cupon = 0.014, cer_emision = 20.1521):
    
    dia_hoy = dia + dt.timedelta(2)
    #Calculo si t + 2 es sabado o domingo
    if (dt.datetime.weekday(dia_hoy) == 5):
        dia_hoy = dia_hoy + dt.timedelta(2)
    elif (dt.datetime.weekday(dia_hoy) == 6):
        dia_hoy = dia_hoy + dt.timedelta(2)
    else:
        dia_hoy
    
    #Describo las fechas de los cupones del titulo
    c_1 = dt.datetime(2020,3,25)
    c_2 = dt.datetime(2020,9,25)
    c_3 = dt.datetime(2021,3,25)
    c_4 = dt.datetime(2021,9,25)
    c_5 = dt.datetime(2022,3,25)
    c_6 = dt.datetime(2022,9,25)
    c_7 = dt.datetime(2023,3,25)
    
    cupones = [c_1,c_2,c_3,c_4,c_5,c_6,c_7]
    dif_cup = [0,(cupones[1] - cupones[0]).days / 365,(cupones[2] - cupones[1]).days / 365,(cupones[3] - cupones[2]).days / 365,(cupones[4] - cupones[3]).days / 365,(cupones[5] - cupones[4]).days / 365,(cupones[6] - cupones[5]).days / 365]
    dif_fechaact = [(cupones[0] - dia_hoy).days / 365,(cupones[1] - dia_hoy).days / 365,(cupones[2] - dia_hoy).days / 365,(cupones[3] - dia_hoy).days / 365,(cupones[4] - dia_hoy).days / 365,(cupones[5] - dia_hoy).days / 365,(cupones[6] - dia_hoy).days / 365]
    dif_fechas = [(cupones[0] - dia_hoy).days,(cupones[1] - dia_hoy).days,(cupones[2] - dia_hoy).days,(cupones[3] - dia_hoy).days,(cupones[4] - dia_hoy).days,(cupones[5] - dia_hoy).days,(cupones[6] - dia_hoy).days]
    
    
    indexacion = cer_hoy / cer_emision
    
    valor_residual = 100
    valor_tecnico = valor_residual * indexacion
    cupon_efectivo = cupon * valor_tecnico / 2
    amort = valor_tecnico
    amortizaciones = [0,0,0,0,0,0,valor_tecnico]
    cashflows = [0,0,0,0,0,0,0]
    
    chequeo_fechas = []
    
    
    datos = pd.DataFrame()
    
    datos["Fechas"] = cupones
    
    for i in cupones:
        if(dia_hoy > i):
            chequeo_fechas.append(1)
        else:
            chequeo_fechas.append(0)
    datos["chequeo"] = chequeo_fechas
    #datos["cupon actual"] = np.where(datos["chequeo"].shift(1) == 1,1,0)
    
    #datos["dif fechas"] = dif_cup
    datos["fechas desc"] = np.where(datos["chequeo"] == 1,0,dif_fechas)
    datos["dif fecha actual"] = np.where(datos["chequeo"] == 1,0,dif_fechaact)
    

    

    datos["VR"] = valor_residual
    datos["valor tecnico"] = np.round(valor_tecnico,2)
    datos["interes efectivo"] = np.round(cupon_efectivo,2)
    datos["amortizaciones"] = np.round(amortizaciones,2)
    datos["cashflows"] = cashflows
    
    datos["cashflows"] = np.where(datos["chequeo"] == 1, 0,datos["interes efectivo"] + datos["amortizaciones"])
    
    
    def vpn(tasa, flujos, dias):  
        return np.sum(flujos / (1 + tasa) ** (dias/360))
    
    def tir(flujos, dias, x0):
        return np.ndarray.item(fsolve(vpn, x0=x0, args=(flujos, dias)))
    
    flujos = datos["cashflows"].values
    flujos=np.insert(flujos,0,-precio)
    
    dias = datos["fechas desc"].values
    dias = np.insert(dias,0,0)
    
    ytm = tir(flujos,dias, x0 = 0.01) * 100
    
    datos["TxC"] = datos["fechas desc"] * datos["cashflows"]
    datos["fact"] = ((1 + ytm / 100) ** datos["dif fecha actual"])
    datos["TxC/fact"] = datos["TxC"] / datos["fact"]
    suma_fact = np.sum(datos["TxC/fact"])
    
    dur_dias = suma_fact / precio
    mac_dur = dur_dias / 365
    mod_dur = mac_dur / (1 + ytm / 100)
    
    pack = {"mac duration": np.round(mac_dur,2),
            "mod duration": np.round(mod_dur,2),
            "ytm": np.round(ytm,2),
            "valor tecnico": np.round(valor_tecnico,2)
        
    } 
    
    return pack["ytm"]
    
    #return datos,ytm
    
def tx24(precio,dia,cer_hoy,cupon = 0.015, cer_emision = 20.1521):
    
    dia_hoy = dia + dt.timedelta(2)
    #Calculo si t + 2 es sabado o domingo
    if (dt.datetime.weekday(dia_hoy) == 5):
        dia_hoy = dia_hoy + dt.timedelta(2)
    elif (dt.datetime.weekday(dia_hoy) == 6):
        dia_hoy = dia_hoy + dt.timedelta(2)
    else:
        dia_hoy
    
    #Describo las fechas de los cupones del titulo
    c_1 = dt.datetime(2020,3,25)
    c_2 = dt.datetime(2020,9,25)
    c_3 = dt.datetime(2021,3,25)
    c_4 = dt.datetime(2021,9,25)
    c_5 = dt.datetime(2022,3,25)
    c_6 = dt.datetime(2022,9,25)
    c_7 = dt.datetime(2023,3,25)
    c_8 = dt.datetime(2023,9,25)
    c_9 = dt.datetime(2024,3,25)
    
    cupones = [c_1,c_2,c_3,c_4,c_5,c_6,c_7,c_8,c_9]
    dif_cup = [0,(cupones[1] - cupones[0]).days / 365,(cupones[2] - cupones[1]).days / 365,(cupones[3] - cupones[2]).days / 365,(cupones[4] - cupones[3]).days / 365,(cupones[5] - cupones[4]).days / 365,(cupones[6] - cupones[5]).days / 365,(cupones[7] - cupones[6]).days / 365,(cupones[8] - cupones[7]).days / 365]
    dif_fechaact = [(cupones[0] - dia_hoy).days / 365,(cupones[1] - dia_hoy).days / 365,(cupones[2] - dia_hoy).days / 365,(cupones[3] - dia_hoy).days / 365,(cupones[4] - dia_hoy).days / 365,(cupones[5] - dia_hoy).days / 365,(cupones[6] - dia_hoy).days / 365,(cupones[7] - dia_hoy).days / 365,(cupones[8] - dia_hoy).days / 365]
    dif_fechas = [(cupones[0] - dia_hoy).days,(cupones[1] - dia_hoy).days,(cupones[2] - dia_hoy).days,(cupones[3] - dia_hoy).days,(cupones[4] - dia_hoy).days,(cupones[5] - dia_hoy).days,(cupones[6] - dia_hoy).days,(cupones[7] - dia_hoy).days,(cupones[8] - dia_hoy).days]
    
    
    indexacion = cer_hoy / cer_emision
    
    valor_residual = 100
    valor_tecnico = valor_residual * indexacion
    cupon_efectivo = cupon * valor_tecnico / 2
    amort = valor_tecnico
    amortizaciones = [0,0,0,0,0,0,0,0,valor_tecnico]
    cashflows = [0,0,0,0,0,0,0,0,0]
    
    chequeo_fechas = []
    
    
    datos = pd.DataFrame()
    
    datos["Fechas"] = cupones
    
    for i in cupones:
        if(dia_hoy > i):
            chequeo_fechas.append(1)
        else:
            chequeo_fechas.append(0)
    datos["chequeo"] = chequeo_fechas
    #datos["cupon actual"] = np.where(datos["chequeo"].shift(1) == 1,1,0)
    
    #datos["dif fechas"] = dif_cup
    datos["fechas desc"] = np.where(datos["chequeo"] == 1,0,dif_fechas)
    datos["dif fecha actual"] = np.where(datos["chequeo"] == 1,0,dif_fechaact)
    

    datos["VR"] = valor_residual
    datos["valor tecnico"] = np.round(valor_tecnico,2)
    datos["interes efectivo"] = np.round(cupon_efectivo,2)
    datos["amortizaciones"] = np.round(amortizaciones,2)
    datos["cashflows"] = cashflows
    
    datos["cashflows"] = np.where(datos["chequeo"] == 1, 0,datos["interes efectivo"] + datos["amortizaciones"])
    
    
    def vpn(tasa, flujos, dias):  
        return np.sum(flujos / (1 + tasa) ** (dias/360))
    
    def tir(flujos, dias, x0):
        return np.ndarray.item(fsolve(vpn, x0=x0, args=(flujos, dias)))
    
    flujos = datos["cashflows"].values
    flujos=np.insert(flujos,0,-precio)
    
    dias = datos["fechas desc"].values
    dias = np.insert(dias,0,0)
    
    ytm = tir(flujos,dias, x0 = 0.01) * 100
    
    datos["TxC"] = datos["fechas desc"] * datos["cashflows"]
    datos["fact"] = ((1 + ytm / 100) ** datos["dif fecha actual"])
    datos["TxC/fact"] = datos["TxC"] / datos["fact"]
    suma_fact = np.sum(datos["TxC/fact"])
    
    dur_dias = suma_fact / precio
    mac_dur = dur_dias / 365
    mod_dur = mac_dur / (1 + ytm / 100)
    
    pack = {"mac duration": np.round(mac_dur,2),
            "mod duration": np.round(mod_dur,2),
            "ytm": np.round(ytm,2),
            "valor tecnico": np.round(valor_tecnico,2)
        
    } 
    
    return pack["ytm"]
    #return datos,ytm
    
    
def tx26(precio,dia,cer_hoy,cupon = 0.02, cer_emision = 22.5303):
    
    dia_hoy = dia + dt.timedelta(2)
    #Calculo si t + 2 es sabado o domingo
    if (dt.datetime.weekday(dia_hoy) == 5):
        dia_hoy = dia_hoy + dt.timedelta(2)
    elif (dt.datetime.weekday(dia_hoy) == 6):
        dia_hoy = dia_hoy + dt.timedelta(2)
    else:
        dia_hoy
    
    #Describo las fechas de los cupones del titulo
    c_1 = dt.datetime(2021,5,4)
    c_2 = dt.datetime(2021,11,4)
    c_3 = dt.datetime(2022,5,4)
    c_4 = dt.datetime(2022,11,4)
    c_5 = dt.datetime(2023,5,4)
    c_6 = dt.datetime(2022,11,4)
    c_7 = dt.datetime(2023,5,4)
    c_8 = dt.datetime(2023,11,4)
    c_9 = dt.datetime(2024,5,4)
    
    c_10 = dt.datetime(2024,11,4)
    c_11 = dt.datetime(2025,5,4)
    c_12 = dt.datetime(2025,11,4)
    c_13 = dt.datetime(2026,5,4)
    c_14 = dt.datetime(2026,11,4)
    
    
    cupones = [c_1,c_2,c_3,c_4,c_5,c_6,c_7,c_8,c_9,c_10,c_11,c_12,c_13,c_14]
    dif_cup = [0,(cupones[1] - cupones[0]).days / 365,(cupones[2] - cupones[1]).days / 365,(cupones[3] - cupones[2]).days / 365,(cupones[4] - cupones[3]).days / 365,(cupones[5] - cupones[4]).days / 365,(cupones[6] - cupones[5]).days / 365,(cupones[7] - cupones[6]).days / 365,
               
               
            (cupones[8] - cupones[7]).days / 365,(cupones[9] - cupones[8]).days / 365,
              (cupones[10] - cupones[9]).days / 365,(cupones[11] - cupones[10]).days / 365,
              (cupones[12] - cupones[11]).days / 365,(cupones[13] - cupones[12]).days / 365]
    
    dif_fechaact = [(cupones[0] - dia_hoy).days / 365,(cupones[1] - dia_hoy).days / 365,(cupones[2] - dia_hoy).days / 365,(cupones[3] - dia_hoy).days / 365,(cupones[4] - dia_hoy).days / 365,(cupones[5] - dia_hoy).days / 365,(cupones[6] - dia_hoy).days / 365,(cupones[7] - dia_hoy).days / 365,(cupones[8] - dia_hoy).days / 365,
                   
                   (cupones[9] - dia_hoy).days / 365,(cupones[10] - dia_hoy).days / 365,(cupones[11] - dia_hoy).days / 365,
                   (cupones[12] - dia_hoy).days / 365,(cupones[13] - dia_hoy).days / 365]
    dif_fechas = [(cupones[0] - dia_hoy).days,(cupones[1] - dia_hoy).days,(cupones[2] - dia_hoy).days,(cupones[3] - dia_hoy).days,(cupones[4] - dia_hoy).days,(cupones[5] - dia_hoy).days,(cupones[6] - dia_hoy).days,(cupones[7] - dia_hoy).days,(cupones[8] - dia_hoy).days,
                 (cupones[9] - dia_hoy).days,(cupones[10] - dia_hoy).days,(cupones[11] - dia_hoy).days,(cupones[12] - dia_hoy).days,(cupones[13] - dia_hoy).days,
                 ]
    
    
    indexacion = cer_hoy / cer_emision
    
    valor_residual = 100
    valor_tecnico = valor_residual * indexacion
    cupon_efectivo = cupon * valor_tecnico / 2
    quita = valor_tecnico * 0.25
    amort = valor_tecnico
    amortizaciones = [0,0,0,0,0,0,0,0,0,0,valor_tecnico * 0.25,valor_tecnico * 0.25,valor_tecnico * 0.25,valor_tecnico * 0.25]
    cashflows = [0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    val_tec = [valor_tecnico,valor_tecnico,valor_tecnico,valor_tecnico,valor_tecnico,valor_tecnico,valor_tecnico,valor_tecnico,
              valor_tecnico,valor_tecnico,valor_tecnico,valor_tecnico - quita,valor_tecnico - quita * 2,valor_tecnico - quita * 3]
    
    int_nuevo = [cupon_efectivo,cupon_efectivo,cupon_efectivo,cupon_efectivo,cupon_efectivo,
                 cupon_efectivo,cupon_efectivo,cupon_efectivo,cupon_efectivo,cupon_efectivo,
                 cupon_efectivo,(valor_tecnico - quita) * cupon /2,(valor_tecnico - quita * 2) * cupon /2,
                 (valor_tecnico - quita * 3) * cupon /2
        
        
    ]
    
    chequeo_fechas = []
    
    
    datos = pd.DataFrame()
    
    datos["Fechas"] = cupones
    
    for i in cupones:
        if(dia_hoy > i):
            chequeo_fechas.append(1)
        else:
            chequeo_fechas.append(0)
    datos["chequeo"] = chequeo_fechas
    #datos["cupon actual"] = np.where(datos["chequeo"].shift(1) == 1,1,0)
    
    #datos["dif fechas"] = dif_cup
    datos["fechas desc"] = np.where(datos["chequeo"] == 1,0,dif_fechas)
    datos["dif fecha actual"] = np.where(datos["chequeo"] == 1,0,dif_fechaact)
    

    datos["VR"] = valor_residual
    datos["valor tecnico"] = np.round(val_tec,2)
    datos["interes efectivo"] = np.round(int_nuevo,2)
    datos["amortizaciones"] = np.round(amortizaciones,2)
    datos["cashflows"] = cashflows
    
    datos["cashflows"] = np.where(datos["chequeo"] == 1, 0,datos["interes efectivo"] + datos["amortizaciones"])
    
    
    def vpn(tasa, flujos, dias):  
        return np.sum(flujos / (1 + tasa) ** (dias/360))
    
    def tir(flujos, dias, x0):
        return np.ndarray.item(fsolve(vpn, x0=x0, args=(flujos, dias)))
    
    flujos = datos["cashflows"].values
    flujos=np.insert(flujos,0,-precio)
    
    dias = datos["fechas desc"].values
    dias = np.insert(dias,0,0)
    
    ytm = tir(flujos,dias, x0 = 0.01) * 100
    
    datos["TxC"] = datos["fechas desc"] * datos["cashflows"]
    datos["fact"] = ((1 + ytm / 100) ** datos["dif fecha actual"])
    datos["TxC/fact"] = datos["TxC"] / datos["fact"]
    suma_fact = np.sum(datos["TxC/fact"])
    
    dur_dias = suma_fact / precio
    mac_dur = dur_dias / 365
    mod_dur = mac_dur / (1 + ytm / 100)
    
    pack = {"mac duration": np.round(mac_dur,2),
            "mod duration": np.round(mod_dur,2),
            "ytm": np.round(ytm,2),
            "valor tecnico": np.round(valor_tecnico,2)
        
    } 
    
    return pack["ytm"]
    #return datos,ytm
    
    
def tc21(precio,dia,cer_hoy,cupon = 0.025, cer_emision = 6.1534):
    
    dia_hoy = dia + dt.timedelta(2)
    #Calculo si t + 2 es sabado o domingo
    if (dt.datetime.weekday(dia_hoy) == 5):
        dia_hoy = dia_hoy + dt.timedelta(2)
    elif (dt.datetime.weekday(dia_hoy) == 6):
        dia_hoy = dia_hoy + dt.timedelta(2)
    else:
        dia_hoy
    
    #Describo las fechas de los cupones del titulo
    c_1 = dt.datetime(2017,1,22)
    c_2 = dt.datetime(2017,7,22)
    c_3 = dt.datetime(2018,1,22)
    c_4 = dt.datetime(2018,7,22)
    c_5 = dt.datetime(2019,1,22)
    c_6 = dt.datetime(2019,7,22)
    c_7 = dt.datetime(2020,1,22)
    c_8 = dt.datetime(2020,7,22)
    c_9 = dt.datetime(2021,1,22)
    c_010 = dt.datetime(2021,7,22)
    
    cupones = [c_1,c_2,c_3,c_4,c_5,c_6,c_7,c_8,c_9,c_010]
    dif_cup = [0,(cupones[1] - cupones[0]).days / 365,(cupones[2] - cupones[1]).days / 365,(cupones[3] - cupones[2]).days / 365,(cupones[4] - cupones[3]).days / 365,(cupones[5] - cupones[4]).days / 365,(cupones[6] - cupones[5]).days / 365,(cupones[7] - cupones[6]).days / 365,(cupones[8] - cupones[7]).days / 365,(cupones[9] - cupones[8]).days / 365]
    dif_fechaact = [(cupones[0] - dia_hoy).days / 365,(cupones[1] - dia_hoy).days / 365,(cupones[2] - dia_hoy).days / 365,(cupones[3] - dia_hoy).days / 365,(cupones[4] - dia_hoy).days / 365,(cupones[5] - dia_hoy).days / 365,(cupones[6] - dia_hoy).days / 365,(cupones[7] - dia_hoy).days / 365,(cupones[8] - dia_hoy).days / 365,(cupones[9] - dia_hoy).days / 365]
    dif_fechas = [(cupones[0] - dia_hoy).days,(cupones[1] - dia_hoy).days,(cupones[2] - dia_hoy).days,(cupones[3] - dia_hoy).days,(cupones[4] - dia_hoy).days,(cupones[5] - dia_hoy).days,(cupones[6] - dia_hoy).days,(cupones[7] - dia_hoy).days,(cupones[8] - dia_hoy).days,(cupones[9] - dia_hoy).days]
    
    
    indexacion = cer_hoy / cer_emision
    
    valor_residual = 100
    valor_tecnico = valor_residual * indexacion
    cupon_efectivo = cupon * valor_tecnico / 2
    amort = valor_tecnico
    amortizaciones = [0,0,0,0,0,0,0,0,0,valor_tecnico]
    cashflows = [0,0,0,0,0,0,0,0,0,0]
    
    chequeo_fechas = []
    
    
    datos = pd.DataFrame()
    
    datos["Fechas"] = cupones
    
    for i in cupones:
        if(dia_hoy > i):
            chequeo_fechas.append(1)
        else:
            chequeo_fechas.append(0)
    datos["chequeo"] = chequeo_fechas
    #datos["cupon actual"] = np.where(datos["chequeo"].shift(1) == 1,1,0)
    
    #datos["dif fechas"] = dif_cup
    datos["fechas desc"] = np.where(datos["chequeo"] == 1,0,dif_fechas)
    datos["dif fecha actual"] = np.where(datos["chequeo"] == 1,0,dif_fechaact)
    

    datos["VR"] = valor_residual
    datos["valor tecnico"] = np.round(valor_tecnico,2)
    datos["interes efectivo"] = np.round(cupon_efectivo,2)
    datos["amortizaciones"] = np.round(amortizaciones,2)
    datos["cashflows"] = cashflows
    
    datos["cashflows"] = np.where(datos["chequeo"] == 1, 0,datos["interes efectivo"] + datos["amortizaciones"])
    
    
    def vpn(tasa, flujos, dias):  
        return np.sum(flujos / (1 + tasa) ** (dias/360))
    
    def tir(flujos, dias, x0):
        return np.ndarray.item(fsolve(vpn, x0=x0, args=(flujos, dias)))
    
    flujos = datos["cashflows"].values
    flujos=np.insert(flujos,0,-precio)
    
    dias = datos["fechas desc"].values
    dias = np.insert(dias,0,0)
    
    ytm = tir(flujos,dias, x0 = 0.01) * 100
    
    datos["TxC"] = datos["fechas desc"] * datos["cashflows"]
    datos["fact"] = ((1 + ytm / 100) ** datos["dif fecha actual"])
    datos["TxC/fact"] = datos["TxC"] / datos["fact"]
    suma_fact = np.sum(datos["TxC/fact"])
    
    dur_dias = suma_fact / precio
    mac_dur = dur_dias / 365
    mod_dur = mac_dur / (1 + ytm / 100)
    
    pack = {"mac duration": np.round(mac_dur,2),
            "mod duration": np.round(mod_dur,2),
            "ytm": np.round(ytm,2),
            "valor tecnico": np.round(valor_tecnico,2)
        
    } 
    
    return pack["ytm"]
    
    #return datos,ytm
    



######
##### ARRANCAN LOS DATAFRAMES CON YIELDS

año = data_tf["año"].to_list()
mes = data_tf["mes"].to_list()
dia = data_tf["dia"].to_list()
px_to21 = data_tf["ARTO213=BA"].to_list()
px_to23 = data_tf["ARTO233=BA"].to_list()
px_to26 = data_tf["ARTO263=BA"].to_list()

ytm_to21 = []
ytm_to23 = []
ytm_to26 = []

for i in range(0,len(px_to21)):
    ytm_to21.append(to21(px_to21[i],dt.datetime(año[i],mes[i],dia[i])))
for j in range(0,len(px_to23)):
    ytm_to23.append(to23(px_to23[j],dt.datetime(año[j],mes[j],dia[j])))
for h in range(0,len(px_to26)):
    ytm_to26.append(to26(px_to26[h],dt.datetime(año[h],mes[h],dia[h])))
    
to21_b["ytm to21"] = ytm_to21
to23_b["ytm to23"] = ytm_to23
to26_b["ytm to26"] = ytm_to26

ytm_tf = pd.DataFrame(index = data_tf["Date"])
ytm_tf["to21"] = ytm_to21
ytm_tf["to23"] = ytm_to23
ytm_tf["to26"] = ytm_to26

px_tf = pd.DataFrame(index = data_tf["Date"])
px_tf["to21"] = px_to21
px_tf["to23"] = px_to23
px_tf["to26"] = px_to26


def spread_bonos(precio1,precio2,bono1,bono2,index = "2020-01-02"):
    data = pd.DataFrame()
    data["precio 1"] = precio1
    data["precio 2"] = precio2
    data["ytm 1"] = bono1
    data["ytm 2"] = bono2
    data = data.loc[index:]
    
    data["spread"] = data["ytm 1"] - data["ytm 2"]
    data["avg spread"] = np.round(np.mean(data["spread"]),2)
    data["-std"] = np.round(data["avg spread"] - 2 * np.std(data["spread"]),2)
    data["+std"] = np.round(data["avg spread"] + 2 * np.std(data["spread"]),2)
    
    
    return data


año_cer = data_cer["año"].to_list()
mes_cer = data_cer["mes"].to_list()
dia_cer = data_cer["dia"].to_list()
cer_hoy = data_cer["CLOSE"].to_list()

px_tx21 = tx21_b["ARTX213=BA"].to_list()
ytm_tx21 = []

for z in range(0,len(px_tx21)):
    ytm_tx21.append(tx21(px_tx21[z],dt.datetime(año_cer[z],mes_cer[z],dia_cer[z]),cer_hoy[z],cupon = 0.01, cer_emision = 19.2205))

tx21_b["ytm tx21"] = ytm_tx21

px_tx22 = tx22_b["ARTX223=BA"].to_list()
ytm_tx22 = []
for q in range(0,len(px_tx22)):
    ytm_tx22.append(tx22(px_tx22[q],dt.datetime(año_cer[q],mes_cer[q],dia_cer[q]),cer_hoy[q]))
tx22_b["ytm tx22"] = ytm_tx22

px_t2x2 = t2x2_b["ART2X23=BA"].to_list()
ytm_t2x2 = []
for e in range(0,len(px_t2x2)):
    ytm_t2x2.append(t2x2(px_t2x2[e],dt.datetime(año_cer[e],mes_cer[e],dia_cer[e]),cer_hoy[e]))
    
t2x2_b["ytm t2x2"] = ytm_t2x2

px_tx23 = tx23_b["ARTX233=BA"].to_list()

ytm_tx23 = []

for r in range(0,len(px_tx23)):
    ytm_tx23.append(tx23(px_tx23[r],dt.datetime(año_cer[r],mes_cer[r],dia_cer[r]),cer_hoy[r]))
    
tx23_b["ytm tx23"] = ytm_tx23

px_tx24 = tx24_b["ARTX243=BA"].to_list()
ytm_tx24 = []

for t in range(0,len(px_tx24)):
    ytm_tx24.append(tx24(px_tx24[t],dt.datetime(año_cer[t],mes_cer[t],dia_cer[t]),cer_hoy[t]))

tx24_b["ytm tx24"] = ytm_tx24

px_tx26 = tx26_b["ARTX263=BA"].to_list()
ytm_tx26 = []
for y in range(0,len(px_tx26)):
    ytm_tx26.append(tx26(px_tx26[y],dt.datetime(año_cer[y],mes_cer[y],dia_cer[y]),cer_hoy[y]))
    
    
tx26_b["ytm tx26"] = ytm_tx26

px_tc21 = tc21_b["ARTC213=BA"].to_list()
ytm_tc21 = []

for u in range(0,len(px_tc21)):
    ytm_tc21.append(tc21(px_tc21[u],dt.datetime(año_cer[u],mes_cer[u],dia_cer[u]),cer_hoy[u]))
    
tc21_b["ytm tc21"] = ytm_tc21


yields_cer = pd.DataFrame(index = data_cer["fecha cer"])
yields_cer["tx21"] = ytm_tx21
yields_cer["tx22"] = ytm_tx22
yields_cer["t2x2"] = ytm_t2x2
yields_cer["tx23"] = ytm_tx23
yields_cer["tx24"] = ytm_tx24
yields_cer["tx26"] = ytm_tx26
yields_cer["tc21"] = ytm_tc21

px_cer = pd.DataFrame(index = data_cer["fecha cer"])
px_cer["tx21"] = px_tx21
px_cer["tx22"] = px_tx22
px_cer["t2x2"] = px_t2x2
px_cer["tx23"] = px_tx23
px_cer["tx24"] = px_tx24
px_cer["tx26"] = px_tx26
px_cer["tc21"] = px_tc21



def spread_cer(bono1,bono2,precio1,precio2,dia = "2020-01-02"):
    df = pd.DataFrame()
    df["yield 1"] = bono1
    df["yield 2"] = bono2
    df["precio 1"] = precio1
    df["precio 2"] = precio2
    df = df.loc[dia:]
    df["spread"] = (df["yield 1"] - df["yield 2"]) * 100
    df["avg spread"] = np.mean(df["spread"])
    
    df["-std"] = df["avg spread"] - np.std(df["spread"])
    df["+std"] = df["avg spread"] + np.std(df["spread"])
    
    
    return df

def tabla_cer(df,target_bond = "tx24",dia = "2020-05-01"):
    df = df.loc[dia:]
    target = df[target_bond]
    spread = pd.DataFrame(index = df.index)
    columns = list(df)
    for i in columns:
        spread[i] = (target - df[i]) * 100
        
        avg_spread = np.mean(spread)
        men_std = avg_spread - np.std(spread)
        mas_std = avg_spread + np.std(spread)
        data = pd.DataFrame()
        data["actual spread"] = spread.iloc[-1]
        data["avg spread"] = np.round(avg_spread,2)
        data["- std"] = np.round(men_std,2)
        data["+ std"] = np.round(mas_std,2)
        
        
        
    
    return data







st.title("Arbitrajes bonos CER")
ticker1 = st.text_input("ticker CER 1:")
ticker2 = st.text_input("Ticker CER 2:")
fecha = st.text_input("fecha comienzo YYYY/MM/DD")
cer_bonds = spread_cer(yields_cer[ticker1],yields_cer[ticker2],px_cer[ticker1],px_cer[ticker2], dia = fecha)


st.write(cer_bonds)
st.line_chart(cer_bonds[["spread","avg spread","-std","+std"]])




st.title("Tabla comparativa de spreads")
bono_target = st.text_input("ingresar bono para calcular spreads:")
fecha_tabla = st.text_input("fecha para calculo:")
tabla = tabla_cer(yields_cer,bono_target, fecha_tabla)
st.write(tabla)



st.title("Arbitraje bonos TF")
tf1 = st.text_input("ticker bono tasa fija 1:")
tf2 = st.text_input("ticker bono tasa fija 2:")
fecha_tf = st.text_input("fecha comienzo tasa fija YYYY/MM/DD")
tf_bonds = spread_bonos(px_tf[tf1], px_tf[tf2], ytm_tf[tf1], ytm_tf[tf2], index = fecha_tf)
st.write(tf_bonds)
st.line_chart(tf_bonds[["spread","avg spread","-std","+std"]])

