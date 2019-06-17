#!/usr/bin/env python3

import glob
import os
import ntpath
from datetime import datetime, timedelta
import time
import dateutil.parser

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
plt.switch_backend('agg')
from astropy.io import fits
from astropy.time import Time
import astropy.coordinates as coord
from astropy.time import TimeDelta
from astropy.coordinates import AltAz
import astropy.units as u
import astropy.time
import pandas as pd
import healpy as hp
import pygal  
from pygal.style import LightSolarizedStyle
import numpy as np
from slackclient import SlackClient


slack_client = SlackClient('BOT TOKEN')



def SUN_POS():
    """ 
    Finds sun position in RA and DEC at sunrise and sunset to find
    which parts of the sky are viewable
    """
    LIST=[]
    loc = coord.EarthLocation(lon=28.7134 * u.deg,
                              lat=17.9058 * u.deg)

    start = Time(datetime.utcnow(), format='datetime', scale='utc') - TimeDelta(1, format='jd')

    for i in range(0,384):
        T = start+TimeDelta(i*225, format='sec')
        LST = T.sidereal_time('apparent', '28.7134d')
        TMP1 = str(LST).split('h')
        TMP2 = TMP1[1].split('m')
        TMP3 = TMP2[1].split('s')
        LST = (15*float(TMP1[0]) + 0.0166667*float(TMP2[0]) + 0.000277778*float(TMP3[0]))
        altazframe = AltAz(obstime=T, location=loc)
        sunaltaz = coord.get_sun(T).transform_to(altazframe)
        if sunaltaz.alt.value <= -47 and len(LIST)==0:
            LIST.append(str(coord.get_sun(T).ra).split('d')[0])
            print('Sunset')
            print(T)

        elif sunaltaz.alt.value >= 48 and len(LIST)==1:
            LIST.append(str(coord.get_sun(T).ra).split('d')[0])
            print('Sunrise')
            print(T)
            break

    UP = 30 + (float(LIST[0])+float(LIST[1]))/2
    LOW = (float(LIST[0])+float(LIST[1]))/2 - 30

    return([UP,LOW])



def get_channels(print_ = False, ch = False):
    """
    Returns list of available channels.
    Will also make sure input token is correct.
    Will return specific channel id associated with channel name
    """
    i = 0

    channels_call = slack_client.api_call("channels.list")
    if channels_call.get('ok'):
        channels = channels_call['channels']
        if print_ == True:
            print("Channels: ")
            for c in channels:
                print(c['name'] + " (" + c['id'] + ")")
        if ch == False:
            return(channels)
        else:
            for c in channels:
                if c['name'] == ch:
                    OUT = i
                else:
                    i+=1
            return(channels[OUT]['id'])

    else:
        print("Unable to authenticate.")


def sub_text(channel_name, message):
    """
    Send a message to requested channel_id
    """
    channel_id = get_channels(ch = channel_name)
    RP = slack_client.api_call(
        "chat.postMessage",
        channel=channel_id,
        text=message,
        username='Nightly_report_generator',
        icon_emoji=':robot_face:'
    )

    return(RP)


def sub_file(channel_name, title1, filename1, message = False):
    """
    Upload a file to requested channel_id
    """

    channel_id = get_channels(ch = channel_name)

    if message == False:
        message = 'No message'

    RP = slack_client.api_call('files.upload',
                               filename=filename1,
                               channels=channel_id,
                               file=open(filename1,'rb'),
                               Title= title1,
                               initial_comment = message
    )

    return(RP)


def night_plot(Event=False, SM =False):
    Revisit = []
    DATE = datetime.strftime(datetime.now() - timedelta(1), '%Y-%m-%d')
    F = glob.glob("LOCATION")

    OF = open(DATE+'_report.csv', 'w') #Output_file
    OF.write('File,NCOADD,CALLIM5,QUALITY,TILENAME,EVENT,RA,DEC\n')

    for fil in F:
        H = fits.open(fil)[1]._header
        FName = ntpath.basename(fil)
        OF.write(FName+','+str(H['NCOADD'])+','+str(H['CALLIM5'])+','+str(H['QUALITY'])+','+str(H['TILENAME'])+','+str(H['EVENT']+','+str(H['RA-TARG'])+','+str(H['DEC-TARG'])+'\n'))
    OF.close()


    DF = pd.read_csv(DATE+'_report.csv')
    num_images = DF.shape[0] #number of images taken that night

    u_tn  = DF.TILENAME.unique() #unique TILENAMES over the evening
    os.remove(DATE+'_report.csv')
    TILES_NUM = len(u_tn)
    if Event == False: #plot just area tiled
        NSIDE= 128 #2048
        OVERLAP = []
        dpp = hp.nside2pixarea(NSIDE, degrees=True)
        m = np.zeros(hp.nside2npix(NSIDE))
        B2, B1 = SUN_POS()
        for i in range(len(m)):
            if hp.pix2ang(NSIDE, i, lonlat=True)[1] > 85 or hp.pix2ang(NSIDE, i, lonlat=True)[1] < -35 :
                 m[i] = 1
            elif B1 < hp.pix2ang(NSIDE, i, lonlat=True)[0] < B2:
                 m[i] = 1

        m[10:100] = 5
        hp.mollview(m, cmap='gist_yarg', coord='C', notext = True, title=DATE, flip = 'geo', cbar=False)
        hp.visufunc.graticule(dpar = 15, dmer = 30)
        PD = np.linspace(-180,180, 7)
        PR = np.linspace(0, 120, 9)
        for i in range (len(PD)-1):
            hp.visufunc.projtext(-180, PR[i], '%s  ' %int(PR[i]), lonlat=True,  fontsize = 10,  horizontalalignment = 'right')
            hp.visufunc.projtext(PD[i]-1, -8.5, '   %s' %int(PD[i]), lonlat=True, fontsize = 10, horizontalalignment = 'center')


        for tile in u_tn:
            try:
                Corn = []
                T = DF.loc[DF['TILENAME'] == tile]
                RA = T.iloc[-1]['RA']
                DEC = T.iloc[-1]['DEC']
                n_R = float(RA.split(':')[0])*15 + float(RA.split(':')[1])*0.25 + float(RA.split(':')[2])*0.00417
                n_D = float(DEC.split(':')[0]) + float(DEC.split(':')[1])*0.0167 + float(DEC.split(':')[2])*0.00028

                Corn.append([n_R+1.85, n_D+2.45])
                Corn.append([n_R-1.85, n_D+2.45])
                Corn.append([n_R-1.85, n_D-2.45])
                Corn.append([n_R+1.85, n_D-2.45])
                Corn.append([n_R+1.85, n_D+2.45]) #corners

                for Xs in range(4):
                    x_lin = np.linspace(Corn[Xs][0], Corn[Xs+1][0]+0.0000000001, 40)
                    y_lin = np.linspace(Corn[Xs][1], Corn[Xs+1][1]+0.0000000001, 40)
                    if T.iloc[-1]['NCOADD'] == 1:
                        hp.visufunc.projplot([x_lin,y_lin], color = 'r', lonlat=True, label='stacks<1')
                    elif T.iloc[-1]['NCOADD'] == 2:
                        hp.visufunc.projplot([x_lin,y_lin], color = 'gold', lonlat=True, label='stacks =2')
                    else:
                        hp.visufunc.projplot([x_lin,y_lin], color = 'green', lonlat=True, label='stacks>=3')

                if T.iloc[-1]['QUALITY'] < 6  and tile not in Revisit:
                    Revisit.append(tile)



                VERT = [hp.pix2vec(NSIDE, hp.ang2pix(NSIDE, n_R+3.16, n_D+2.09, lonlat=True)),
                        hp.pix2vec(NSIDE, hp.ang2pix(NSIDE, n_R+3.16, n_D-2.09, lonlat=True)),
                        hp.pix2vec(NSIDE, hp.ang2pix(NSIDE, n_R-3.16, n_D-2.09, lonlat=True)),
                        hp.pix2vec(NSIDE, hp.ang2pix(NSIDE, n_R-3.16, n_D+2.09, lonlat=True))]


                TMP_List = hp.query_polygon(NSIDE,VERT)
                OVERLAP.extend(TMP_List)

            except:
                print('error in RA/DEC')


        Area = len(list(set(OVERLAP)))*dpp
        hp.visufunc.projtext(0 ,-90 , '                             Area = %s$°^2$ ' %round(Area,2)  ,color = 'black', fontsize = 15, lonlat=True)
        red = mpatches.Patch(color='red', label='1 stack')
        gold = mpatches.Patch(color='gold', label='2 stacks')
        green = mpatches.Patch(color='green', label='3 or more stack')
        plt.legend(handles=[red, gold, green], loc= 'lower center')
        plt.savefig(DATE+'.png')
        return(DATE+'.png',round(Area,2), TILES_NUM, Revisit)


    else: #plot sky map, area tiled, percentage region covered
        #Event = skymap#
        prob = 0
        OVERLAP = []
        C_nr = []
        C_nd = []
        skymap,  header = hp.read_map(SM, h= True)
        header = dict(header)
        NSIDE = header['NSIDE']
        max1 = len(skymap)
        dpp = hp.nside2pixarea(NSIDE, degrees=True) #degrees-sqaured per pixel
        st = header['DATE-OBS']

        hp.mollview(skymap, title = DATE, coord = 'C', cmap = 'ocean_r', notext = True, xsize = 1000, flip = 'geo')
        hp.visufunc.graticule(dpar = 15, dmer = 30)
        PD = np.linspace(-180,180, 7)
        PR = np.linspace(0, 120, 9)
        for i in range (len(PD)-1):
            hp.visufunc.projtext(-180, PR[i], '%s  ' %int(PR[i]) , lonlat= True, fontsize = 10,  horizontalalignment = 'right')
            hp.visufunc.projtext(PD[i]-1, -8.5, '   %s' %int(PD[i]) , lonlat= True, fontsize = 10, horizontalalignment = 'center')

        TILES_NUM = 0
        for tile in u_tn:
            T = DF.loc[DF['TILENAME'] == tile]
            if T.iloc[-1]['EVENT'] == Event:
                try:
                    TILES_NUM += 1
                    Corn = []
                    RA = T.iloc[-1]['RA']
                    DEC = T.iloc[-1]['DEC']
                    n_R = float(RA.split(':')[0])*15 + float(RA.split(':')[1])*0.25 + float(RA.split(':')[2])*0.00417
                    n_D = float(DEC.split(':')[0]) + float(DEC.split(':')[1])*0.0167 + float(DEC.split(':')[2])*0.00028

                    C_nr.append([n_R+1.85, n_R-1.85])
                    C_nd.append([n_D+2.45, n_D-2.45])

                    Corn.append([n_R+1.85, n_D+2.45])
                    Corn.append([n_R-1.85, n_D+2.45])
                    Corn.append([n_R-1.85, n_D-2.45])
                    Corn.append([n_R+1.85, n_D-2.45])
                    Corn.append([n_R+1.85, n_D+2.45]) #corners

                    for Xs in range(4):
                        x_lin = np.linspace(Corn[Xs][0], Corn[Xs+1][0]+0.0000000001, 40)
                        y_lin = np.linspace(Corn[Xs][1], Corn[Xs+1][1]+0.0000000001, 40)
                        hp.visufunc.projplot([x_lin,y_lin], lonlat = True, color = 'orange')

                    VERT = [hp.pix2vec(NSIDE, hp.ang2pix(NSIDE, n_R+3.16, n_D+2.09, lonlat=True)),
                            hp.pix2vec(NSIDE, hp.ang2pix(NSIDE, n_R+3.16, n_D-2.09, lonlat=True)),
                            hp.pix2vec(NSIDE, hp.ang2pix(NSIDE, n_R-3.16, n_D-2.09, lonlat=True)),
                            hp.pix2vec(NSIDE, hp.ang2pix(NSIDE, n_R-3.16, n_D+2.09, lonlat=True))]



                    TMP_List = hp.query_polygon(NSIDE,VERT)
                    prob += sum(skymap[TMP_List])
                    skymap[TMP_List] = 0
                    OVERLAP.extend(TMP_List)

                except:
                    print('error in RA/DEC')



        ###Ugly layout here for better plotting later
        Area = len(list(set(OVERLAP)))*dpp
        hp.visufunc.projtext(-15, -70 , '                                          Prob = %s %%' %round(prob*100,2)  ,color = 'black', fontsize = 15, lonlat=True)
        hp.visufunc.projtext(0 ,-90 , '                             Area = %s$°^2$ ' %round(Area,2)  ,color = 'black', fontsize = 15, lonlat=True)
        plt.savefig(Event+'_'+DATE+'.png')
        hp.projaxes.GnomonicAxes.legend( ('3 or more stacks', '2 stacks', '1 stack'), loc='upper left')
        return(Event+'_'+DATE+'.png', round(Area,2), round(prob*100,2), TILES_NUM)


#night_plot(Event = 'LVC_S190426c', SM = 'bayestar_no_virgo.fits')

def event_night_report(EVENT, skymap, GCN1, GCN2, DT1, DT2):
    DATE = datetime.strftime(datetime.now() - timedelta(1), '%Y-%m-%d')
    OUT = open(EVENT+'_'+DATE+'_GCN.txt', 'w')
    Template = open('template.txt','r')
    PNG, AREA, PROB, TILES = night_plot(Event = EVENT, SM = skymap)



    for lin in Template:
        if 'AREA' in lin:
            lin = lin.replace('AREA', str(AREA))
        if 'EVENT_NAME' in lin:
            lin = lin.replace('EVENT_NAME', EVENT)
        if 'PERCENTAGE' in lin:
            lin = lin.replace('PERCENTAGE', str(PROB))
        if 'DT1' in lin:
            lin =lin.replace('DT1', DT1)
        if 'DT2' in lin:
            lin = lin.replace('DT2', DT2)
        if 'GCN_NUM2' in lin:
            lin = lin.replace('GCN_NUM2', GCN2)
        if 'GCN_NUM' in lin:
            lin = lin.replace('GCN_NUM', GCN1)

        OUT.write(lin)
    OUT.close()
    M = 'Last night GOTO covered '+str(AREA)+' deg sq. \n GOTO has a '+str(PROB)+'%% of seeing the counterpart \n'+'in '+str(TILES)+'pointings'
    RP = sub_file('test1234', PNG.replace('.png','') , PNG, message=M)
    return(RP)

def night_report_gen():
    PNG, AREA, TILES, Revisit = night_plot()
    if len(Revisit) > 0:
        H = 'Suggested tiles to Revisit:\n'
        for Re in Revisit:
            H += Re+'\n'
        M = 'Last night GOTO covered '+str(AREA)+' deg squared in '+ str(TILES)+' tiles \n'+H+' \nThis message was sent from the nightly report generator'

    else:
        M = 'Last night GOTO covered '+str(AREA)+' deg squared in '+ str(TILES)+' tiles \nThis message was sent from the nightly report generator'

    RP = sub_file('channel', PNG.replace('.png','') , PNG, message=M)
    return(RP)

#print(event_night_report('LVC_S190426c', 'LALInference1.fits', '#2345', '#12343', 'fol', 'lol'))
print(night_report_gen())




