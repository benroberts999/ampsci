#!/usr/bin/python3

def readFile(fname):
    from operator import itemgetter
    ifile = open(fname, "r")
    data = []
    for line in ifile.readlines():
        point = line.split(' ')
        try:
            float(point[0])
        except ValueError:
            continue
        x = float(point[0])
        y = float(point[1])
        data.append([x, y])
    return(data)
#

def write3D(filename, data):
    with open(filename,'w') as f:
        for xy in data:
            o=''
            for val in xy:
                o = o + str(val) + " "
            print(o, file=f)
#

def writeFavourdedRegion(filename, data0, data_min, data_max, nsig):
    import math as m

    # DAMA observed in the 1--2 keV region
    Obs = 0.0198
    dObs= 0.0067  #1.65, for 90% C.L.
    se0 = 1.e-37

    out = []
    oa_min = 1
    for i in range(len(data0)): #blah
        Mx = data0[i][0]
        S0 = data0[i][1]
        dS_max = abs(data_max[i][1] - S0)
        dS_min = abs(S0 - data_min[i][1])
        #if data_min[i][1] == 0:
        if S0 == 0:
            #out.append([Mx, 0, 0, 0])
            continue
        #if data_min[i][1] > 0 and data_min[i][1] < oa_min:
        #    oa_min = data_min[i][1]
        if data_min[i][1] == 0:
            data_min[i][1] = 1.e-7
            # just a v. small number.. [destructive!!] ??
        se = se0*Obs/S0
        arg_min = pow(dObs/Obs,2) + pow(dS_min/S0,2)
        arg_max = pow(dObs/Obs,2) + pow(dS_max/S0,2)
        er_min = m.sqrt(arg_min)*se
        er_max = m.sqrt(arg_max)*se
        se_min = se0*(Obs-dObs)/(data_max[i][1])
        se_max = se0*(Obs+dObs)/data_min[i][1]
        #out.append([Mx, se, se-er_min, se+er_max])
        out.append([Mx, se, se_min, se_max])
        #out.append([Mx, se, se_min, se+er_max])
    write3D(filename,out)
#



data0=readFile("Sm_mx-ak-I_S1-Sm0_h.out")
data_min=readFile("Sm_mx-ak-I_S1-SmMin_h.out")
data_max=readFile("Sm_mx-ak-I_S1-SmMax_h.out")
nsig = 1.6
writeFavourdedRegion("heavy-favoured.txt",data0,data_min,data_max,nsig)
