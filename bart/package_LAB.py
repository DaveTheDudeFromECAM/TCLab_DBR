import numpy as np

import matplotlib.pyplot as plt
from IPython.display import display, clear_output

#-----------------------------------        
def LeadLag_RT(MV, Kp, Tlead, Tlag, Ts, PV, PVInit=0, method='EBD'):
    
    """
    The function "Lid-lag_RT" needs to be included in a "for or while loop".
    
    :MV: input vector
    :Kp: process gain
    
    :Tlead: lead time constant [s]
    :Tlag: lag time constant [s]

    :Ts: sampling period [s]
    :PV: output vector
    :PVInit: (optional: default value is 0)
    :method: discretisation method (optional: default value is 'EBD')
        EBD: Euler Backward difference
        EFD: Euler Forward difference
        TRAP: Trapezoïdal method
    
    The function "FO_RT" appends a value to the output vector "PV".
    The appended value is obtained from a recurrent equation that depends on the discretisation method.
    """    
    
    if (Tlead != 0):
        K = Ts/Tlag
        if len(PV) == 0:
            PV.append(PVInit)
        else: # MV[k+1] is MV[-1] and MV[k] is MV[-2]
            if method == 'EBD':
                #PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*MV[-1])
                
                PV.append( (1/(1+K))*PV[-1] + (K*Kp/(1+K)) * ( (1+(Tlead/Ts))*MV[-1]-(Tlead/Ts)* MV[-2] ) )
                
            elif method == 'EFD':
                PV.append((1-K)*PV[-1] + K*Kp*MV[-2])
            elif method == 'TRAP':
                PV.append((1/(2*T+Ts))*((2*T-Ts)*PV[-1] + Kp*Ts*(MV[-1] + MV[-2])))            
            else:
                PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*MV[-1])
    else:
        PV.append(Kp*MV[-1])

#-----------------------------------
def PID_RT(SP, PV, Man, MVMan, MVFF, Kc, Ti, Td, alpha, Ts, MVMin, MVMax, MV, MVP, MVI, MVD, E, ManFF=False, PVInit=0, method="EBD-EBD"):
    
    
    methodI= method.split("-")[0]
    methodD = method.split("-")[1]
    
    # Initialisation of E
    if len(PV) == 0 :
        E.append(SP[-1] - PVInit)
    else:
        E.append(SP[-1] - PV[-1])
    
# slide 194

    #Proportional Action
    MVP.append(Kc*E[-1])
        
    #Integral Action
    if len(MVI) == 0:
        MVI.append((Kc*Ts/Ti)*E[-1])
    else:
        if methodI == 'TRAP':
            MVI.append(MVI[-1] + (0.5*Kc*Ts/Ti)*(E[-1] + E[-2]))
        else:
            MVI.append(MVI[-1] + (Kc*Ts/Ti)*E[-1])
                
                
    #Derivating Action
    Tfd = alpha*Td
    if Td > 0:
        if len(MVD) ==0:
            if len(E) ==1:
                MVD.append((Kc*Td/(Tfd+Ts))*(E[-1]))
            else:
                MVD.append((Kc*Td/(Tfd+Ts))*(E[-1]-E[-2]))
        else:
            if len(E) == 1:
                MVD.append((Tfd/(Tfd+Ts))*MVD[-1]+(Kc*Td/(Tfd+Ts))*(E[-1]))
            else:
                MVD.append((Tfd/(Tfd+Ts))*MVD[-1]+(Kc*Td/(Tfd+Ts))*(E[-1]-E[-2]))
                          
    #value = MVP[-1]+MVI[-1]+MVD[-1]
            
    #Mode manuel et anti wind_up  
    #init MVFF
    if ManFF == True:
        MVFFi = MVFF[-1]
    else:
        MVFFi = 0
        
                   
    if Man[-1]:       
        if ManFF:
            MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1]
        else:
            MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1] - MVFFi
    
    value = MVP[-1]+MVI[-1]+MVD[-1]+MVFFi
                   
    #Min et Max de MV
    if value >= MVMax:
        MVI[-1] = MVMax - MVP[-1] - MVD[-1] - MVFFi
        value = MVMax
                   
    if value <= MVMin:
        MVI[-1] = MVMin - MVP[-1] - MVD[-1] - MVFFi
        value = MVMin
    
    #Update MV
    MV.append(value)

#-----------------------------------
def IMC_tuning(Kp, Tlag1 = 0, Tlag2=0, theta=0, gamma=0, process="FOPDT", model="broida_simple", Tg=0, Tu=0, a=0, t1=0, t2=0):
    
# """
# The function "imc_tuning" is only for first and second order systems.
# :Kp: process gain
# :Tlag1: first (or main) lag time constant [s] used in your process
# :Tlag2: second lag time constant [s] used in your process
# :theta: delay [s] used in your process
# :gamma : constant used to get the closed loop time constant
# :process: process order (ex : FOPDT first order system wuth delay)
# :model: broida_simple or broida_complex for FOPDT
# :Tg:
# :Tu:
# :a:
# :t1: time for 28% of PV (100% being the steady state)
# :t2: time for 44% of PV
# :return: imc tuning parameters respectively:
# - Kc: controller gain
# - Ti: reset time
# - Td: derivative time
# The function "imc_tuning" returns the parameteres that you will use in your PID depending on your process parameters
# """
    if (process == "FOPDT"):
        if (model == "broida_simple"):
            Tlag1 = Tg
            theta = Tu
        elif (model == "broida_complex"):
            Tlag1 = 5.5*(t2 - t1)
            theta = (2.8*t1) - (1.8*t2)

        Tc = gamma * Tlag1
        Kc = ((Tlag1 + theta/2) / (Tc + theta/2)) / Kp
        Ti = Tlag1 + theta/2
        Td = (Tlag1*theta) / (2*Tlag1 + theta)

    elif (process == "SOPDT"):
        if (model == "vdG"):
            Tlag1 = Tg * ((3*a*math.exp(1) - 1) / (1 + a*math.exp(1)))
            Tlag2 = Tg * ((1 - a*math.exp(1)) / (1 + a*math.exp(1)))
            theta = Tu - ((Tlag1*Tlag2) / (Tlag1 + 3*Tlag2))

        Tc = gamma * Tlag1
        Kc = ((Tlag1 + Tlag2) / (Tc + theta)) / Kp
        Ti = Tlag1 + Tlag2
        Td = (Tlag1*Tlag2) / (Tlag1 + Tlag2)

    # else :
    # an = Tu / Tg
    # table = [ [0.0, 1.0], [0.10, 2.72], [0.22, 3.69],[0.32, 4.46],[0.41, 5.12],[0.49, 5.70],[0.57, 6.23] ]

    # for i in range(len(table)):
    # if ( table[i][0] <= an < table[i+1][0]):
    # n = i + 1
    # bn = table[i][1]

    # Tlag1 = Tg / bn
    # Tuth = an * Tg
    # theta = Tu - Tuth
    return Kc, Ti, Td
        