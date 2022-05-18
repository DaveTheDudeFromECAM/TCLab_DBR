#FINAL VERSION

import math
import numpy as np
from package_DBR import *
import matplotlib.pyplot as plt
from IPython.display import display, clear_output

#-----------------------------------        
def Lead_Lag_RT(MV, Kp, Tlead, Tlag, Ts, PV, PVInit=0, method='EBD'):
    
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
class Controller:
    
    def __init__(self, parameters):
        
        self.parameters = parameters
        self.parameters['Kc'] = parameters['Kc'] if 'Kc' in parameters else 2.0
        self.parameters['Ti'] = parameters['Ti'] if 'Ti' in parameters else 100.0
        self.parameters['Td'] = parameters['Td'] if 'Td' in parameters else 10.0
        self.parameters['alpha'] = parameters['alpha'] if 'aplha' in parameters else 1.0

#----------------------------------- 
def PID_RT(SP, PV, Man, MVMan, MVFF, Kc, Ti, Td, alpha, Ts, MVMin, MVMax, MV, MVP, MVI, MVD, E, ManFF=False, PVInit=0, method="EBD-EBD"):

    """
    Help on function PID_RT in module package_LAB:

    PID_RT (SP, PV, Man, MVMan, MVFF, Kc, Ti, Td, alpha, Ts, MVMin, MVMax, MV, MVP, MVI, MVD, E, ManFF=False, PVInit=@, method='EBD-EBD')
    The function "PID RT"
    needs to be included in a
    "for or while loop"
    
    SP: SP (or SetPoint) vector
    PV: PV (or Process Value) vector
    Man: Man (or Manual controller mode) vector (True or False)
    MVMan: MVMan (or Manual value for MV) vector
    MVFF: MVFF (or Feedforward) vector

    Kc: controller gain
    Ti: integral time constant [s]
    Td: derivative time constant [s]
    alpha: Tfd = alpha*Td where Tfd is the derivative filter time constant [s]
    Ts: sampling period [s]

    MVMin: minimum value for MV (used for saturation and anti wind-up)
    MVMax: maximum value for MV (used for saturation and anti wind-up)

    MV: MV (or Manipulated Value) vector

    MVP: MVP (or Propotional part of MV) vector
    MVI: MVI (or Integral part of MV) vector
    MVD: MVD (or Derivative part of MV) vector

    E: E (or control Error) vector
    ManFF: Activated FF in manual mode (optional: default boolean value is False)
    PVInit: Initial value for PV (optional: default value is e): used if PID_RT is ran first in the squence and no value of PV is available yet.
    
    method: discretisation method (optional: default value is 'EBD')
    EBD-EBD: EBD for integral action and EBD for derivative action
    EBD-TRAP: EBD for integral action and TRAP for derivative action
    TRAP-EBD: TRAP for integral action and EBD for derivative action
    TRAP-TRAP: TRAP for integral action and TRAP for derivative action
    
    The function "PID_ RT"
    appends new values to the vectors "MV",
    "MVP". "MVI".
    and "MVD"
    The appended values are based on the PID algorithm, the controller mode, and feedforward.
    Note that saturation of "MV" within the limits [MVMin MVMax] is implemented with anti wind-up.
    """
    
    

    
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
def IMC_tuning(Kp, Tlag1 = 0, Tlag2=0, theta=0, gamma=0, process="FOPDT"):  

    """
    The function "imc_tuning" is only for first and second order systems.
    :Kp: process gain
    :Tlag1: first (or main) lag time constant [s] used in your process
    :Tlag2: second lag time constant [s] used in your process
    :theta: delay [s] used in your process
    :gamma : constant used to get the closed loop time constant
    :process: process order (ex : FOPDT first order system wuth delay)
    :model: broida_simple or broida_complex for FOPDT
    :Tg:
    :Tu:
    :a:
    :t1: time for 28% of PV (100% being the steady state)
    :t2: time for 44% of PV
    :return: imc tuning parameters respectively:
    - Kc: controller gain
    - Ti: reset time
    - Td: derivative time
    The function "imc_tuning" returns the parameteres that you will use in your PID depending on your process parameters
    """
    
    if (process == "FOPDT"):

        Tc = gamma * Tlag1
        Kc = ((Tlag1 + theta/2) / (Tc + theta/2)) / Kp
        Ti = Tlag1 + theta/2
        Td = (Tlag1*theta) / (2*Tlag1 + theta)

    elif (process == "SOPDT"):

        Tc = gamma * Tlag1
        Kc = ((Tlag1 + Tlag2) / (Tc + theta)) / Kp
        Ti = Tlag1 + Tlag2
        Td = (Tlag1*Tlag2) / (Tlag1 + Tlag2)

    return Kc, Ti, Td

#----------------------------------- 

def StabilityMargins(P: Process, C: Controller, omega):
    """
    The function "stability_margins" needs to have 2 processes object in paramaters.
    
    :P: the system process
    :C: the controller 
    :omega: frequency vector (rad/s); generated by a command of the type "omega = np.logspace(-2, 2, 10000)".
        
    The function "stability_margins" generates the bodes plots of the Loop gain and gives the phase and gain margins.
    """    
    
    L = Process({})
    L.parameters['Kp'] = P.parameters['Kp']  + C.parameters['Kc'] 
    L.parameters['Tlag1'] = P.parameters['Tlag1']  + C.parameters['Ti'] 
    L.parameters['Tlag2'] = P.parameters['Tlag2']  + C.parameters['Td'] 
    L.parameters['theta'] = P.parameters['theta']  + C.parameters['alpha'] 
    
    Ls = Bode(L,omega,Show=False)
    
    gain_values = 20*np.log10(np.abs(Ls))    
    phase_values = (180/np.pi)*np.unwrap(np.angle(Ls))
    
    for i in range(len(gain_values)):
        if gain_values[i] > 0 and gain_values[i + 1] < 0: 
            x_gain = i
            y_gain = gain_values[i]
            Wc = round(omega[i],2)
            phase_margin = round(abs(phase_values[i] + 180),2)

    
    for i in range(len(phase_values)):
        if phase_values[i] > -180 and phase_values[i + 1] < -180:     
            x_phase = i
            y_phase = phase_values[i]
            Wu = round(omega[i],2)
            gain_margin = round(abs(gain_values[i]), 2)
    
        
    fig, (ax_gain, ax_phase) = plt.subplots(2,1)
    fig.set_figheight(12)
    fig.set_figwidth(22)

    # Gain part
    ax_gain.axhline(y =0, color = 'b', linestyle = '-')
    ax_gain.semilogx(omega,20*np.log10(np.abs(Ls)),label='L(s)')
    gain_min = np.min(20*np.log10(np.abs(Ls)/5))
    gain_max = np.max(20*np.log10(np.abs(Ls)*5))
    ax_gain.vlines(omega[x_phase],gain_min,gain_values[x_phase], color = 'r', linestyle = '--')
    ax_gain.vlines(omega[x_gain],gain_min,0, color = 'blue', linestyle = '--')
    
    #ax_gain.vlines(omega[x_phase],gain_values[x_phase],0, color = 'g')
    
    ax_gain.annotate(
    '', xy=(omega[x_phase], gain_values[x_phase]), xycoords='data',
    xytext=(omega[x_phase],0), textcoords='data',
    arrowprops={'arrowstyle': '<->'})
    
    ax_gain.annotate(
    f"Gain margin = {gain_margin}dB at {Wu}rad/s ", xy=(omega[x_phase], gain_values[x_phase] /2 ), xycoords='data',
    xytext=(5, 0), textcoords='offset points')
    
    ax_gain.set_xlim([np.min(omega), np.max(omega)])
    ax_gain.set_ylim([gain_min, gain_max])
    ax_gain.set_ylabel('Amplitude |L| [db]')
    ax_gain.set_title('Bode plot of L')
    ax_gain.legend(loc='best')
    
    # Phase part
    ax_phase.axhline(y =-180, color = 'b', linestyle = '-')
    ax_gain.vlines(y_gain,-180,y_phase)
    ax_phase.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(Ls)),label='L(s)')
    ax_phase.set_xlim([np.min(omega), np.max(omega)])
    ph_min = np.min((180/np.pi)*np.unwrap(np.angle(Ls))) - 10
    ph_max = np.max((180/np.pi)*np.unwrap(np.angle(Ls))) + 10
    ax_phase.vlines(omega[x_gain],phase_values[x_gain] ,ph_max, color = 'blue', linestyle = '--')
    ax_phase.vlines(omega[x_phase],-180,ph_max,color = 'r', linestyle = '--')
    #ax_phase.vlines(omega[x_gain],-180 ,phase_values[x_gain], color = 'g')
    
    ax_phase.annotate(
    '', xy=(omega[x_gain], -180), xycoords='data',
    xytext=(omega[x_gain], phase_values[x_gain]), textcoords='data',
    arrowprops={'arrowstyle': '<->'})
    
    ax_phase.annotate(
    f"Phase margin = {phase_margin}° at {Wc}rad/s ", xy=(omega[x_gain], (-180 + phase_values[x_gain]) /2 ), xycoords='data',
    xytext=(5, 0), textcoords='offset points')
    
    ax_phase.set_ylim([np.max([ph_min, -200]), ph_max])
    ax_phase.set_ylabel(r'Phase $\angle L$ [°]')
    ax_phase.legend(loc='best')