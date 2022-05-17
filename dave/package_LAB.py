import numpy as np

import matplotlib.pyplot as plt
from IPython.display import display, clear_output

#-----------------------------------
def myRound(x, base=5):
    
    """
    Returns a float that is the closest multiple of "base" near "x"
    Based on: https://stackoverflow.com/questions/2272149/round-to-5-or-other-number-in-python
    
    :x: parameter that is rounded to a multiple of "base"
    :base: "base" parameter (optional: default value is 5)
    
    :return: rounded parameter
    """
    
    return float(base * round(float(x)/base))

#-----------------------------------
def SelectPath_RT(path,time,signal):
    
    """
    The function "SelectPath_RT" needs to be included in a "for or while loop".
    
    :path: dictionary input describing a path in time. Example: path = {0: 0, 5: 1, 50: 2, 80: 3, 100: 3}
    :time: time vector.
    :signal: signal vector that is being constructed using the input "path" and the vector "time".
    
    The function "SelectPath_RT" takes the last element in the vector "time" and, given the input "path", it appends the correct value to the vector "signal".
    """    
    
    for timeKey in path:
        if(time[-1] >= timeKey):
            timeKeyPrevious = timeKey    
    
    value = path[timeKeyPrevious]
    signal.append(value)

#-----------------------------------
def Delay_RT(MV,theta,Ts,MV_Delay,MVInit=0):
    
    """
    The function "Delay_RT" needs to be included in a "for or while loop".
    
    :MV: input vector
    :theta: delay [s]
    :Ts: sampling period [s]
    :MV_Delay: delayed input vector
    :MVInit: (optional: default value is 0)
    
    The function "Delay_RT" appends a value to the vector "MV_Delay".
    The appended value corresponds to the value in the vector "MV" "theta" seconds ago.
    If "theta" is not a multiple of "Ts", "theta" is replaced by Ts*int(np.ceil(theta/Ts)), i.e. the closest multiple of "Ts" larger than "theta".
    If the value of the vector "input" "theta" seconds ago is not defined, the value "MVInit" is used.
    """
    
    NDelay = int(np.ceil(theta/Ts))
    if NDelay > len(MV)-1:
        MV_Delay.append(MVInit)
    else:    
        MV_Delay.append(MV[-NDelay-1])

#-----------------------------------        
def FO_RT(MV,Kp,T,Ts,PV,PVInit=0,method='EBD'):
    
    """
    The function "FO_RT" needs to be included in a "for or while loop".
    
    :MV: input vector
    :Kp: process gain
    :T: lag time constant [s]
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
    
    if (T != 0):
        K = Ts/T
        if len(PV) == 0:
            PV.append(PVInit)
        else: # MV[k+1] is MV[-1] and MV[k] is MV[-2]
            if method == 'EBD':
                PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*MV[-1])
            elif method == 'EFD':
                PV.append((1-K)*PV[-1] + K*Kp*MV[-2])
            elif method == 'TRAP':
                PV.append((1/(2*T+Ts))*((2*T-Ts)*PV[-1] + Kp*Ts*(MV[-1] + MV[-2])))            
            else:
                PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*MV[-1])
    else:
        PV.append(Kp*MV[-1])

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
def FOPDT(MV,Kp,T,theta,Ts,MVInit=0,PVInit=0,method='EBD'):
    
    """
    The function "FOPDT" DOES NOT need to be included in a "for or while loop": this block is for offline use.
    
    :MV: input vector
    :Kp: process gain
    :T: lag time constant [s]
    :theta: delay [s]
    :Ts: sampling period [s]
    :MVInit: (optional: default value is 0)    
    :PVInit: (optional: default value is 0)
    :method: discretisation method (optional: default value is 'EBD')
        EBD: Euler Backward difference
        EFD: Euler Forward difference
        TRAP: Trapezoïdal method
        
    :return: simulated FOPDT output vector         
    
    The function "FOPDT" returns the simulated output FOPDT vector from the input vector "MV" and the input parameters.
    """    
    
    MVDelay = []
    MVTemp = []
    PVSim = []    
    
    for i in range(0,len(MV)):
        MVTemp.append(MV[i])
        Delay_RT(MVTemp,theta,Ts,MVDelay,MVInit)
        FO_RT(MVDelay,Kp,T,Ts,PVSim,PVInit,method)
            
    return PVSim

#-----------------------------------
def SOPDT(MV,Kp,T1,T2,theta,Ts,MVInit=0,PVInit=0,method='EBD'):
    
    """
    The function "SOPDT" DOES NOT need to be included in a "for or while loop": this block is for offline use.
    
    :MV: input vector
    :Kp: process gain
    :T1: first (or main) lag time constant [s]
    :T2: second lag time constant [s]    
    :theta: delay [s]
    :Ts: sampling period [s]
    :MVInit: (optional: default value is 0)    
    :PVInit: (optional: default value is 0)
    :method: discretisation method (optional: default value is 'EBD')
        EBD: Euler Backward difference
        EFD: Euler Forward difference
        TRAP: Trapezoïdal method
        
    :return: simulated SOPDT output vector         
    
    The function "SOPDT" returns the simulated SOPDT output vector from the input vector "MV" and the input parameters.
    """     
    
    MVDelay = []
    MVTemp = []
    PV1 = []
    PVSim = []    
    
    for i in range(0,len(MV)):
        MVTemp.append(MV[i])
        Delay_RT(MVTemp,theta,Ts,MVDelay,MVInit)
        FO_RT(MVDelay,Kp,T1,Ts,PV1,PVInit,method)
        FO_RT(PV1,1,T2,Ts,PVSim,PVInit,method)
            
    return PVSim

#-----------------------------------
def FOPDT_cost(p,MV,PV,Ts,*args):
    
    """
    :p: parameter vector:
        Kp = p[0]: process gain
        T = p[1]: lag time constant [s]
        theta = p[2]: delay [s]
    :MV: input vector used during experimentation
    :PV: experimental output vector obtained in response to the input vector MV
    :Ts: sampling period [s]
    :args: object, axes and line handles for representing PV and the simulated PV at each function call (optional)
        fig: figure object
        ax1: axes object
        l1: line object for PV
        l2: line object for simulated PV
        
    :return: identification cost with FOPDT model
    
    The function "FOPDT_cost" returns the identification cost, i.e. the sum of the model errors squared.
    The model error is the difference between the experimental output "PV" and the simulated output with a FOPDT model and the parameter set "p".
    
    The assumption is that MVInit and PVInit are zero for the use of Delay_RT and FO_RT in the code.
    The 'EBD' discretisation method is used.
    This assumption is met after a "cleaning operation" after which MV[0] = 0 and PV[0] are 0.
    """     
    
    Kp = p[0]
    T = p[1]
    theta = np.max((0,p[2]))
    
    MVDelay = []
    MVTemp = []
    PVSim = []
    t = []
    
    objective = 0
    
    for i in range(0,len(MV)):
        t.append(i*Ts)
        MVTemp.append(MV[i])
        Delay_RT(MVTemp,theta,Ts,MVDelay)
        FO_RT(MVDelay,Kp,T,Ts,PVSim)
        objective = objective + (PV[i] - PVSim[i])**2      
    
    for fig, ax1, l1, l2 in args:
        l1.set_data(t,PV)
        l2.set_data(t,PVSim)
        ax1.set_xlim(0, t[-1]+1)
        ax1.set_ylim(myRound(np.min(PV),5)-1, myRound(np.max(PV),5)+1)        
        # Comment following line otherwise optimisation too slow
        # ax1.text(100, 0.1, 'Kp = {0:3.2f}, T = {1:3.2f}, theta = {2:3.2f}'.format(Kp, T, theta), bbox={'facecolor': 'white'})
        clear_output(wait=True)
        display(fig)
            
    return objective

#-----------------------------------
def SOPDT_cost(p,MV,PV,Ts,*args):
    
    """
    :p: parameter vector:
        Kp = p[0]: process gain
        T1 = p[1]: first or main lag time constant [s]
        T2 = p[2]: second lag time constant [s]    
        theta = p[3]: delay [s]
    :MV: input vector used during experimentation
    :PV: experimental output vector obtained in response to the input vector MV
    :Ts: sampling period [s]
    :args: object, axes and line handles for representing PV and the simulated PV at each function call (optional)
        fig: figure object
        ax1: axes object
        l1: line object for PV
        l2: line object for simulated PV
        
    :return: identification cost with SOPDT model
    
    The function "SOPDT_cost" returns the identification cost, i.e. the sum of the model errors squared.
    The model error is the difference between the experimental output "PV" and the simulated output with a SOPDT model and the parameter set "p".
    
    The assumption is that MVInit and PVInit are zero for the use of Delay_RT and FO_RT in the code.
    The 'EBD' discretisation method is used.
    This assumption is met after a "cleaning operation" after which MV[0] = 0 and PV[0] are 0.
    """    
    
    Kp = p[0]
    T1 = p[1]
    T2 = p[2]
    theta = np.max((0,p[3]))    
    
    MVDelay = []
    MVTemp = []
    PV1 = []
    PVSim = []
    t = []
    
    objective = 0
    
    for i in range(0,len(MV)):
        t.append(i*Ts)
        MVTemp.append(MV[i])
        Delay_RT(MVTemp,theta,Ts,MVDelay)
        FO_RT(MVDelay,Kp,T1,Ts,PV1)
        FO_RT(PV1,1,T2,Ts,PVSim)
        objective = objective + (PV[i] - PVSim[i])**2      
    
    for fig, ax1, l1, l2 in args:
        l1.set_data(t,PV)
        l2.set_data(t,PVSim)
        ax1.set_xlim(0, t[-1]+1)
        ax1.set_ylim(myRound(np.min(PV),5)-1, myRound(np.max(PV),5)+1)        
        # Comment following line otherwise optimisation too slow
        # ax1.text(100, 0.1, 'Kp = {0:3.2f}, T1 = {1:3.2f}, T2 = {2:3.2f}, theta = {3:3.2f}'.format(Kp, T1, T2, theta), bbox={'facecolor': 'white'})
        clear_output(wait=True)
        display(fig)
            
    return objective

#-----------------------------------        
class Process:
    
    def __init__(self, parameters):
        
        self.parameters = parameters
        self.parameters['Kp'] = parameters['Kp'] if 'Kp' in parameters else 1.0
        self.parameters['theta'] = parameters['theta'] if 'theta' in parameters else 0.0
        self.parameters['Tlead1'] = parameters['Tlead1'] if 'Tlead1' in parameters else 0.0
        self.parameters['Tlead2'] = parameters['Tlead2'] if 'Tlead2' in parameters else 0.0
        self.parameters['Tlag1'] = parameters['Tlag1'] if 'Tlag1' in parameters else 0.0
        self.parameters['Tlag2'] = parameters['Tlag2'] if 'Tlag2' in parameters else 0.0

#-----------------------------------         
class Controller:
    
    def __init__(self, parameters):
        
        self.parameters = parameters
        self.parameters['Kc'] = parameters['Kc'] if 'Kc' in parameters else 2.0
        self.parameters['Ti'] = parameters['Ti'] if 'Ti' in parameters else 100.0
        self.parameters['Td'] = parameters['Td'] if 'Td' in parameters else 10.0
        self.parameters['alpha'] = parameters['alpha'] if 'aplha' in parameters else 1.0
        
#-----------------------------------  

def Bode(P,omega, Show = True):
    
    """
    :P: Process as defined by the class "Process".
        Use the following command to define the default process which is simply a unit gain process:
            P = Process({})
        
        A delay, two lead time constants and 2 lag constants can be added.
        
        Use the following commands for a SOPDT process:
            P.parameters['Kp'] = 1.1
            P.parameters['Tlag1'] = 10.0
            P.parameters['Tlag2'] = 2.0
            P.parameters['theta'] = 2.0
        
        Use the following commands for a unit gain Lead-lag process:
            P.parameters['Tlag1'] = 10.0        
            P.parameters['Tlead1'] = 15.0     

    :C: PID controller as defined by the class "Controller".
        Use the following command to define the default PID :
            P = Process({})
        
    :omega: frequency vector (rad/s); generated by a command of the type "omega = np.logspace(-2, 2, 10000)".
    :Show: boolean value (optional: default value = True). If Show = True, the Bode diagram is shown. Otherwise Ps (P(j omega)) (vector of complex numbers) is returned.
    
    The function "Bode" generates the Bode diagram of the process P
    """     
    
    s = 1j*omega
    
    Ptheta = np.exp(-P.parameters['theta']*s)
    PGain = P.parameters['Kp']*np.ones_like(Ptheta)
    PLag1 = 1/(P.parameters['Tlag1']*s + 1)
    PLag2 = 1/(P.parameters['Tlag2']*s + 1)
    PLead1 = P.parameters['Tlead1']*s + 1
    PLead2 = P.parameters['Tlead2']*s + 1
    
    Ps = np.multiply(Ptheta,PGain)
    Ps = np.multiply(Ps,PLag1)
    Ps = np.multiply(Ps,PLag2)
    Ps = np.multiply(Ps,PLead1)
    Ps = np.multiply(Ps,PLead2)
    
    if Show == True:
    
        fig, (ax_gain, ax_phase) = plt.subplots(2,1)
        fig.set_figheight(12)
        fig.set_figwidth(22)

        # Gain part
        ax_gain.semilogx(omega,20*np.log10(np.abs(Ps)),label='P(s)')
        ax_gain.semilogx(omega,20*np.log10(np.abs(PGain)),label='Pgain')
        if P.parameters['theta'] > 0:
            ax_gain.semilogx(omega,20*np.log10(np.abs(Ptheta)),label='Ptheta(s)')
        if P.parameters['Tlag1'] > 0:
            ax_gain.semilogx(omega,20*np.log10(np.abs(PLag1)),label='PLag1(s)')
        if P.parameters['Tlag2'] > 0:        
            ax_gain.semilogx(omega,20*np.log10(np.abs(PLag2)),label='PLag2(s)')
        if P.parameters['Tlead1'] > 0:        
            ax_gain.semilogx(omega,20*np.log10(np.abs(PLead1)),label='PLead1(s)')
        if P.parameters['Tlead2'] > 0:    
            ax_gain.semilogx(omega,20*np.log10(np.abs(PLead2)),label='PLead2(s)')    
        gain_min = np.min(20*np.log10(np.abs(Ps)/5))
        gain_max = np.max(20*np.log10(np.abs(Ps)*5))
        ax_gain.set_xlim([np.min(omega), np.max(omega)])
        ax_gain.set_ylim([gain_min, gain_max])
        ax_gain.set_ylabel('Amplitude |P| [db]')
        ax_gain.set_title('Bode plot of P')
        ax_gain.legend(loc='best')
    
        # Phase part
        ax_phase.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(Ps)),label='P(s)')
        ax_phase.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(PGain)),label='Pgain')
        if P.parameters['theta'] > 0:    
            ax_phase.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(Ptheta)),label='Ptheta(s)')
        if P.parameters['Tlag1'] > 0:        
            ax_phase.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(PLag1)),label='PLag1(s)')
        if P.parameters['Tlag2'] > 0:        
            ax_phase.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(PLag2)),label='PLag2(s)')
        if P.parameters['Tlead1'] > 0:        
            ax_phase.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(PLead1)),label='PLead1(s)')
        if P.parameters['Tlead2'] > 0:        
            ax_phase.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(PLead2)),label='PLead2(s)')    
        ax_phase.set_xlim([np.min(omega), np.max(omega)])
        ph_min = np.min((180/np.pi)*np.unwrap(np.angle(Ps))) - 10
        ph_max = np.max((180/np.pi)*np.unwrap(np.angle(Ps))) + 10
        ax_phase.set_ylim([np.max([ph_min, -200]), ph_max])
        ax_phase.set_ylabel(r'Phase $\angle P$ [°]')
        ax_phase.legend(loc='best')
    else:
        return Ps

#----------------------------------- 
def PID_RT(SP, PV, Man, MVMan, MVFF, Kc, Ti, Td, alpha, Ts, MVMin, MVMax, MV, MVP, MVI, MVD, E, ManFF=False, PVInit=0, method='EBD-EBD'):

    """
    Help on function PID_RT in module package_DBR_Advanced:
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
    
    #initialisation of E
    if len(PV) == 0 :
        E.append(SP[-1] - PVInit)
    else:
        E.append(SP[-1] - PV[-1])

    if Man[-1] == False:
        
        #Proportional
        MVP.append(Kc*E[-1])

        #Integral
        if len(MVI)==0:
            MVI.append((Kc*Ts/Ti)*E[-1])
        else:
            MVI.append(MVI[-1]+(Kc*Ts/Ti)*E[-1])

        #Derivating
        Tfd = alpha*Td

        if len(MVD) == 0 and len(E)>1:
            MVD.append((Kc*Td/(Tfd+Ts))*(E[-1]-E[-2]))
        elif len(E) == 1:
            MVD.append(0)
        else:
            MVD.append((Tfd/(Tfd+Ts))*MVD[-1]+(Kc*Td/(Tfd+Ts))*(E[-1]-E[-2]))

        MV.append(MVP[-1]+MVI[-1]+MVD[-1])

    else:


        """
        MVI.append(0)
        MVP.append(0)
        MVD.append(0)
        if len(MVMan)==0:
            MV.append(0)
        else:
            MV.append(MVMan[-1])
        """

#----------------------------------- 

def StabilityMargins(P: Process, C: PID, omega):
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