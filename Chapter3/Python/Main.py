"""
# Main code for Sufficient Statistic project
"""
import sys 
sys.path.insert(0,'OneAssetBond/')

import numpy as np
import defineSSParametersIOUsBond as Params
from copy import copy
import pickle
from SteadyState_fsolve import SteadyState_fsolve
from FluctuationsOneAssetIOUsBond import FluctuationsOneAssetIOUsBond, Fsys, SGU_solver, plot_IRF
import matplotlib.pyplot as plt

EconomyParams = copy(Params.parm_one_asset_IOUsBond)
  

SSEconomy = SteadyState_fsolve(**EconomyParams)


##### Choose whether calculate new steady state or use ond one

SSS = SSEconomy.SolveSteadyState() # New steady state
pickle.dump(SSS, open("btoy5_tau_1.p", "wb"))

SSS=pickle.load(open("btoy5_tau_1.p", "rb")) # Use old steady state

#############3#################################################3



##### Preparation: Calculation of sufficient statistics

#   MPCs
par = SSS['par']
mpar = SSS['mpar']
targets = SSS['targets']
grid = SSS['grid']
c_guess = SSS['c_policy']
inc = SSS['inc']
jd = SSS['joint_distr'].copy()
joint_distr = jd.copy()
NW=(targets['N']/par['H'])*targets['W']
WW=NW*np.ones((mpar['nm'],mpar['nh'])) #Wages
WW[:,-1]=SSS['Profits_fc']*par['profitshare']

meshes_m,meshes_h = np.meshgrid(grid['m'],grid['h'], indexing='ij')
aux_x = par['tau']*targets['N']/par['H']*targets['W']*meshes_h/(1+par['gamma'])
aux_x[:,-1]=0
C_ind = c_guess + aux_x
Output = SSS['Output']
C_agg = SSS['C_agg']
B = SSS['targets']['B']
RB = SSS['par']['RB']

# MPC
WW_h=WW[0,:]
WW_h_mesh=np.multiply(WW,meshes_h)
grid_h_aux=grid['h']
MPC_m = np.zeros((mpar['nm'],mpar['nh']))

for hh in range(0,mpar['nh']):
#    MPC_m[:,hh]=np.gradient(np.squeeze(C_ind[:,hh]))/np.gradient(grid['m'])  # MPC_m_ is same with MPC_m
    MPC_m[:,hh]=np.gradient(np.squeeze(c_guess[:,hh]))/np.gradient(grid['m'])  # MPC_m_ is same with MPC_m
MPC_h = np.zeros((mpar['nm'],mpar['nh']))

for mm in range(0,mpar['nm']):
   MPC_h[mm,:]=np.gradient(C_ind[mm,:])/np.gradient(np.multiply(WW_h,grid_h_aux))

NNP = meshes_m
URE = inc['labor'] + meshes_m - c_guess + inc['profits']

# MPC_m=np.min(MPC_m,1)

##### Calculation of sufficient statistics

# 1. Aggregate income channel: EI[Yi/Y MPC_m]
Inc_wt_MPC = np.sum(np.sum( np.multiply(np.multiply(WW_h_mesh/np.sum(np.sum(np.multiply(jd.copy(),WW_h_mesh))), MPC_m), jd.copy()) ))

# 3. Fisher channel
Redist_elas_P = np.sum(np.sum(np.multiply(np.multiply(MPC_m.copy(), NNP.copy()), jd.copy()))) - np.sum(np.sum(np.multiply(MPC_m.copy(),jd.copy())))*np.sum(np.sum(np.multiply(NNP.copy(),jd.copy())))

# 4. Interest rate exposure channel
Redist_elas_R = np.sum(np.sum(np.multiply(np.multiply(MPC_m.copy(),URE.copy()),jd.copy()))) - np.sum(np.sum(np.multiply(MPC_m.copy(),jd.copy())))*np.sum(np.sum(np.multiply(URE.copy(),jd.copy())))

# 5. substitution channel
#sig_i = par['xi']**(-1)*np.multiply(c_guess,1/C_ind)
#Hick_scaling = np.sum(np.sum( np.multiply(np.multiply(np.multiply(sig_i, (1-MPC_m.copy())), C_ind.copy()), jd.copy()) ))

sig_i = par['xi']**(-1)
Hick_scaling = np.sum(np.sum( np.multiply(np.multiply(np.multiply(sig_i, (1-MPC_m.copy())), c_guess.copy()), jd.copy()) ))


# 6. additional term for GHH preference



#### changes of policy parameters: irrelevant of a steady state

#SSS['par']['rho_B'] = 0.99
#SSS['par']['gamma_pi'] = 1.25
SSS['par']['theta_pi'] = 2.0
##############################################################################



EX2SR=FluctuationsOneAssetIOUsBond(**SSS)

SR=EX2SR.StateReduc()

State       = np.zeros((SR['mpar']['numstates'],1))
State_m     = State
Contr       = np.zeros((SR['mpar']['numcontrols'],1))
Contr_m     = Contr

Fsysresult = Fsys(State, State_m, Contr, Contr_m, SR['Xss'], SR['Yss'], 
         SR['Gamma_state'], SR['Gamma_control'], SR['InvGamma'], SR['Copula'],SR['par'], 
         SR['mpar'], SR['grid'], SR['targets'],SR['P_H'],SR['aggrshock'],SR['oc'])
    
    
SGUresult=SGU_solver(SR['Xss'],SR['Yss'],SR['Gamma_state'],SR['Gamma_control'],SR['InvGamma'],SR['Copula'],
                     SR['par'],SR['mpar'],SR['grid'],SR['targets'],SR['P_H'],SR['aggrshock'],SR['oc'])


plotresult = plot_IRF(SR['mpar'],SR['par'],SGUresult['gx'],SGUresult['hx'],SR['joint_distr'],
         SR['Gamma_state'],SR['grid'],SR['targets'],SR['os'],SR['oc'],SR['Output'],C_ind,c_guess,
         WW_h_mesh,MPC_m,Inc_wt_MPC,Redist_elas_P,Redist_elas_R,Hick_scaling,URE,NNP)

os = SR['os']
oc = SR['oc']
IRF_state_sparse = plotresult['IRF_state_sparse']
IRF_RB = 100*IRF_state_sparse[mpar['numstates']-os,1:]
IRF_Y=100*IRF_state_sparse[-1-oc+2, :-1]
IRF_W=100*IRF_state_sparse[-1-oc+3, :-1]
PI=1+IRF_state_sparse[-1-oc+1, 1:]
RB=par['RB']*(1+IRF_state_sparse[mpar['numstates']-os,1:])
IRF_RBREAL=100*100*(RB/PI-par['RB'])
IRF_RB=100*100*(RB-par['RB'])
IRF_PI=100*100*IRF_state_sparse[-1-oc+1, :-1]
IRF_Cagg = 100*IRF_state_sparse[-1-oc+8, :-1]
IRF_Xagg = 100*IRF_state_sparse[-1-oc+9, :-1]
    
print(round(plotresult['IRF_Xagg'][0,0],3))
print(round(plotresult['IRF_X_by_suff'][0,0],3))

factor1percent = -100.0/(plotresult['IRF_RB'][0,0])

plt.figure()
line1,=plt.plot(range(1,mpar['maxlag']),factor1percent*np.squeeze(np.asarray(plotresult['IRF_Cagg'])),label='Consumption IRF')
plt.legend(handles=[line1])
plt.plot(range(0,mpar['maxlag']-1),np.zeros((mpar['maxlag']-1)),'k--' )
plt.xlabel('Quarter')
plt.ylabel('Percent') 
plt.title('Impulse Response to Consumption')
plt.savefig('./Figures/IRF_C.pdf',format='pdf')

plt.figure()
line1,=plt.plot(range(1,mpar['maxlag']),factor1percent/100*np.squeeze(np.asarray(plotresult['IRF_RB'])),label='Real Interest Rate IRF')
plt.plot(range(0,mpar['maxlag']-1),np.zeros((mpar['maxlag']-1)),'k--' )
plt.legend(handles=[line1])
plt.xlabel('Quarter')
plt.ylabel('Percent') 
plt.title('Impulse Response for Real Interest Rate')
plt.savefig('./Figures/IRF_RB.pdf',format='pdf')

plt.figure()
line1,=plt.plot(range(1,mpar['maxlag']),factor1percent*np.squeeze(np.asarray(plotresult['IRF_W'])),label='Real Wages IRF')
plt.plot(range(0,mpar['maxlag']-1),np.zeros((mpar['maxlag']-1)),'k--' )
plt.legend(handles=[line1])
plt.xlabel('Quarter')
plt.ylabel('Percent') 
plt.title('Impulse Response for Real Wages')
plt.savefig('./Figures/IRF_W.pdf',format='pdf')

plt.figure()
line1,=plt.plot(range(1,mpar['maxlag']),factor1percent/100*np.squeeze(np.asarray(plotresult['IRF_PI'])),label='Inflation IRF')
plt.plot(range(0,mpar['maxlag']-1),np.zeros((mpar['maxlag']-1)),'k--' )
plt.legend(handles=[line1])
plt.xlabel('Quarter')
plt.ylabel('Percent') 
plt.title('Impulse Response for Inflation')
plt.savefig('./Figures/IRF_PI.pdf',format='pdf')

decom_table = "    \\begin{table}  \n"
decom_table += '\\begin{center} \n'
decom_table +=  '\\caption{Transmission Channel Importance}\label{table:trans_channel}'

decom_table += "\\begin{tabular}{lc}  \n"
decom_table += "\\toprule \n"

decom_table += "Aggregate Income & {:.1f}".format(100*plotresult['IRF_AggInc'][0,0]/plotresult['IRF_Cagg'][0,0]) +"\\% \n"
decom_table += "\\\\ Earnings Heterogeneity & {:.1f}".format(100*plotresult['IRF_EarnHet'][0,0]/plotresult['IRF_Cagg'][0,0]) +"\\% \n"
decom_table += "\\\\ Interest Rate Exposure & {:.1f}".format(100*plotresult['IRF_IRE'][0,0]/plotresult['IRF_Cagg'][0,0]) +"\\% \n"
decom_table += "\\\\ Fisher & {:.1f}".format(100*plotresult['IRF_Fisher'][0,0]/plotresult['IRF_Cagg'][0,0]) +"\\% \n"
decom_table += "\\\\ Intertemporal Substitution & {:.1f}".format(100*plotresult['IRF_IntSubs'][0,0]/plotresult['IRF_Cagg'][0,0]) +"\\% \n"
decom_table += "\\\\ GHH Channel & {:.1f}".format(100*plotresult['IRF_GHH'][0,0]/plotresult['IRF_Cagg'][0,0]) +"\\% \n"
decom_table += "\\\\ Error & {:.1f}".format(100*(plotresult['IRF_Cagg'][0,0]-plotresult['IRF_C_by_suff'][0,0])/plotresult['IRF_Cagg'][0,0]) +"\\% \n"

decom_table += "\\\\ \\bottomrule  \n"
decom_table += " \end{tabular}   \n"
decom_table += '\\end{center} \n'
decom_table += "\\end{table}  \n"

with open('./Tables/Decomposition.tex','w') as f:
    f.write(decom_table)
    f.close()

