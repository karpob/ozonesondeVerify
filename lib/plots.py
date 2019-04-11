import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np

def plotSondeAndAnalysisStats(p, o, idx,top, bottom, controlName, experimentName,  tag):

    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col',sharey='row')
    f.suptitle('Ozone Sonde Count:{:d}'.format(int(max(o['count_sonde']) )) ) 
    ax1.plot( o['av_sonde'][idx], p[idx],'k', label='Sondes' )
    ax1.plot( o['av_ana1'][idx], p[idx],'g', label = controlName )
    ax1.set_yscale('log')
    ax1.set_ylim([bottom,top])
    ax1.invert_yaxis()
    ax1.set_yticks(np.array([1000.0, 100.0, 10.0,1.0]))
    ax1.set_yticklabels(['1000.0','100.0','10.0','1.0'])
    ax1.legend(prop={'size':7})

    ax3.plot( o['av_sonde'][idx], p[idx],'k', label='Sondes' )
    ax3.plot( o['av_ana2'][idx], p[idx],'g', label = experimentName )
    ax3.set_yscale('log')
    ax3.set_ylim([bottom,top])
    ax3.invert_yaxis()
    ax3.set_yticks(np.array([1000.0, 100.0, 10.0,1.0]))
    ax3.set_yticklabels(['1000.0','100.0','10.0','1.0'])
    ax3.legend(prop={'size':7})

    d = o['av_ana1'][idx] - o['av_sonde'][idx]
    d1max = max(abs(d)) + max(o['std_ana1'][idx])
    ax2.plot( d,p[idx],'k')
    ax2.plot( d,p[idx],'ko',markersize=2.5)
    ax2.set_ylim([bottom,top])
    ax2.fill_betweenx(p[idx],d+o['std_ana1'][idx], d-o['std_ana1'][idx], color='grey' )
    ax2.axvline(x=0, color='black')

    d = o['av_ana2'][idx] - o['av_sonde'][idx]
    dmax = max(abs(d)) 
    dmax = max(abs(d)) + max(o['std_ana2'][idx])
    if(d1max>dmax):dmax = d1max

    ax4.plot( d,p[idx],'k')
    ax4.plot( d,p[idx],'ko', markersize=2.5)
    ax4.fill_betweenx(p[idx],d+o['std_ana2'][idx], d-o['std_ana2'][idx], color='grey' )
    ax4.set_ylim([bottom,top])
    ax4.set_xlim((-1.0*dmax,dmax))
    ax4.axvline(x=0, color='black')
 
    ax1.set_ylabel('Pressure [hPa]')
    ax3.set_ylabel('Pressure [hPa]')
    ax3.set_xlabel('Ozone [mPa]')
    ax4.set_xlabel('Difference Ozone [mPa]')
    print("Saving plot: {}".format(tag+'_profileStats.png'))
    plt.savefig(tag+'_profileStats.png')
    plt.close()

